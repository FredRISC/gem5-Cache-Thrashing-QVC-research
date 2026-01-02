// src/mem/VictimL1/VictimL1.cc
#include "mem/VictimL1/VictimL1.hh"
#include "base/trace.hh"
#include "debug/VictimL1.hh"
#include <iomanip>

namespace gem5 {
namespace memory {

// for stats file
VictimL1::VictimL1Stats::VictimL1Stats(statistics::Group *parent)
    : statistics::Group(parent),
      hits(this, "hits", "Number of read hits in VictimL1"),
      misses(this, "misses", "Number of read misses forwarded to memory"),
      writeUpdates(this, "write_updates", "Number of write updates performed"),
      forwardReqs(this, "forwards", "Total requests forwarded to memory")
{
}


// =========================================================================
// Global Instance & Boilerplate
// =========================================================================


// Define a file-local pointer and implement the public accessor function
static VictimL1* g_victim_l1_ptr = nullptr;
VictimL1* VictimL1::getInstance() { return g_victim_l1_ptr; }

//Constructor of the L1 Victim Cache
VictimL1::VictimL1(const VictimL1Params &params) :
    ClockedObject(params),
    cpuPort(params.name + ".cpu_side", this),
    memPort(params.name + ".mem_side", this),
    size(params.size),
    assoc(params.assoc),
    blockSize(params.block_size),
    pendingResp(nullptr),
    sendEvent([this]{ processSendEvent(); }, name()),
    cpuSideBlocked(false),
    stats(this)
{
    // Initialize cache storage
    numSets = size / (assoc * blockSize);
    sets.resize(numSets);
    for (auto &set : sets) {
        for (int i = 0; i < assoc; ++i) {
            set.emplace_back(blockSize);
        }
    }
		g_victim_l1_ptr = this;
}


Port &VictimL1::getPort(const std::string &if_name, PortID idx)
{
    if (if_name == "cpu_side") return cpuPort;
    if (if_name == "mem_side") return memPort;
    return ClockedObject::getPort(if_name, idx);
}


// =========================================================================
// Software Control Logic
// =========================================================================
void VictimL1::addProtectedRange(Addr start, uint64_t size) {
    DPRINTF(VictimL1, "Adding protected range: %#x - %#x\n", start, start + size);
    //std::cout << "DEBUG_RANGE: " << std::hex << start << " - " << (start + size) << std::dec << std::endl;
    protectedRanges.emplace_back(start, start + size);
}

void VictimL1::clearProtectedRanges() {
    DPRINTF(VictimL1, "Clearing all protected ranges.\n");
    protectedRanges.clear();
}

bool VictimL1::isProtected(Addr addr) {
    for (const auto &r : protectedRanges) {
        if (r.contains(addr)) return true;
    }
    return false;
}


// =========================================================================
// Cache Management Logic
// =========================================================================

// Helper: Just find the line by tag. Used for Writes and Allocations.
VictimL1::CacheLine* VictimL1::findLine(Addr addr) {
    Addr blkAddr = addr & ~(blockSize - 1);
    unsigned setIdx = (blkAddr / blockSize) % numSets;
    Addr tag = blkAddr / (blockSize * numSets);

    for (auto &line : sets[setIdx]) {
        if (line.valid && line.tag == tag) { // check at cacheline granularity
            return &line;
        }
    }
    return nullptr;
}

// Helper: Find line by tag AND verify data validity. Used for READS.
VictimL1::CacheLine* VictimL1::accessLine(Addr addr, unsigned size) {
    CacheLine* line = findLine(addr);
    if (line) {
        // We found the tag, but do we have the specific bytes?
        unsigned offset = addr & (blockSize - 1);
        if (line->isRangeValid(offset, size)) {
            line->lastTouch = curTick();
            return line;
        }
    }
    return nullptr;
}

// Helper: Allocate or Reuse a line. Used on Fill.
VictimL1::CacheLine* VictimL1::allocateLine(Addr addr) {
    // 1. Check if we already have the line (Merge case)
    CacheLine* existing = findLine(addr);
    if (existing) {
        existing->lastTouch = curTick();
        return existing;
    }

    // 2. Need new slot. Calculate Set.
    Addr blkAddr = addr & ~(blockSize - 1);
    unsigned setIdx = (blkAddr / blockSize) % numSets;
    Addr tag = blkAddr / (blockSize * numSets);
    
    auto &set = sets[setIdx];
    CacheLine* victim = &set[0];

    // 3. Search for invalid or LRU
    for (auto &line : set) {
        if (!line.valid) {
            victim = &line;
            break;
        }
        if (line.lastTouch < victim->lastTouch) {
            victim = &line;
        }
    }

    victim->valid = true;
    victim->tag = tag;
    victim->lastTouch = curTick();
    // CRITICAL: Reset the byte mask because this is a new allocation (or eviction reuse)
    std::fill(victim->byteValid.begin(), victim->byteValid.end(), false);
    return victim;
}



// =========================================================================
// Port Logic
// =========================================================================

// Track packets that missed in QVC so we fill on return
// We use a sender state to tag these packets
struct VictimL1SenderState : public Packet::SenderState {
    bool allocateOnReturn;
    VictimL1SenderState(bool alloc) : allocateOnReturn(alloc) {}
};


// Request Logic (Write-Through)
bool VictimL1::CPUSidePort::recvTimingReq(PacketPtr pkt) {
    bool protectedAddr = owner->isProtected(pkt->getAddr());

    // 1. READ HIT
    // Must be Protected AND Have valid data bytes
    if (pkt->isRead() && protectedAddr) {
        auto line = owner->accessLine(pkt->getAddr(), pkt->getSize());
        if (line) {
            // Check flow control
            if (owner->pendingResp != nullptr) {
                owner->cpuSideBlocked = true;
                return false; 
            }

            // INCREMENT HIT
            owner->stats.hits++;

						// Response logic
            pkt->makeResponse();
            pkt->setData(line->data.data() + pkt->getOffset(owner->blockSize));
            
            // Debug print
            uint64_t *debug_val = (uint64_t*)(line->data.data() + pkt->getOffset(owner->blockSize));
            DPRINTF(VictimL1, "READ HIT Addr: %#x | Val: %#016x\n", pkt->getAddr(), *debug_val);
            
            // Schedule response next cycle
            owner->pendingResp = pkt;
            owner->schedule(owner->sendEvent, owner->clockEdge(Cycles(1)));
            return true;
        }
        else {
            // INCREMENT MISS (Protected but not found/valid)
            owner->stats.misses++;
        }
    }

    // 2. WRITE UPDATE
    // If we have the line (even partial), update it to keep it fresh.
    if (pkt->isWrite() && protectedAddr) {
        auto line = owner->findLine(pkt->getAddr());
        if (line) {
            // Write logic
            pkt->writeData(line->data.data() + pkt->getOffset(owner->blockSize));
            // Mark these bytes as valid!
            line->validateRange(pkt->getOffset(owner->blockSize), pkt->getSize());
            line->lastTouch = curTick();
            
            // INCREMENT WRITE UPDATE
            owner->stats.writeUpdates++;
        }
    }

    // 3. FORWARDING
    // INCREMENT FORWARD
    owner->stats.forwardReqs++;
    
    bool allocOnReturn = (pkt->isRead() && protectedAddr);
    pkt->pushSenderState(new VictimL1SenderState(allocOnReturn));
    
    if (owner->memPort.sendTimingReq(pkt)) {
        return true;
    } else {
        Packet::SenderState *popped = pkt->popSenderState();
        delete popped;
        return false;
    }
}

// --- Response Logic ---
bool VictimL1::MemSidePort::recvTimingResp(PacketPtr pkt) {
    VictimL1SenderState *state = dynamic_cast<VictimL1SenderState*>(pkt->popSenderState());
    
    if (state && state->allocateOnReturn) {
        // Allocate (or find existing) line
        auto line = owner->allocateLine(pkt->getAddr());
        
        // Write data
        pkt->writeData(line->data.data() + pkt->getOffset(owner->blockSize));
        
        // Mark bytes valid
        line->validateRange(pkt->getOffset(owner->blockSize), pkt->getSize());
        
        DPRINTF(VictimL1, "FILL Addr: %#x | Validated %d bytes\n", pkt->getAddr(), pkt->getSize());
    }
    
    delete state;
    
    // Forward response to CPU (handling flow control)
    if (owner->pendingResp) {
        // Collision: pendingResp is waiting, but now MemSide brings another response.
        // In simple CPU models, we can send it immediately.
        // But strictly we should return false if we can't send.
        // O3CPU usually accepts responses unless full.
        if (!owner->cpuPort.sendTimingResp(pkt)) {
            panic("VictimL1: Failed to forward MemSide response! Queue full logic required.");
        }
        return true;
    }
    
    return owner->cpuPort.sendTimingResp(pkt);
}

void VictimL1::processSendEvent() {
    if (pendingResp) {
        if (cpuPort.sendTimingResp(pendingResp)) {
            pendingResp = nullptr;
            if (cpuSideBlocked) {
                cpuSideBlocked = false;
                cpuPort.sendRetryReq();
            }
        }
    }
}


// If we tried to send a Response to the CPU and it was blocked, the CPU calls this later to say "I'm ready now." We retry sending.
void VictimL1::CPUSidePort::recvRespRetry() {
    // 1. Retry Pending
    if (owner->pendingResp) {
        if (owner->cpuPort.sendTimingResp(owner->pendingResp)) {
            owner->pendingResp = nullptr;
            if (owner->cpuSideBlocked) {
                owner->cpuSideBlocked = false;
                owner->cpuPort.sendRetryReq();
            }
        }
    }
    // 2. Retry Memory Side
    owner->memPort.sendRetryResp();
}

// If we tried to send a Request to Ruby and it was blocked, Ruby calls this later. We retry sending.
void VictimL1::MemSidePort::recvReqRetry() { 
	owner->cpuPort.sendRetryReq(); 
}


// Boilerplate mandatory overrides

// Used for debugging/syscalls (instant data access, zero time). We just forward it to Ruby so the debugger sees coherent memory.
void VictimL1::CPUSidePort::recvFunctional(PacketPtr pkt) { owner->memPort.sendFunctional(pkt); }
// What addresses do we handle? We ask Ruby (memPort) and pass that answer up to the CPU.
AddrRangeList VictimL1::CPUSidePort::getAddrRanges() const { return owner->memPort.getAddrRanges(); }
// Just forward to memory. We don't cache in atomic mode for simplicity.
Tick VictimL1::CPUSidePort::recvAtomic(PacketPtr pkt){ return owner->memPort.sendAtomic(pkt);}

} // namespace memory
} // namespace gem5

