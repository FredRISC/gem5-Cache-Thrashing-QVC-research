// src/mem/VictimL1/VictimL1.hh
#ifndef __MEM_QVC_VICTIM_L1_HH__
#define __MEM_QVC_VICTIM_L1_HH__

#include "mem/port.hh"
#include "mem/packet.hh"
#include "sim/clocked_object.hh"
#include "params/VictimL1.hh"
#include <vector>
#include <list>
#include <unordered_map>
#include "sim/eventq.hh"
#include <deque>
#include "base/statistics.hh"

namespace gem5 {
namespace memory {

class VictimL1 : public ClockedObject {
  private:
  
  
    // --- Cache Storage Structures ---
    struct CacheLine {
        Addr tag;
        bool valid; // Is the tag valid (cacheline granularity)
        std::vector<uint8_t> data;
        std::vector<bool> byteValid; // Is the data byte valid in the cacheline (byte granularity)
        Tick lastTouch; // For LRU
        
        CacheLine(unsigned blkSize) : valid(false), tag(0), lastTouch(0) {
            data.resize(blkSize);
            byteValid.resize(blkSize, false);
        }
        
        // Check if a specific range of bytes is valid
        bool isRangeValid(unsigned offset, unsigned size) const {
            if (offset + size > byteValid.size()) return false;
            for (unsigned i = 0; i < size; ++i) {
                if (!byteValid[offset + i]) return false;
            }
            return true;
        }

        // Mark a range of bytes as valid
        void validateRange(unsigned offset, unsigned size) {
            if (offset + size > byteValid.size()) return;
            for (unsigned i = 0; i < size; ++i) {
                byteValid[offset + i] = true;
            }
        }
    };
  
  
    // --- Port Definitions ---
    
    //ResponsePort (CPUSide): Faces the CPU. It receives "recvTimingReq" Requests and sends "sendTimingResp" Responses.
    class CPUSidePort : public ResponsePort {
      private:
        VictimL1 *owner;
      public:
        CPUSidePort(const std::string& name, VictimL1 *owner) :
            ResponsePort(name), owner(owner) {}
        bool recvTimingReq(PacketPtr pkt) override;
        void recvFunctional(PacketPtr pkt) override;
        AddrRangeList getAddrRanges() const override;
        void recvRespRetry() override;
        Tick recvAtomic(PacketPtr pkt) override;
    };
		
		//RequestPort (MemSide): Faces the rest of the memory system (Ruby). It sends "sendTimingReq" Requests and receives "recvTimingResp" Responses.
    class MemSidePort : public RequestPort {
      private:
        VictimL1 *owner;
      public:
        MemSidePort(const std::string& name, VictimL1 *owner) :
            RequestPort(name), owner(owner) {}
        bool recvTimingResp(PacketPtr pkt) override;
        void recvReqRetry() override;
    };


    // --- Class Members ---
    CPUSidePort cpuPort;
    MemSidePort memPort;
    
    // Cache State
    unsigned size;
    unsigned assoc;
    unsigned blockSize;
    unsigned numSets;
    std::vector<std::vector<CacheLine>> sets; // Simple Set-Associative Cache: vector of sets
    
    // Software Control
		std::vector<AddrRange> protectedRanges; // Store a list of ranges
    
    // Flow Control & Scheduling
    PacketPtr pendingResp; // Buffer for 1 response
    EventFunctionWrapper sendEvent;
		bool cpuSideBlocked; //Did we tell CPU "false" recently?
    void processSendEvent();
        
    // --- Internal Helpers ---
    CacheLine* findLine(Addr addr); // Just finds tag (for Writes/Alloc)
    CacheLine* accessLine(Addr addr, unsigned size); // Finds tag AND checks data valid (for Reads)
    CacheLine* allocateLine(Addr addr); // Finds space or evicts

  public:
    VictimL1(const VictimL1Params &params);
    Port &getPort(const std::string &if_name, PortID idx=InvalidPortID) override;
    // --- The Software Control API---
    // These will be called by the Pseudo_Inst
    static VictimL1* getInstance();
    void addProtectedRange(Addr start, uint64_t size);
    void clearProtectedRanges();
    bool isProtected(Addr addr);
    
    
    // For stats file
    struct VictimL1Stats : public statistics::Group {
        statistics::Scalar hits;
        statistics::Scalar misses;
        statistics::Scalar writeUpdates;
        statistics::Scalar forwardReqs;
        
        VictimL1Stats(statistics::Group *parent);
    } stats;
    
};

} // namespace memory
} // namespace gem5

#endif // __MEM_QVC_VICTIM_L1_HH__

