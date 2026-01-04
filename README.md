# gem5-Cache-Thrashing-Mitigation-and-Software-Defined-QoS-Victim-Filter-Cache-research

**Author:** Fredrick Chiu  
**Affiliation:** Rice University ECE  
**Context:** Advanced Computer Architecture Research (Projects A, B, C)

## Overview
This repository (`gem5-Cache-Thrashing-QVC-research`) is a three-part investigation of cache thrashing in gem5, progressing from *measurement* to *stress* to *mitigation*.

It is organized as:
1. **Project A (Single-core self-contention):** Establishes how sensitive a pointer-chasing graph workload is to cache/memory behavior on a single core, and shows the need of protecting the important data.
2. **Project B (Multicore contention):** Demonstrates the scaling problem once multiple cores share the cache hierarchy: shared cache capacity/bandwidth becomes a bottleneck under mixed workloads, and the victim’s effective share can collapse as core count rises.   
3. **Project C (Mitigation via QVC):** Implements and evaluates a QoS mechanism (QVC) as an intercepting filter that consults a software-managed range table; protected requests can be served locally while unprotected requests are forwarded to the baseline Ruby hierarchy. 

Across these phases, the core message is that (1) modern cores can partially mask latency, (2) shared-cache contention creates disproportionate slowdowns under workload mixing, and (3) a software-defined protection mechanism can claw back performance when critical pointer-chasing regions are kept close to the core. 

## Key findings
- Although Out-of-order execution can hide some miss latency, the performance degradations measured were still siginificant in both single-core and multi-core experiments. 
- In multicore scaling, shared-cache pressure and bandwidth division make the victim’s per-core cache service effectively shrink with core count, motivating per-application QoS rather than one-size-fits-all sharing  
- The QVC's filtering approach is simpler to prototype than a deep Ruby coherence-protocol integration, but it can incur simulation overhead and could break inherent optimization between the CPU and Ruby; despite that, the report shows a measured 28% gain for the protected server workload.

Please investigate the report to understand the detail.
---

## 🔬 Project Phases

### Phase 1: The Vulnerability (Project A & B)
We characterized a Dijkstra Single-Source Shortest Path (SSSP) kernel on a modeled RISC-V Out-of-Order (O3) CPU.
*   **Finding:** Victim workloads (e.g. pointer-chasing programs) have high locality and their data will be evicted from the shared L2 or L3. L3 misses are the primary bottleneck.
*   **The Aggressor:** We introduced a synthetic Vector Aggressor workload (SAXPY with RVV intrinsics). Note that in project A, the aggressor lies in the same program, modeling self-polluting programs.
*   **Result:** When sharing an L3 cache, the Vector Aggressor creates bursty, high-bandwidth traffic that evicts the workload's hot working set, demonstrating siginificant performance degradation (e.g. 90% IPC drop)

### Phase 2: The Solution (Project C)
We implemented a **Software-Defined QoS Victim Filter Cache (QVC)**, or "Victim Filter." This is a sidecar hardware unit that sits between the CPU Core and the L1 Cache. It talks to Ruby MESI Three-level Cache Hierarchy and the CPU.

#### Architecture: Software-Defined Filter
Instead of complex hardware replacement policies, we let the software define what matters.
1.  **Identify:** The application marks its "Hot Set" (e.g., the `dist` array and adjacency list headers) using a custom instruction.
2.  **Filter:** The QVC intercepts memory requests.
    *   **Protected Range?** Service locally (1-cycle latency).
    *   **Unprotected?** Forward immediately to the standard cache hierarchy.
3.  **Result:** Even when the L3 is thrashing, the QVC guarantees hits for the critical data, recovering the performacne (IPC).

#### Implementation Highlights
*   **New SimObject:** `src/mem/VictimL1/` contains the C++ logic for the filter cache.
*   **Byte-Granularity Validity:** Implemented bit-masks to handle sub-block fills from the Ruby memory system, preventing null-pointer faults.
*   **Protected Data Registration:** The pseudo-instruction (`qvcCtrl`) registers the protected data ranges.
*   **Virtual-to-Physical Translation:** The pseudo-instruction (`qvcCtrl`) includes a page-table walker to robustly register physical frames in fragmented memory.

---

## 📂 Repository Structure

### Core Gem5 Modifications
*   `src/mem/VictimL1/`: **[NEW]** The Victim Cache SimObject (C++ logic and Python parameters).
*   `src/sim/pseudo_inst.cc`: **[MODIFIED]** Added `qvcCtrl` logic with protected data registration and page table translation.
*   `include/gem5/m5ops.h`: **[MODIFIED]** Added opcode definition (`0x58`) for the custom instruction.

### Configuration & Scripts
*   `configs/Server_VictimL1_1C_config.py`: The master simulation script. Configures the RISC-V O3CPU, Ruby MESI_Three_Level protocol, and dynamically inserts the QVC for Core 0.

### Workloads
*   `workloads/workload_magic_wand_1A.cpp`: Project A's self-polluted victim workload. It models a common image editor scenario.
*   `workloads/barnes_hut_1A.cpp`: Project A's self-polluted victim workload. It runs a Barnes-Hut N-body simulation for interparticle or interstellar interactions. 
*   `workloads/Dijkstra_BFS_Server_1B.cpp`: Project B's victim workload. The BFS Point-of-Interest Server Victim workload. It first runs an dijkstra's algorithm to find the path and then uses BFS to serve user queries.
*   `workloads/Server_VictimL1.cpp`: Same as project B but now with protected data. Instrumented with `m5_qvc_ctrl` to register the adjacency list headers.
*   `workloads/saxpy_vector_aggressor.cpp`: The Aggressor workload using RISC-V Vector intrinsics.
---

## 🚀 How to Build & Run

### 1. Build the Simulator
Requires dependencies (scons, python3, etc.) compatible with gem5 v24+.
```bash
# Build the RISC-V Optimized binary
scons build/RISCV/gem5.opt -j$(nproc)
```
### 2. Workload Compilation
A RISC-V toolchain that supports RVV is necessary. We also need to link gem5 source library for compilation.
```bash
riscv64-unknown-elf-g++ -O3 -march=rv64gcv -static \
    -I$GEM5_PATH/include -DGEM5_M5OPS_H \
    -L$GEM5_PATH/util/m5/build/riscv/out \
    Server_VictimL1.cpp \
    -lm5 \
    -o $1
```
### 3. Run with the Python Configuration File
```bash
./build/RISCV/gem5.opt  $PATH_TO_CONFIG/Server_VictimL1_1C_config.py --mode=multicore-contention --cores=16
```

Please check out the Python Configuration files to understand the command line options
