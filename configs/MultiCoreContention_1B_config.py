# MultiCoreContention.py
# Project 1B: Multi-Core Spatial Cache Contention
# gem5 v24 component-based configuration for dual-core RISC-V system

import m5
from m5.objects import *
import argparse

# Import gem5 v24 components
from gem5.components.boards.simple_board import SimpleBoard
from gem5.components.processors.base_cpu_core import BaseCPUCore
from gem5.components.processors.base_cpu_processor import BaseCPUProcessor
from gem5.components.cachehierarchies.ruby.mesi_three_level_cache_hierarchy import (
    MESIThreeLevelCacheHierarchy,
)
from gem5.components.memory.single_channel import SingleChannelDDR3_1600, DIMM_DDR5_6400
from gem5.simulate.simulator import Simulator
from gem5.utils.requires import requires
from gem5.resources.resource import BinaryResource, obtain_resource
from gem5.isas import ISA

# Ensure RISC-V ISA
requires(isa_required=ISA.RISCV)

# Command-line arguments
parser = argparse.ArgumentParser(
    description="Project 1B: Multi-Core Cache Contention"
)

# Workload selection
parser.add_argument(
    "--victim",
    type=str,
    default="core0_150G",
    #required=True,
    help="Path to Dijkstra victim binary (Core 0)"
)
parser.add_argument(
    "--core1_workload",
    type=str,
    default="saxpy_vector_aggressor",
    help="Path to vector streaming aggressor binary (Core 1)"
)
parser.add_argument(
    "--scalar-aggressor",
    type=str,
    default="saxpy_scalar_aggressor",
    help="Path to scalar streaming aggressor binary (Core 1, for control)"
)
parser.add_argument(
    "--mode",
    type=str,
    choices=["baseline", "multicore-contention", "scalar-contention"],
    default="baseline",
    help="Experiment mode: baseline, multicore-contention, or scalar-contention"
)

# Architecture parameters
parser.add_argument(
    "--l1i-size",
    type=str,
    default="32kB",
    help="L1 instruction cache size (default: 32kB)"
)
parser.add_argument(
    "--l1d-size",
    type=str,
    default="32kB",
    help="L1 data cache size (default: 32kB)"
)
parser.add_argument(
    "--l2-size",
    type=str,
    default="256kB",
    help="private L2 cache size (default: 256kB)"
)
parser.add_argument(
    "--l3-size",
    type=str,
    default="2MB",
    help="Shared L3 cache size (default: 2MB)"
)
parser.add_argument(
    "--vlen",
    type=int,
    default=128,
    help="Vector register length in bits (default: 128)"
)
parser.add_argument(
    "--elen",
    type=int,
    default=64,
    help="Max element width in bits (default: 64)"
)

parser.add_argument(
    "--cores",
    type=int,
    default=1,
    help="Number of Cores"
)

args = parser.parse_args()


print(f"Experiment Mode: {args.mode}")


# --- Custom RVV Core (gem5 v24 style) ---
class RVVCore(BaseCPUCore):
    def __init__(self, elen, vlen, cpu_id):
        # Initialize with RiscvO3CPU and RISC-V ISA
        super().__init__(core=RiscvO3CPU(cpu_id=cpu_id), isa=ISA.RISCV)
        
        # Configure RISC-V Vector Extension
        self.core.isa[0].elen = elen
        self.core.isa[0].vlen = vlen
        
        # Enlarge LSQ for vector operations
        self.core.LQEntries = 256
        self.core.SQEntries = 256

# --- System Components ---

# Processor: Cores with RVV support
# Each core will be assigned its own process/binary
if args.mode == "multicore-contention":
    cores = [RVVCore(args.elen, args.vlen, cpu_id=i) for i in range(args.cores)]
    processor = BaseCPUProcessor(cores=cores)
else:
    processor = BaseCPUProcessor(
        cores=[RVVCore(args.elen, args.vlen, cpu_id=0)]
    )    


# Memory: Single-channel DDR3-1600 OR Dual-Channel DIMM_DDR5_6400
memory = DIMM_DDR5_6400(size="4GiB")


# Cache Hierarchy: MESI coherence with shared L2
# The shared L2 is where spatial contention occurs!
cache_hierarchy = MESIThreeLevelCacheHierarchy(
    l1d_size=args.l1d_size,
    l1d_assoc=8,
    l1i_size=args.l1i_size,
    l1i_assoc=8,
    l2_size=args.l2_size,
    l2_assoc=16,
    l3_size=args.l3_size,
    l3_assoc=128, 
    num_l3_banks=1   
)


# Board: Assemble components
board = SimpleBoard(
    clk_freq="3GHz",
    processor=processor,
    memory=memory,
    cache_hierarchy=cache_hierarchy,
)

# --- Workload Assignment ---

board.set_workload(obtain_resource("riscv-gapbs-bfs-run"))


aggressor_binary_path_list = []
aggressor_process_list = []
PID = 200

#core0 workload setup
victim_binary = BinaryResource(local_path="./configs/project/" + args.victim)
victim_binary_path = victim_binary.get_local_path()
victim_process = Process()
victim_process.executable = victim_binary_path
victim_process.cmd = [victim_binary_path] #+ arguments
victim_process.pid = PID


#victim_binary = BinaryResource(local_path="./configs/project/" + args.victim)


if args.mode == "multicore-contention":
    # Vector Contention: Core 0 victim + Core 1 memory-intensive workload + Core N random
    if not args.core1_workload:
        raise ValueError("--core1_workload required for multicore-contention mode")
    #core1 runs custom process
    for i in range(1,args.cores):
        aggressor_binary = BinaryResource(local_path="./configs/project/" + args.core1_workload)# if (i == 1) else obtain_resource("riscv-gapbs-bfs-run")
        aggressor_binary_path = aggressor_binary.get_local_path()
        aggressor_process = Process()
        aggressor_process.executable = aggressor_binary_path
        aggressor_process.cmd = [aggressor_binary_path] #+ arguments

        PID = PID + 1
        aggressor_process.pid = PID
        aggressor_process_list.append(aggressor_process)     
        aggressor_binary_path_list.append(aggressor_binary_path)    
   

elif args.mode == "scalar-contention":
    # Scalar Contention: Core 0 victim + Core 1 scalar aggressor (control)
    if not args.scalar_aggressor:
        raise ValueError("--scalar-aggressor required for scalar-contention mode")
    
    #aggressor_binary = BinaryResource(local_path="./configs/project/" + args.scalar_aggressor)
     
    aggressor_process = Process()
    aggressor_binary_path = "./configs/project/" + args.scalar_aggressor
    aggressor_process.executable = aggressor_binary_path
    aggressor_process.cmd = [aggressor_binary_path] #+ arguments
    PID = PID + 1
    aggressor_process.pid = PID



board.set_is_workload_set(True)
board.workload=SEWorkload.init_compatible(victim_binary_path)
board.processor.cores[0].set_workload(victim_process)
if args.mode == "multicore-contention":
    for i in range(1,args.cores):
        board.workload=SEWorkload.init_compatible(aggressor_binary_path_list[i-1])
        board.processor.cores[i].set_workload(aggressor_process_list[i-1])

#board.set_se_binary_workload(victim_binary)
#processor.cores[0].set_workload(process_list[0])
#processor.cores[1].set_workload(process_list[1])

board.memory.mem_ctrl[0].dram.range = AddrRange('0', '2GiB')
board.memory.mem_ctrl[1].dram.range = AddrRange('2GiB', '4GiB')


# --- Print System Configuration ---
print("")
print("System Configuration")
print("-" * 70)
print(f"  CPU Model: RiscvO3CPU (Out-of-Order)")
print(f"  Number of Cores: {args.cores}")
print(f"  Clock Frequency: 3GHz")
print(f"  Memory: DDR5-6400, 4GiB")
print("")
print("Cache Hierarchy (MESI Three-Level):")
print(f"  L1I: {args.l1i_size}, 8-way (private per core)")
print(f"  L1D: {args.l1d_size}, 8-way (private per core)")
print(f"  L2: {args.l2_size}, 16-way (private per core))")
print(f"  L3: {args.l3_size}, 128-way (shared among cores))")

print(f"  Coherence: MESI protocol via Ruby")
print("")
print("RISC-V Vector Extension:")
print(f"  VLEN: {args.vlen} bits")
print(f"  ELEN: {args.elen} bits")
print(f"  vl setting: vlmax (for maximum memory pressure)")
print("")

# --- Run Simulation ---

simulator = Simulator(board=board, full_system=False)
simulator.run()

print("")
print("=" * 70)
print("Simulation Complete!")
print(f"Simulated ticks: {simulator.get_current_tick()}")
print(f"Exit cause: {simulator.get_last_exit_event_cause()}")
print("=" * 70)
print("")
