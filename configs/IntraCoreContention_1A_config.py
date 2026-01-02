# project1A IntraCoreContention.py
#
# A gem5 configuration script for a simple two-level cache hierarchy
# targeting the RISC-V ISA with the Vector (RVV) extension.
# This script is designed for a "Project 1A" in a computer
# architecture course, following gem5 v24 coding conventions.
#
# It models a system with a single RISC-V core, separate L1 instruction and
# data caches, a unified L2 cache, and a main memory bus.
#
# Usage (with a custom binary):
# gem5.opt project1a_riscv_rvv.py --cmd=<path_to_riscv_binary>
#
# Usage (with a gem5-resource binary):
# gem5.opt project1a_riscv_rvv.py --resource=rvv-matmul

import m5
from m5.objects import *
import argparse

# Import the gem5 components.
from gem5.components.boards.simple_board import SimpleBoard

from gem5.components.processors.base_cpu_core import BaseCPUCore
from gem5.components.processors.base_cpu_processor import BaseCPUProcessor

from gem5.components.cachehierarchies.ruby.mesi_two_level_cache_hierarchy import (
    MESITwoLevelCacheHierarchy,
)
from gem5.components.cachehierarchies.classic.no_cache import NoCache

from gem5.components.memory.single_channel import SingleChannelDDR3_1600

#Import simulator and utilities
from gem5.simulate.simulator import Simulator
from gem5.utils.requires import requires

# This utility is used to download and cache pre-built binaries.
from gem5.resources.resource import obtain_resource, BinaryResource

# The ISA we are targeting.
from gem5.isas import ISA


# Configuration Starts
# Check that we are running on a RISC-V host.
requires(isa_required=ISA.RISCV)

# --- 1. Define Command-Line Options ---
parser = argparse.ArgumentParser()
parser.add_argument("-resource", required=False, type=str, default="rvv-matmul")
parser.add_argument("-cmd", "--cmd", required=False, type=str, default="")
parser.add_argument("-l1i-size", "--l1i_size", required=False, type=str, default="32kB")
parser.add_argument("-l1d-size", "--l1d_size", required=False, type=str, default="32kB")
parser.add_argument("-l2-size", "--l2_size", required=False, type=str, default="256kB")
parser.add_argument("-v", "--vlen", required=False, type=int, default=128)
parser.add_argument("-e", "--elen", required=False, type=int, default=64)
parser.add_argument("-c", "--cores", required=False, type=int, default=1)


# --- 2. Setup the System Components ---

# Parse the command-line arguments
args = parser.parse_args()

# Use a simple DDR3 memory controller.
memory = SingleChannelDDR3_1600(size="2GiB")

# Enable the RISC-V Vector Extension (RVV) by setting vlen and elen.
# This is done on the RISCV CPU core object, which inherits from BaseCPUCore.
class RVVCore(BaseCPUCore):
    def __init__(self, elen, vlen, cpu_id):
        super().__init__(core=RiscvO3CPU(cpu_id=cpu_id), isa=ISA.RISCV)
        self.core.isa[0].elen = elen
        self.core.isa[0].vlen = vlen
        self.core.LQEntries = 256
        self.core.SQEntries = 256
'''
        self.core.fetchWidth = width
        self.core.decodeWidth = width
        self.core.renameWidth = width
        self.core.issueWidth = width
        self.core.wbWidth = width
        self.core.commitWidth = width

        self.core.numROBEntries = rob_size

        self.core.numPhysIntRegs = num_int_regs
        self.core.numPhysFloatRegs = num_fp_regs

        self.core.branchPred = TournamentBP()
'''


# Instantiating a list of RVV-enabled cores on the processor
processor = BaseCPUProcessor(
    cores=[RVVCore(args.elen, args.vlen, i) for i in range(args.cores)]
)


# --- 3. Define and Connect the Cache Hierarchy ---
#cache_hierarchy = NoCache()
cache_hierarchy = MESITwoLevelCacheHierarchy(
    l1d_size=args.l1d_size,
    l1d_assoc=8,
    l1i_size=args.l1i_size,
    l1i_assoc=8,
    l2_size=args.l2_size,
    l2_assoc=16,
    num_l2_banks=1,
)

# --- 4. Assemble the Board and Set the Workload ---

board = SimpleBoard(
    clk_freq="3GHz",
    processor=processor,
    memory=memory,
    cache_hierarchy=cache_hierarchy,
)

# Set the workload for the simulation.
# We support both custom binaries (--cmd) and gem5-resources (--resource).
if args.cmd:
    #binary = obtain_resource("x86-pattern-print")
    binary = BinaryResource(local_path="./configs/project/" + args.cmd)
elif args.resource:
    binary = obtain_resource(args.resource)
    
board.set_se_binary_workload(binary)

# --- 5. Run the Simulation ---

print("--- Starting gem5 simulation (RISC-V with RVV) ---")
print(f"  L1I Cache Size: {args.l1i_size}")
print(f"  L1D Cache Size: {args.l1d_size}")
print(f"  L2 Cache Size: {args.l2_size}")
print(f"  RVV VLEN: {args.vlen}, ELEN: {args.elen}")

simulator = Simulator(board=board, full_system=False)
simulator.run()

print("--- Simulation Finished! ---")