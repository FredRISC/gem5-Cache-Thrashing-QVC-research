# src/mem/VictimL1/VictimL1.py
from m5.objects.ClockedObject import ClockedObject
from m5.params import *
from m5.proxy import *

class VictimL1(ClockedObject):
    type = 'VictimL1'
    cxx_header = "mem/VictimL1/VictimL1.hh"
    cxx_class = 'gem5::memory::VictimL1'

    # Ports
    cpu_side = ResponsePort("Port towards the CPU")
    mem_side = RequestPort("Port towards the memory system (Ruby)")

    # Configuration
    size = Param.MemorySize("128kB", "Size of the victim cache")
    assoc = Param.Int(4, "Associativity") # 4-way set associative
    block_size = Param.Int(64, "Cache line size")
    
    system = Param.System(Parent.any, "System this cache is part of")

