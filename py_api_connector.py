import sys; sys.path.append("build")
import nd_compress
import igraph as ig
from op_add_lvl_lib import op_add_lvl

def compress(graph_path, c_lvl, num_proc):
    assert(c_lvl>0 and c_lvl<=3, "c_lvl needs to be between 1 and 3")
    assert(num_proc>0, "num_proc needs to be more than 0")

    



