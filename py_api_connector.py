import numpy as np
import igraph as ig
import sys; sys.path.append("build")
import nd_compress
import meta_graph_compress_lib

def compress(graph_path, c_lvl, num_proc, meta=True):
    assert c_lvl>0 and c_lvl<=3, "c_lvl needs to be between 1 and 3"
    assert num_proc>0, "num_proc needs to be more than 0"
    edges_l, (n_attr_dict, s_attr_dict) = nd_compress.compress(graph_path, c_lvl, num_proc)
    edges = np.array(edges_l, dtype=int)
    edges = edges.reshape((-1, 2))

    g = ig.Graph(edges=edges, directed=True)
    for attr in n_attr_dict:
        g.vs[attr] = n_attr_dict[attr]
    for attr in s_attr_dict:
        g.vs[attr] = s_attr_dict[attr]
    if meta:
        g = meta_graph_compress_lib.get_meta_graph(g, num_proc)
        g = meta_graph_compress_lib.meta_graph_add_lvl(g, num_proc)
        g = meta_graph_compress_lib.meta_graph_inject_proc_edges(g, num_proc)
    return g