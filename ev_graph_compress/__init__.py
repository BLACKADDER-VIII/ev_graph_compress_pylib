from . import nd_compress
import igraph as ig

def compress(graph_path, c_lvl, num_proc):
    assert 0 < c_lvl <= 3, "c_lvl must be between 1 and 3"
    assert num_proc > 0, "num_proc must be greater than 0"
    edges, (n_attr_dict, s_attr_dict) = nd_compress.compress(graph_path, c_lvl, num_proc)
    g = ig.Graph(edges, directed=True)
    for attr in n_attr_dict:
        g.vs[attr] = n_attr_dict[attr]
    for attr in s_attr_dict:
        g.vs[attr] = s_attr_dict[attr]
    return g
