import igraph as ig

def get_meta_graph(g, NUM_PROCS):
    g.vs['orig_id'] = list(range(g.vcount()))

    # Extracting subgraphs for each process
    proc_nodes = [[] for i in range(NUM_PROCS+2)]
    proc_graphs = []
    for v in g.vs:
        proc_nodes[int(v['process_id'])].append(v)
    for node_list in proc_nodes:
        proc_graphs.append(g.subgraph(node_list))


    # Making meta nodes from connected components

    graph_to_meta_node_map = [-1 for i in range(g.vcount())]
    meta_node_to_graph_map = []
    meta_node_id = 0
    meta_node_process_id_map = []

    for pid, pi in enumerate(proc_graphs):
        if pid == NUM_PROCS:    # For clique nodes, want to separate the components for diff. All-to-all functions | Also leaves alltoallv as checkpoints for miniAMR benchmarking
            mpi_func_list = set()
            for v in pi.vs:
                mpi_func_list.add(v['mpi_function'])
            mpi_func_based_nodes = {f: [] for f in mpi_func_list}
            for v in pi.vs:
                mpi_func_based_nodes[v['mpi_function']].append(v.index)
            # Make diff subgraphs based on MPI funcs and make meta nodes
            for func in mpi_func_based_nodes:
                sg_nodes = mpi_func_based_nodes[func]
                sg = pi.subgraph(sg_nodes)
                comps = sg.components(mode='weak')
                for c in comps:
                    curr_meta_node_mems = []
                    for v in c:
                        graph_to_meta_node_map[sg.vs[v]['orig_id']] = meta_node_id
                        curr_meta_node_mems.append(sg.vs[v]['orig_id'])
                    meta_node_to_graph_map.append(curr_meta_node_mems)
                    meta_node_process_id_map.append(pid)
                    meta_node_id+=1
            continue

        comps = pi.components(mode = 'weak')
        for c in comps:
            curr_meta_node_mems = []
            for v in c:
                graph_to_meta_node_map[pi.vs[v]['orig_id']] = meta_node_id
                curr_meta_node_mems.append(pi.vs[v]['orig_id'])
            meta_node_to_graph_map.append(curr_meta_node_mems)
            meta_node_process_id_map.append(pid)
            meta_node_id+=1

    mpi_enc_map = {'MPI_Allreduce': 'a',
    'MPI_Barrier': 'b',
    'MPI_Bcast': 'c',
    'MPI_Finalize': 'd',
    'MPI_Init': 'e',
    'MPI_Isend': 'f',
    'MPI_Isend_Waitsome': 'g',
    'MPI_Recv': 'h',
    'MPI_Reduce': 'i',
    'MPI_Waitsome': 'j',
    'MPI_Waitany': 'k',
    'MPI_Waitall': 'l',
    'MPI_Send': 'm',
    'MPI_Wait': 'n',
    'MPI_Irecv': 'o',
    'MPI_Alltoall': 'p',
    'MPI_Alltoallv': 'q'}

    # Making mpi feature string

    meta_node_mpi_string_map = []
    meta_node_isend_seq_map = []
    meta_node_lvl_map = []
    meta_node_evg_node_map = []
    for m in meta_node_to_graph_map:
        sg = g.subgraph(m)
        order = sg.topological_sorting()    # To get the right order of MPI func as IDs are not topological order preserved
        curr_str = ""
        isend_str = ""
        evg_node_str = ""
        for n in order:
            vid = sg.vs[n]['orig_id']
            curr_str += mpi_enc_map[g.vs[vid]['mpi_function']]
            curr_str += str(int(g.vs[vid]['num_main_func']))
            evg_node_str += sg.vs[n]['ev_nodes'] if evg_node_str == "" else "_" + sg.vs[n]['ev_nodes']
            if len(sg.vs[n]['isend_seq'])>0:
                if len(isend_str) > 0:
                    isend_str += "_" + sg.vs[n]['isend_seq']
                else:
                    isend_str = sg.vs[n]['isend_seq']
        meta_node_mpi_string_map.append(curr_str)
        meta_node_isend_seq_map.append(isend_str)
        meta_node_lvl_map.append(sg.vs[order[0]]['lvl'])
        meta_node_evg_node_map.append(evg_node_str)

    # Drawing Edges

    edges = []

    for u in range(g.vcount()):
        for v in g.neighbors(u, mode='out'):
            edges.append((graph_to_meta_node_map[u], graph_to_meta_node_map[v]))

    meta_g = ig.Graph(edges = edges, directed = True)
    meta_g.simplify()

    meta_g.vs['mpi_str'] = meta_node_mpi_string_map
    meta_g.vs['process_id'] = meta_node_process_id_map
    meta_g.vs['isend_seq'] = meta_node_isend_seq_map
    meta_g.vs['lvl'] = meta_node_lvl_map
    meta_g.vs['ev_nodes'] = meta_node_evg_node_map
    return meta_g

def meta_graph_add_lvl(g, NUM_PROC):


    proc_nodes = [[] for i in range(NUM_PROC+2)]
    for v in g.vs:
        proc_nodes[int(v['process_id'])].append(v)

    for pn in proc_nodes[:NUM_PROC]:
        lvl_ord = sorted(pn, key=lambda x: x['lvl'])
        for lvl, v in enumerate(lvl_ord):
            g.vs[v.index]['lvl'] = lvl
        
    return g

def meta_graph_inject_proc_edges(g, NUM_PROC):

    proc_nodes = [[] for i in range(NUM_PROC+2)]

    for v in g.vs:
        proc_nodes[int(v['process_id'])].append(v)

    for i in range(len(proc_nodes)):
        proc_nodes[i] = sorted(proc_nodes[i], key = lambda x: int(x['lvl']))


    new_edges = []
    for pid in range(NUM_PROC):
        for i in range(len(proc_nodes[pid])-1):
            new_edges.append((proc_nodes[pid][i], proc_nodes[pid][i+1]))

    g.add_edges(new_edges)

    return g
