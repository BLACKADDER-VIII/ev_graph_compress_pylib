#include <vector>
extern "C" {
#include <igraph.h>
#include <igraph_attributes.h>
}

// Assigns a topological-sort level ('lvl') to every vertex within its process subgraph.
// Mirrors op_add_lvl() from op_add_lvl_lib.py.
static void op_add_lvl(igraph_t *g, int num_proc) {
    igraph_integer_t n = igraph_vcount(g);

    igraph_vector_t process_id_map;
    igraph_vector_init(&process_id_map, 0);
    VANV(g, "process_id", &process_id_map);

    std::vector<igraph_vector_int_t> proc_nodes(num_proc);
    for (int p = 0; p < num_proc; p++)
        igraph_vector_int_init(&proc_nodes[p], 0);

    for (igraph_integer_t i = 0; i < n; i++) {
        int pid = (int)VECTOR(process_id_map)[i];
        igraph_vector_int_push_back(&proc_nodes[pid], i);
    }

    igraph_vector_t lvls;
    igraph_vector_init(&lvls, n);
    for (igraph_integer_t i = 0; i < n; i++)
        VECTOR(lvls)[i] = -1;

    for (int p = 0; p < num_proc; p++) {
        igraph_t sg;
        igraph_vs_t vs;
        igraph_vs_vector(&vs, &proc_nodes[p]);
        igraph_induced_subgraph(g, &sg, vs, IGRAPH_SUBGRAPH_AUTO);
        igraph_vs_destroy(&vs);

        igraph_vector_int_t topol;
        igraph_vector_int_init(&topol, 0);
        igraph_topological_sorting(&sg, &topol, IGRAPH_OUT);

        igraph_integer_t sz = igraph_vector_int_size(&topol);
        for (igraph_integer_t lvl = 0; lvl < sz; lvl++) {
            igraph_integer_t sg_v   = VECTOR(topol)[lvl];
            igraph_integer_t orig_v = VECTOR(proc_nodes[p])[sg_v];
            VECTOR(lvls)[orig_v] = (double)lvl;
        }

        igraph_vector_int_destroy(&topol);
        igraph_destroy(&sg);
    }

    SETVANV(g, "lvl", &lvls);

    igraph_vector_destroy(&lvls);
    igraph_vector_destroy(&process_id_map);
    for (int p = 0; p < num_proc; p++)
        igraph_vector_int_destroy(&proc_nodes[p]);
}
