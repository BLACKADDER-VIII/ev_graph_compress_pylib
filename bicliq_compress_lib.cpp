#include <omp.h>
#include <vector>
#include <unordered_map>
#include <iostream>
extern "C"{
#include <igraph.h>
#include <igraph_attributes.h>
}
#include <set>
#include <stdio.h>

class bicliq_compressor {
public:
    enum MPI_FUNCTION {
        MPI_Isend,
        MPI_Send,
        MPI_Recv,
        MPI_Irecv,
        MPI_Allreduce,
        MPI_Bcast,
        MPI_Reduce,
        MPI_Alltoall,
        MPI_Alltoallv,
        MPI_Wait,
        MPI_Waitall,
        MPI_Waitany,
        MPI_Waitsome,
        MPI_Barrier,
        MPI_Test,
        MPI_Init,
        MPI_Finalize
    };

private:
    int num_procs;
    igraph_integer_t graph_num_nodes;
    igraph_integer_t num_comps;

    igraph_t graph, compressed_graph, induced_sg;
    std::vector<igraph_t> comps;

    igraph_vector_t attr_proc_id;
    igraph_strvector_t attr_mpi_func;
    igraph_strvector_t attr_isend_seq;
    igraph_vector_t attr_num_main_func;
    igraph_vector_t attr_lvl;
    igraph_strvector_t attr_ev_nodes;

    igraph_vector_t cmpr_attr_proc_id;
    igraph_strvector_t cmpr_attr_mpi_func;
    igraph_vector_t cmpr_attr_num_main_funcs;
    igraph_strvector_t cmpr_attr_isend_seq;
    igraph_vector_t cmpr_attr_lvl;
    igraph_strvector_t cmpr_attr_ev_nodes;

    igraph_vector_int_t graph_to_sg_map, sg_to_graph_map;
    igraph_vector_int_t member_to_comp_map, comp_size_map;
    std::vector<std::vector<igraph_integer_t>> comp_to_member_map;
    std::vector<std::pair<igraph_vector_int_t, igraph_vector_int_t>> comp_inner_maps;
    std::vector<int> comp_is_bicliq_map;
    std::vector<igraph_integer_t> graph_to_compressed_graph_map;
    std::vector<std::vector<igraph_integer_t>> compressed_graph_to_graph_map;

    igraph_vector_int_t compressed_edges;

    std::vector<MPI_FUNCTION> mpi_func_map;

    const static std::unordered_map<std::string, MPI_FUNCTION> mpi_str_to_enum;
    const static std::unordered_map<MPI_FUNCTION, char*> mpi_enum_to_str;

    void load_attributes() {
        VASV(&graph, "mpi_function", &attr_mpi_func);
        VASV(&graph, "isend_seq", &attr_isend_seq);
        VANV(&graph, "process_id", &attr_proc_id);
        VANV(&graph, "num_main_func", &attr_num_main_func);
        VANV(&graph, "lvl", &attr_lvl);
        VASV(&graph, "ev_nodes", &attr_ev_nodes);
        mpi_func_map.resize(graph_num_nodes);
        for (igraph_integer_t i = 0; i < graph_num_nodes; i++) {
            mpi_func_map[i] = mpi_str_to_enum.at(igraph_strvector_get(&attr_mpi_func, i));
        }
    }

    void extract_induced_subgraph() {
        igraph_vs_t sg_vs;
        igraph_vector_int_t isend_wait_nodes;
        igraph_vector_int_init(&isend_wait_nodes, 0);
        for (igraph_integer_t i = 0; i < graph_num_nodes; i++) {
            if (mpi_func_map[i] == MPI_Isend   || mpi_func_map[i] == MPI_Waitsome ||
                mpi_func_map[i] == MPI_Wait     || mpi_func_map[i] == MPI_Waitall  ||
                mpi_func_map[i] == MPI_Waitany) {
                igraph_vector_int_push_back(&isend_wait_nodes, i);
            }
        }
        igraph_vs_vector(&sg_vs, &isend_wait_nodes);
        igraph_vector_int_init(&graph_to_sg_map, 0);
        igraph_vector_int_init(&sg_to_graph_map, 0);
        igraph_induced_subgraph_map(
            &graph,
            &induced_sg,
            sg_vs,
            IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH,
            &graph_to_sg_map,
            &sg_to_graph_map
        );
        for (igraph_integer_t i = 0; i < graph_num_nodes; i++) {
            VECTOR(graph_to_sg_map)[i] -= 1;
        }
    }

    void populate_components() {
        igraph_vector_int_init(&member_to_comp_map, 0);
        igraph_vector_int_init(&comp_size_map, 0);
        igraph_connected_components(
            &induced_sg,
            &member_to_comp_map,
            &comp_size_map,
            &num_comps,
            IGRAPH_WEAK
        );
        std::cout << "Number of components in induced subgraph: " << num_comps << std::endl;
        comp_to_member_map.resize(num_comps);
        for (igraph_integer_t i = 0; i < igraph_vcount(&induced_sg); i++) {
            comp_to_member_map[VECTOR(member_to_comp_map)[i]].push_back(i);
        }
        comps.resize(num_comps);
        comp_inner_maps.resize(num_comps);
        std::cout << "Maps resized..." << std::endl;
        for (igraph_integer_t c = 0; c < num_comps; c++) {
            igraph_vector_int_init(&(comp_inner_maps[c].first), 0);
            igraph_vector_int_init(&(comp_inner_maps[c].second), 0);
            if (VECTOR(comp_size_map)[c] == 1) {
                igraph_empty(&comps[c], 0, IGRAPH_DIRECTED);
                VECTOR(graph_to_sg_map)[VECTOR(sg_to_graph_map)[comp_to_member_map[c][0]]] = -1;
                continue;
            }
            igraph_vs_t vs;
            igraph_vector_int_t mems;
            igraph_vector_int_init(&mems, 0);
            for (igraph_integer_t i : comp_to_member_map[c]) {
                igraph_vector_int_push_back(&mems, i);
            }
            igraph_vs_vector(&vs, &mems);
            igraph_induced_subgraph_map(
                &induced_sg,
                &comps[c],
                vs,
                IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH,
                &comp_inner_maps[c].first,
                &comp_inner_maps[c].second
            );
        }
    }

    bool is_biclique(igraph_t& g, igraph_vector_bool_t& bi_sets) {
        igraph_bool_t is_b;
        igraph_is_bipartite(&g, &is_b, &bi_sets);
        if (is_b) {
            igraph_integer_t set_1 = 0, set_2 = 0;
            for (igraph_integer_t j = 0; j < igraph_vector_bool_size(&bi_sets); j++) {
                (VECTOR(bi_sets)[j]) ? set_1++ : set_2++;
            }
            if (igraph_ecount(&g) == set_1 * set_2 && set_1 * set_2 > 1)
                return true;
        }
        return false;
    }

    void compress_bicliq(igraph_integer_t comp_id) {
        igraph_vector_bool_t bi_sets;
        igraph_vector_bool_init(&bi_sets, 0);
        bool is_bc = is_biclique(comps[comp_id], bi_sets);
        if (is_bc)
            comp_is_bicliq_map[comp_id] = true;
    }

    void make_compressed_nodes() {
        graph_to_compressed_graph_map.resize(graph_num_nodes, -1);
        igraph_integer_t curr_node = 0;
        std::vector<igraph_integer_t> curr_vec;
        igraph_integer_t non_sg_nodes = 0, non_bq = 0, bq = 0;
        for (igraph_integer_t i = 0; i < graph_num_nodes; i++) {
            if (VECTOR(graph_to_sg_map)[i] == -1) {
                non_sg_nodes++;
                curr_vec.push_back(i);
                compressed_graph_to_graph_map.push_back(curr_vec);
                graph_to_compressed_graph_map[i] = curr_node;
                curr_node++;
                curr_vec.clear();
            }
        }
        for (igraph_integer_t i = 0; i < num_comps; i++) {
            if (igraph_vcount(&comps[i]) == 0)
                continue;
            if (!comp_is_bicliq_map[i]) {
                for (igraph_integer_t j = 0; j < (igraph_integer_t)comp_to_member_map[i].size(); j++) {
                    curr_vec.push_back(VECTOR(sg_to_graph_map)[comp_to_member_map[i][j]]);
                    compressed_graph_to_graph_map.push_back(curr_vec);
                    graph_to_compressed_graph_map[VECTOR(sg_to_graph_map)[comp_to_member_map[i][j]]] = curr_node;
                    curr_node++;
                    curr_vec.clear();
                    non_bq++;
                }
                continue;
            }
            for (igraph_integer_t j = 0; j < (igraph_integer_t)comp_to_member_map[i].size(); j++) {
                curr_vec.push_back(VECTOR(sg_to_graph_map)[comp_to_member_map[i][j]]);
                graph_to_compressed_graph_map[VECTOR(sg_to_graph_map)[comp_to_member_map[i][j]]] = curr_node;
                bq++;
            }
            compressed_graph_to_graph_map.push_back(curr_vec);
            curr_node++;
            curr_vec.clear();
        }
        std::cout << "Nodes not in induced sg: " << non_sg_nodes << " Non bicliq: " << non_bq << " Bicliq: " << bq << std::endl;
        int inv_nodes = 0;
        for (igraph_integer_t j = 0; j < graph_num_nodes; j++) {
            if (graph_to_compressed_graph_map[j] == -1)
                inv_nodes++;
        }
        std::cout << "#Invalid nodes: " << inv_nodes << std::endl;
    }

    void draw_compressed_edges() {
        igraph_vector_int_init(&compressed_edges, 0);
        igraph_integer_t u, v;
        for (igraph_integer_t e = 0; e < igraph_ecount(&graph); e++) {
            igraph_edge(&graph, e, &u, &v);
            igraph_vector_int_push_back(&compressed_edges, graph_to_compressed_graph_map[u]);
            igraph_vector_int_push_back(&compressed_edges, graph_to_compressed_graph_map[v]);
        }
    }

    void populate_compressed_attr(igraph_integer_t num_cmp_nodes) {
        igraph_vector_init(&cmpr_attr_proc_id, num_cmp_nodes);
        igraph_strvector_init(&cmpr_attr_mpi_func, num_cmp_nodes);
        igraph_vector_init(&cmpr_attr_num_main_funcs, num_cmp_nodes);
        igraph_strvector_init(&cmpr_attr_isend_seq, num_cmp_nodes);
        igraph_vector_init(&cmpr_attr_lvl, num_cmp_nodes);
        igraph_strvector_init(&cmpr_attr_ev_nodes, num_cmp_nodes);

        #pragma omp parallel for
        for (igraph_integer_t i = 0; i < num_cmp_nodes; i++) {
            if (compressed_graph_to_graph_map[i].size() == 1) {
                VECTOR(cmpr_attr_proc_id)[i] = VECTOR(attr_proc_id)[compressed_graph_to_graph_map[i][0]];
                igraph_strvector_set(&cmpr_attr_mpi_func, i, igraph_strvector_get(&attr_mpi_func, compressed_graph_to_graph_map[i][0]));
                VECTOR(cmpr_attr_num_main_funcs)[i] = VECTOR(attr_num_main_func)[compressed_graph_to_graph_map[i][0]];
                igraph_strvector_set(&cmpr_attr_isend_seq, i, igraph_strvector_get(&attr_isend_seq, compressed_graph_to_graph_map[i][0]));
                VECTOR(cmpr_attr_lvl)[i] = VECTOR(attr_lvl)[compressed_graph_to_graph_map[i][0]];
                igraph_strvector_set(&cmpr_attr_ev_nodes, i, igraph_strvector_get(&attr_ev_nodes, compressed_graph_to_graph_map[i][0]));
            } else {
                VECTOR(cmpr_attr_proc_id)[i] = num_procs + 1;
                igraph_strvector_set(&cmpr_attr_mpi_func, i, "MPI_Isend_Waitsome");
                int num_main = 0;
                std::string isend_str = "";
                std::string ev_node_str = "";
                for (int k = 0; k < (int)compressed_graph_to_graph_map[i].size(); k++) {
                    num_main += VECTOR(attr_num_main_func)[compressed_graph_to_graph_map[i][k]];
                    isend_str += igraph_strvector_get(&attr_isend_seq, compressed_graph_to_graph_map[i][k]);
                    ev_node_str += igraph_strvector_get(&attr_ev_nodes, compressed_graph_to_graph_map[i][k]);
                    if (k != (int)compressed_graph_to_graph_map[i].size() - 1) {
                        isend_str += "_";
                        ev_node_str += "_";
                    }
                }
                igraph_strvector_set(&cmpr_attr_isend_seq, i, isend_str.c_str());
                VECTOR(cmpr_attr_num_main_funcs)[i] = num_main;
                VECTOR(cmpr_attr_lvl)[i] = -1;
                igraph_strvector_set(&cmpr_attr_ev_nodes, i, ev_node_str.c_str());
            }
        }
    }

public:
    bicliq_compressor(igraph_t graph, int num_procs)
        : graph(graph), num_procs(num_procs) {}

    ~bicliq_compressor() {
        igraph_destroy(&graph);
        // igraph_destroy(&compressed_graph);
        igraph_destroy(&induced_sg);
        igraph_vector_destroy(&attr_proc_id);
        igraph_strvector_destroy(&attr_mpi_func);
        igraph_strvector_destroy(&attr_isend_seq);
        igraph_vector_destroy(&attr_num_main_func);
        igraph_vector_destroy(&attr_lvl);
        igraph_strvector_destroy(&attr_ev_nodes);
        igraph_vector_destroy(&cmpr_attr_proc_id);
        igraph_strvector_destroy(&cmpr_attr_mpi_func);
        igraph_vector_destroy(&cmpr_attr_num_main_funcs);
        igraph_strvector_destroy(&cmpr_attr_isend_seq);
        igraph_vector_destroy(&cmpr_attr_lvl);
        igraph_strvector_destroy(&cmpr_attr_ev_nodes);
        igraph_vector_int_destroy(&graph_to_sg_map);
        igraph_vector_int_destroy(&sg_to_graph_map);
        igraph_vector_int_destroy(&member_to_comp_map);
        igraph_vector_int_destroy(&comp_size_map);
        igraph_vector_int_destroy(&compressed_edges);
    }

    igraph_t compress() {
        igraph_set_attribute_table(&igraph_cattribute_table);
        igraph_strvector_init(&attr_mpi_func, 0);
        igraph_vector_init(&attr_proc_id, 0);
        igraph_vector_init(&attr_num_main_func, 0);
        igraph_strvector_init(&attr_isend_seq, 0);
        igraph_vector_init(&attr_lvl, 0);
        igraph_strvector_init(&attr_ev_nodes, 0);

        graph_num_nodes = igraph_vcount(&graph);

        load_attributes();
        std::cout << "Loaded attributes..." << std::endl;

        extract_induced_subgraph();
        std::cout << "Subgraph extracted..." << std::endl;

        populate_components();
        std::cout << "Components populated..." << std::endl;

        comp_is_bicliq_map.resize(num_comps, false);
        #pragma omp parallel for
        for (igraph_integer_t j = 0; j < num_comps; j++) {
            compress_bicliq(j);
        }
        std::cout << "Bicliques consolidated..." << std::endl;

        make_compressed_nodes();
        std::cout << "Compressed node maps made..." << std::endl;
        igraph_integer_t num_compr_nodes = compressed_graph_to_graph_map.size();
        std::cout << "#Compressed nodes: " << num_compr_nodes << std::endl;

        populate_compressed_attr(num_compr_nodes);
        std::cout << "Attributes added..." << std::endl;

        draw_compressed_edges();
        std::cout << "Compressed edges drawn..." << std::endl;

        igraph_empty(&compressed_graph, num_compr_nodes, IGRAPH_DIRECTED);
        igraph_add_edges(&compressed_graph, &compressed_edges, 0);
        igraph_simplify(&compressed_graph, 1, 1, 0);

        SETVANV(&compressed_graph, "process_id", &cmpr_attr_proc_id);
        SETVASV(&compressed_graph, "mpi_function", &cmpr_attr_mpi_func);
        SETVANV(&compressed_graph, "num_main_func", &cmpr_attr_num_main_funcs);
        SETVASV(&compressed_graph, "isend_seq", &cmpr_attr_isend_seq);
        SETVANV(&compressed_graph, "lvl", &cmpr_attr_lvl);
        SETVASV(&compressed_graph, "ev_nodes", &cmpr_attr_ev_nodes);

        return compressed_graph;
    }
};

const std::unordered_map<std::string, bicliq_compressor::MPI_FUNCTION> bicliq_compressor::mpi_str_to_enum = {
    {"MPI_Isend",    MPI_Isend},
    {"MPI_Send",     MPI_Send},
    {"MPI_Recv",     MPI_Recv},
    {"MPI_Irecv",    MPI_Irecv},
    {"MPI_Allreduce",MPI_Allreduce},
    {"MPI_Bcast",    MPI_Bcast},
    {"MPI_Reduce",   MPI_Reduce},
    {"MPI_Alltoall", MPI_Alltoall},
    {"MPI_Alltoallv",MPI_Alltoallv},
    {"MPI_Wait",     MPI_Wait},
    {"MPI_Waitall",  MPI_Waitall},
    {"MPI_Waitany",  MPI_Waitany},
    {"MPI_Waitsome", MPI_Waitsome},
    {"MPI_Barrier",  MPI_Barrier},
    {"MPI_Test",     MPI_Test},
    {"MPI_Init",     MPI_Init},
    {"MPI_Finalize", MPI_Finalize}
};

const std::unordered_map<bicliq_compressor::MPI_FUNCTION, char*> bicliq_compressor::mpi_enum_to_str = {
    {bicliq_compressor::MPI_Isend,    "MPI_Isend"},
    {bicliq_compressor::MPI_Send,     "MPI_Send"},
    {bicliq_compressor::MPI_Recv,     "MPI_Recv"},
    {bicliq_compressor::MPI_Irecv,    "MPI_Irecv"},
    {bicliq_compressor::MPI_Allreduce,"MPI_Allreduce"},
    {bicliq_compressor::MPI_Bcast,    "MPI_Bcast"},
    {bicliq_compressor::MPI_Reduce,   "MPI_Reduce"},
    {bicliq_compressor::MPI_Alltoall, "MPI_Alltoall"},
    {bicliq_compressor::MPI_Alltoallv,"MPI_Alltoallv"},
    {bicliq_compressor::MPI_Wait,     "MPI_Wait"},
    {bicliq_compressor::MPI_Waitall,  "MPI_Waitall"},
    {bicliq_compressor::MPI_Waitany,  "MPI_Waitany"},
    {bicliq_compressor::MPI_Waitsome, "MPI_Waitsome"},
    {bicliq_compressor::MPI_Barrier,  "MPI_Barrier"},
    {bicliq_compressor::MPI_Test,     "MPI_Test"},
    {bicliq_compressor::MPI_Init,     "MPI_Init"},
    {bicliq_compressor::MPI_Finalize, "MPI_Finalize"}
};
