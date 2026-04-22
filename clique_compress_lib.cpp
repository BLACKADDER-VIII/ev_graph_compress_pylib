#include <omp.h>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
extern "C"{
#include <igraph.h>
#include <igraph_attributes.h>
}
#include <set>
#include <stdio.h>

class clique_compressor {
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

    igraph_t graph, compressed_graph, alltoall_sg;

    igraph_strvector_t mpi_function_map;
    igraph_vector_t process_id_map;
    igraph_strvector_t isend_seq;
    igraph_vector_t num_main_func;
    igraph_vector_t lvl;
    igraph_strvector_t ev_nodes;

    igraph_integer_t num_nodes;
    igraph_integer_t num_proc;

    igraph_vector_int_t graph_to_sg_map;
    igraph_vector_int_t sg_to_graph_map;

    igraph_vector_int_t membership;
    igraph_vector_int_t csize;
    igraph_integer_t num_comps;

    igraph_vector_int_t alltoall_nodes;
    igraph_vs_t alltoall_vs;

    igraph_integer_t sg_num_nodes;
    igraph_integer_t num_compr_nodes;

    std::vector<std::vector<igraph_integer_t>> comp_members;
    std::vector<igraph_t> components;
    std::vector<std::pair<igraph_vector_int_t, igraph_vector_int_t>> comp_inner_maps;

    std::vector<std::pair<igraph_integer_t, igraph_integer_t>> sg_to_clique;
    std::vector<igraph_integer_t> clique_to_sg;
    std::vector<std::vector<std::vector<igraph_integer_t>>> clique_members;

    std::vector<std::string> compr_ev_nodes;

    igraph_vector_int_t compressed_edges;
    std::vector<igraph_integer_t> graph_to_compressed_nodes_map, compressed_nodes_to_graph_map;

    std::vector<MPI_FUNCTION> orig_graph_mpi_map;

    igraph_strvector_t attr_mpi_function;
    igraph_vector_t attr_proc_id;
    igraph_vector_t attr_num_main_func;
    igraph_strvector_t attr_isend_seq;
    igraph_vector_t attr_lvl;
    igraph_strvector_t attr_compr_ev_nodes;

    const static std::unordered_map<std::string, MPI_FUNCTION> mpi_str_to_enum;
    const static std::unordered_map<MPI_FUNCTION, char*> mpi_enum_to_str;

    void load_components(igraph_integer_t comp_id) {
        igraph_vector_int_init(&(comp_inner_maps[comp_id].first), 0);
        igraph_vector_int_init(&(comp_inner_maps[comp_id].second), 0);
        igraph_vector_int_t vertices;
        igraph_vs_t v_s;
        igraph_vector_int_init(&vertices, comp_members[comp_id].size());
        for (igraph_integer_t i = 0; i < (igraph_integer_t)comp_members[comp_id].size(); i++) {
            VECTOR(vertices)[i] = comp_members[comp_id][i];
        }
        igraph_vs_vector(&v_s, &vertices);
        igraph_induced_subgraph_map(
            &alltoall_sg,
            &components[comp_id],
            v_s,
            IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH,
            &comp_inner_maps[comp_id].first,
            &comp_inner_maps[comp_id].second
        );
    }

    void process_component(igraph_integer_t comp_id) {
        igraph_vector_int_t neighbors;
        igraph_vector_int_init(&neighbors, 0);
        igraph_t& curr_comp = components[comp_id];
        std::vector<igraph_integer_t> curr_cliq_mems;
        igraph_integer_t remain_mems = comp_members[comp_id].size();
        std::vector<bool> is_cliqued(remain_mems, false);
        igraph_integer_t curr_v = 0;
        igraph_integer_t num_neigh, curr_neigh;
        std::vector<igraph_integer_t> curr_cliq_sg_nodes;
        while (remain_mems > 0) {
            curr_cliq_mems.clear();
            curr_cliq_sg_nodes.clear();
            curr_cliq_mems.push_back(curr_v);
            igraph_neighbors(&curr_comp, &neighbors, curr_v, IGRAPH_OUT);
            num_neigh = igraph_vector_int_size(&neighbors);
            int curr_proc = VECTOR(process_id_map)[VECTOR(sg_to_graph_map)[VECTOR(comp_inner_maps[comp_id].second)[curr_v]]];
            for (igraph_integer_t i = 0; i < num_neigh; i++) {
                curr_neigh = VECTOR(neighbors)[i];
                if (VECTOR(process_id_map)[VECTOR(sg_to_graph_map)[VECTOR(comp_inner_maps[comp_id].second)[curr_neigh]]] != curr_proc)
                    curr_cliq_mems.push_back(curr_neigh);
            }
            for (auto v : curr_cliq_mems) {
                is_cliqued[v] = true;
                curr_cliq_sg_nodes.push_back(VECTOR(comp_inner_maps[comp_id].second)[v]);
            }
            clique_members[comp_id].push_back(curr_cliq_sg_nodes);
            remain_mems -= curr_cliq_mems.size();
            for (igraph_integer_t j = 0; j < (igraph_integer_t)comp_members[comp_id].size(); j++) {
                if (!is_cliqued[j]) {
                    curr_v = j;
                    break;
                }
            }
            if (remain_mems < 0)
                std::cerr << "Remain mems became negative" << std::endl;
        }
    }

    void populate_compressed_nodes() {
        igraph_integer_t num_compressed_nodes = 0;
        std::cout << "Num nodes: " << num_nodes << std::endl;
        int nodes_processed = 0;
        for (igraph_integer_t i = 0; i < num_nodes; i++) {
            igraph_integer_t n = VECTOR(graph_to_sg_map)[i];
            if (n == -1) {
                nodes_processed++;
                graph_to_compressed_nodes_map[i] = num_compressed_nodes;
                compressed_nodes_to_graph_map.push_back(i);
                igraph_strvector_push_back(&attr_mpi_function, mpi_enum_to_str.at(orig_graph_mpi_map[i]));
                igraph_vector_push_back(&attr_proc_id, VECTOR(process_id_map)[i]);
                igraph_strvector_push_back(&attr_isend_seq, igraph_strvector_get(&isend_seq, i));
                igraph_vector_push_back(&attr_num_main_func, VECTOR(num_main_func)[i]);
                igraph_vector_push_back(&attr_lvl, VECTOR(lvl)[i]);
                num_compressed_nodes++;
            }
        }
        std::cout << "# nodes outside cliques: " << nodes_processed << std::endl;
        int num_cliques_processed = 0;
        int clique_nodes_compressed = 0;
        for (auto kv : clique_members) {
            for (auto v : kv) {
                num_cliques_processed++;
                for (auto mem : v) {
                    clique_nodes_compressed++;
                    nodes_processed++;
                    graph_to_compressed_nodes_map[VECTOR(sg_to_graph_map)[mem]] = num_compressed_nodes;
                }
                compressed_nodes_to_graph_map.push_back(VECTOR(sg_to_graph_map)[v[0]]);
                igraph_strvector_push_back(&attr_mpi_function, mpi_enum_to_str.at(orig_graph_mpi_map[VECTOR(sg_to_graph_map)[v[0]]]));
                igraph_vector_push_back(&attr_proc_id, num_proc);
                igraph_strvector_push_back(&attr_isend_seq, "");
                igraph_vector_push_back(&attr_num_main_func, 1);
                igraph_vector_push_back(&attr_lvl, -1);
                num_compressed_nodes++;
            }
        }
        std::cout << "#Cliques: " << num_cliques_processed << " #nodes compressed: " << clique_nodes_compressed << std::endl;
        std::cout << "Total nodes processed: " << nodes_processed << std::endl;
    }

    void write_compressed_edges() {
        igraph_vector_int_t neighbors;
        igraph_integer_t num_edges = 0;
        igraph_vector_int_init(&neighbors, 0);
        igraph_integer_t u, v;
        igraph_integer_t num_neigh;
        for (igraph_integer_t i = 0; i < num_nodes; i++) {
            igraph_neighbors(&graph, &neighbors, i, IGRAPH_OUT);
            u = graph_to_compressed_nodes_map[i];
            num_neigh = igraph_vector_int_size(&neighbors);
            for (igraph_integer_t j = 0; j < num_neigh; j++) {
                v = graph_to_compressed_nodes_map[VECTOR(neighbors)[j]];
                igraph_vector_int_push_back(&compressed_edges, u);
                igraph_vector_int_push_back(&compressed_edges, v);
                num_edges++;
            }
        }
        std::cout << num_edges << " edges written..." << std::endl;
    }

    std::vector<std::string> get_compressed_node_to_ev_node_mapping(igraph_integer_t num_compressed_nodes) {
        std::vector<std::string> mapping(num_compressed_nodes, "");
        for (igraph_integer_t i = 0; i < num_nodes; i++) {
            igraph_integer_t compressed_node = graph_to_compressed_nodes_map[i];
            std::string ev_node = igraph_strvector_get(&ev_nodes, i);
            if (mapping[compressed_node] == "") {
                mapping[compressed_node] = ev_node;
            } else {
                mapping[compressed_node] += "_" + ev_node;
            }
        }
        return mapping;
    }

public:
    clique_compressor(igraph_t graph, int num_procs)
        : graph(graph), num_proc(num_procs) {}

    ~clique_compressor() {
        igraph_destroy(&graph);
        igraph_destroy(&compressed_graph);
        igraph_destroy(&alltoall_sg);
        igraph_strvector_destroy(&mpi_function_map);
        igraph_vector_destroy(&process_id_map);
        igraph_strvector_destroy(&isend_seq);
        igraph_vector_destroy(&num_main_func);
        igraph_vector_destroy(&lvl);
        igraph_strvector_destroy(&ev_nodes);
        igraph_vector_int_destroy(&graph_to_sg_map);
        igraph_vector_int_destroy(&sg_to_graph_map);
        igraph_vector_int_destroy(&membership);
        igraph_vector_int_destroy(&csize);
        igraph_vector_int_destroy(&alltoall_nodes);
        igraph_vector_int_destroy(&compressed_edges);
        igraph_strvector_destroy(&attr_mpi_function);
        igraph_vector_destroy(&attr_proc_id);
        igraph_vector_destroy(&attr_num_main_func);
        igraph_strvector_destroy(&attr_isend_seq);
        igraph_vector_destroy(&attr_lvl);
        igraph_strvector_destroy(&attr_compr_ev_nodes);
    }

    igraph_t compress() {
        igraph_set_attribute_table(&igraph_cattribute_table);
        igraph_strvector_init(&mpi_function_map, 0);
        igraph_vector_init(&process_id_map, 0);
        igraph_vector_int_init(&sg_to_graph_map, 0);
        igraph_strvector_init(&isend_seq, 0);
        igraph_strvector_init(&attr_isend_seq, 0);
        igraph_vector_init(&num_main_func, 0);
        igraph_vector_init(&lvl, 0);
        igraph_strvector_init(&ev_nodes, 0);
        igraph_strvector_init(&attr_mpi_function, 0);
        igraph_vector_init(&attr_proc_id, 0);
        igraph_vector_init(&attr_num_main_func, 0);
        igraph_vector_init(&attr_lvl, 0);
        igraph_strvector_init(&attr_compr_ev_nodes, 0);

        num_nodes = igraph_vcount(&graph);
        VASV(&graph, "mpi_function", &mpi_function_map);
        VANV(&graph, "process_id", &process_id_map);
        VASV(&graph, "isend_seq", &isend_seq);
        VANV(&graph, "num_main_func", &num_main_func);
        VANV(&graph, "lvl", &lvl);
        VASV(&graph, "ev_nodes", &ev_nodes);
        std::cout << "Loaded attributes..." << std::endl;

        igraph_vector_int_init(&membership, 0);
        igraph_vector_int_init(&csize, 0);
        igraph_vector_int_init(&graph_to_sg_map, 0);

        orig_graph_mpi_map.resize(num_nodes);
        #pragma omp parallel for
        for (igraph_integer_t i = 0; i < num_nodes; i++) {
            orig_graph_mpi_map[i] = mpi_str_to_enum.at(igraph_strvector_get(&mpi_function_map, i));
        }
        std::cout << "MPI functions enumerated..." << std::endl;

        igraph_vector_int_init(&alltoall_nodes, 0);
        #pragma omp parallel for
        for (igraph_integer_t i = 0; i < num_nodes; i++) {
            if (orig_graph_mpi_map[i] == MPI_Allreduce || orig_graph_mpi_map[i] == MPI_Alltoall || orig_graph_mpi_map[i] == MPI_Alltoallv) {
                #pragma omp critical
                igraph_vector_int_push_back(&alltoall_nodes, i);
            }
        }
        std::cout << "Alltoall communication nodes isolated..." << std::endl;

        igraph_vs_vector(&alltoall_vs, &alltoall_nodes);
        sg_num_nodes = igraph_vector_int_size(&alltoall_nodes);
        sg_to_clique.resize(sg_num_nodes);
        std::cout << "Subgraph identity map initialized..." << std::endl;

        igraph_induced_subgraph_map(
            &graph,
            &alltoall_sg,
            alltoall_vs,
            IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH,
            &graph_to_sg_map,
            &sg_to_graph_map
        );

        for (igraph_integer_t i = 0; i < num_nodes; i++) {
            VECTOR(graph_to_sg_map)[i] -= 1;
        }
        std::cout << "Subgraph extracted..." << std::endl;

        igraph_connected_components(
            &alltoall_sg,
            &membership,
            &csize,
            &num_comps,
            IGRAPH_STRONG
        );
        std::cout << "Connected components isolated..." << std::endl;

        comp_members.resize(num_comps);
        comp_inner_maps.resize(num_comps);
        components.resize(num_comps);
        clique_members.resize(num_comps);
        std::cout << "Component vectors initialized..." << std::endl;

        for (igraph_integer_t i = 0; i < sg_num_nodes; i++) {
            comp_members[VECTOR(membership)[i]].push_back(i);
        }
        std::cout << "Component memberships populated..." << std::endl;

        for (igraph_integer_t i = 0; i < num_comps; i++) {
            load_components(i);
        }
        std::cout << "Loaded components..." << std::endl;

        #pragma omp parallel for
        for (igraph_integer_t i = 0; i < num_comps; i++) {
            process_component(i);
        }
        std::cout << "Processed components..." << std::endl;

        igraph_vector_int_init(&compressed_edges, 0);
        graph_to_compressed_nodes_map.resize(num_nodes);
        populate_compressed_nodes();
        write_compressed_edges();

        num_compr_nodes = compressed_nodes_to_graph_map.size();

        compr_ev_nodes = get_compressed_node_to_ev_node_mapping(num_compr_nodes);
        igraph_strvector_init(&attr_compr_ev_nodes, num_compr_nodes);
        for (igraph_integer_t i = 0; i < num_compr_nodes; i++) {
            igraph_strvector_set(&attr_compr_ev_nodes, i, compr_ev_nodes[i].c_str());
        }

        igraph_empty(&compressed_graph, num_compr_nodes, IGRAPH_DIRECTED);
        igraph_add_edges(&compressed_graph, &compressed_edges, 0);

        SETVANV(&compressed_graph, "process_id", &attr_proc_id);
        SETVASV(&compressed_graph, "mpi_function", &attr_mpi_function);
        SETVASV(&compressed_graph, "isend_seq", &attr_isend_seq);
        SETVANV(&compressed_graph, "num_main_func", &attr_num_main_func);
        SETVANV(&compressed_graph, "lvl", &attr_lvl);
        SETVASV(&compressed_graph, "ev_nodes", &attr_compr_ev_nodes);

        igraph_simplify(&compressed_graph, 1, 1, 0);
        return compressed_graph;
    }
};

const std::unordered_map<std::string, clique_compressor::MPI_FUNCTION> clique_compressor::mpi_str_to_enum = {
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

const std::unordered_map<clique_compressor::MPI_FUNCTION, char*> clique_compressor::mpi_enum_to_str = {
    {clique_compressor::MPI_Isend,    "MPI_Isend"},
    {clique_compressor::MPI_Send,     "MPI_Send"},
    {clique_compressor::MPI_Recv,     "MPI_Recv"},
    {clique_compressor::MPI_Irecv,    "MPI_Irecv"},
    {clique_compressor::MPI_Allreduce,"MPI_Allreduce"},
    {clique_compressor::MPI_Bcast,    "MPI_Bcast"},
    {clique_compressor::MPI_Reduce,   "MPI_Reduce"},
    {clique_compressor::MPI_Alltoall, "MPI_Alltoall"},
    {clique_compressor::MPI_Alltoallv,"MPI_Alltoallv"},
    {clique_compressor::MPI_Wait,     "MPI_Wait"},
    {clique_compressor::MPI_Waitall,  "MPI_Waitall"},
    {clique_compressor::MPI_Waitany,  "MPI_Waitany"},
    {clique_compressor::MPI_Waitsome, "MPI_Waitsome"},
    {clique_compressor::MPI_Barrier,  "MPI_Barrier"},
    {clique_compressor::MPI_Test,     "MPI_Test"},
    {clique_compressor::MPI_Init,     "MPI_Init"},
    {clique_compressor::MPI_Finalize, "MPI_Finalize"}
};
