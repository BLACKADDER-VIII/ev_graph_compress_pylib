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

class op_compressor {
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

    std::vector<std::vector<igraph_integer_t>> process_to_node_vec;
    std::vector<std::vector<igraph_integer_t>> compressed_node_vec;
    std::vector<int> node_to_process_map;
    std::vector<igraph_integer_t> node_to_compressed_node_map;

    igraph_t graph, compressed_graph;
    igraph_strvector_t mpi_function_map;
    igraph_vector_t process_id_map;

    igraph_vector_t atr_proc_id;
    igraph_vector_t atr_num_funcs;
    igraph_strvector_t atr_mpi_func;
    igraph_strvector_t atr_isend_seq;
    igraph_strvector_t atr_compr_ev_nodes;

    std::vector<std::vector<char*>> compr_func_name;
    std::vector<std::vector<int>> compr_num_funcs;
    std::vector<std::vector<int>> compr_proc_id;
    std::vector<std::vector<std::vector<int>>> compr_isend_seq;
    std::vector<std::string> compr_ev_nodes;

    igraph_vector_int_t compressed_edges;
    igraph_integer_t graph_num_nodes;

    std::vector<MPI_FUNCTION> mpi_function_map_encoded;

    const static std::unordered_map<std::string, MPI_FUNCTION> mpi_func_encode_map;
    const static std::unordered_map<MPI_FUNCTION, char*> mpi_encode_to_string_map;

    bool check_is_collective(MPI_FUNCTION func) {
        return func == MPI_Allreduce || func == MPI_Reduce || func == MPI_Alltoall ||
               func == MPI_Alltoallv || func == MPI_Bcast  || func == MPI_Barrier  ||
               func == MPI_Init      || func == MPI_Finalize;
    }

    void populate_process_nodes() {
        process_to_node_vec.resize(num_procs);
        compressed_node_vec.resize(num_procs);
        for (igraph_integer_t i = 0; i < graph_num_nodes; i++) {
            process_to_node_vec[VECTOR(process_id_map)[i]].push_back(i);
        }
    }

    void populate_compressed_nodes(int proc) {
        std::vector<igraph_integer_t> node_chain;
        std::vector<igraph_integer_t> curr_node_grouping;
        bool is_grouped = false;
        MPI_FUNCTION curr_func;
        int vid;

        auto flush_curr_group = [&]() {
            vid = node_chain.size();
            node_chain.push_back(vid);
            compr_func_name[proc].push_back(mpi_encode_to_string_map.at(curr_func));
            compr_proc_id[proc].push_back(proc);
            compr_num_funcs[proc].push_back(curr_node_grouping.size());
            std::vector<int> isend_seq;
            is_grouped = false;
            for (igraph_integer_t i : curr_node_grouping) {
                node_to_compressed_node_map[i] = vid;
            }
            if (curr_func == MPI_Isend) {
                igraph_vector_int_t isend_neigh;
                igraph_vector_int_init(&isend_neigh, 0);
                for (igraph_integer_t i : curr_node_grouping) {
                    #pragma omp critical
                    {
                        igraph_neighbors(&graph, &isend_neigh, i, IGRAPH_OUT);
                    }
                    for (int j = 0; j < igraph_vector_int_size(&isend_neigh); j++) {
                        if (VECTOR(process_id_map)[VECTOR(isend_neigh)[j]] != proc) {
                            isend_seq.push_back(VECTOR(process_id_map)[VECTOR(isend_neigh)[j]]);
                            break;
                        }
                    }
                }
                igraph_vector_int_destroy(&isend_neigh);
            }
            compr_isend_seq[proc].push_back(isend_seq);
            curr_node_grouping.clear();
        };

        for (auto v : process_to_node_vec[proc]) {
            if (is_grouped) {
                if (mpi_function_map_encoded[v] == curr_func) {
                    curr_node_grouping.push_back(v);
                    continue;
                }
                flush_curr_group();
            }
            if (check_is_collective(mpi_function_map_encoded[v])) {
                vid = node_chain.size();
                node_chain.push_back(vid);
                compr_func_name[proc].push_back(mpi_encode_to_string_map.at(mpi_function_map_encoded[v]));
                compr_proc_id[proc].push_back(proc);
                compr_num_funcs[proc].push_back(1);
                std::vector<int> isend_seq;
                compr_isend_seq[proc].push_back(isend_seq);
                node_to_compressed_node_map[v] = vid;
                continue;
            }
            curr_func = mpi_function_map_encoded[v];
            is_grouped = true;
            curr_node_grouping.push_back(v);
        }
        compressed_node_vec[proc] = node_chain;
    }

    void adjust_compressed_node_labels(int proc) {
        igraph_integer_t shifter = 0;
        for (int i = 0; i < proc; i++) {
            shifter += compressed_node_vec[i].size();
        }
        for (auto& i : compressed_node_vec[proc]) {
            i += shifter;
        }
        for (igraph_integer_t i : process_to_node_vec[proc]) {
            node_to_compressed_node_map[i] += shifter;
        }
        for (igraph_integer_t i = 0; i < (igraph_integer_t)compr_func_name[proc].size(); i++) {
            igraph_strvector_set(&atr_mpi_func, i + shifter, compr_func_name[proc][i]);
            igraph_vector_set(&atr_num_funcs, i + shifter, compr_num_funcs[proc][i]);
            igraph_vector_set(&atr_proc_id, i + shifter, proc);

            std::vector<int> isend_seq = compr_isend_seq[proc][i];
            std::string isend_str = "";
            if (isend_seq.size() > 0) {
                isend_str += std::to_string(isend_seq[0]);
                for (int j = 1; j < (int)isend_seq.size(); j++) {
                    isend_str += "_" + std::to_string(isend_seq[j]);
                }
            }
            igraph_strvector_set(&atr_isend_seq, i + shifter, isend_str.c_str());
        }
    }

    void populate_compressed_edges(int proc) {
        igraph_vector_int_t neighbors;
        igraph_vector_int_init(&neighbors, 0);
        for (auto i : process_to_node_vec[proc]) {
            #pragma omp critical
            igraph_neighbors(&graph, &neighbors, i, IGRAPH_OUT);
            igraph_integer_t u = node_to_compressed_node_map[i];
            int num_n = igraph_vector_int_size(&neighbors);
            for (int j = 0; j < num_n; j++) {
                igraph_integer_t v = node_to_compressed_node_map[VECTOR(neighbors)[j]];
                #pragma omp critical
                {
                    igraph_vector_int_push_back(&compressed_edges, u);
                    igraph_vector_int_push_back(&compressed_edges, v);
                }
            }
        }
    }

    std::vector<std::string> get_compressed_node_to_ev_node_mapping(igraph_integer_t num_compressed_nodes) {
        std::vector<std::string> mapping(num_compressed_nodes, "");
        for (igraph_integer_t i = 0; i < graph_num_nodes; i++) {
            igraph_integer_t compressed_node = node_to_compressed_node_map[i];
            std::string ev_node = std::to_string(i);
            if (mapping[compressed_node] == "") {
                mapping[compressed_node] = ev_node;
            } else {
                mapping[compressed_node] += "_" + ev_node;
            }
        }
        return mapping;
    }

public:
    op_compressor(igraph_t graph, int num_procs)
        : graph(graph), num_procs(num_procs) {}

    ~op_compressor() {
        igraph_destroy(&graph);
        igraph_destroy(&compressed_graph);
        igraph_strvector_destroy(&mpi_function_map);
        igraph_vector_destroy(&process_id_map);
        igraph_vector_destroy(&atr_proc_id);
        igraph_vector_destroy(&atr_num_funcs);
        igraph_strvector_destroy(&atr_mpi_func);
        igraph_strvector_destroy(&atr_isend_seq);
        igraph_strvector_destroy(&atr_compr_ev_nodes);
        igraph_vector_int_destroy(&compressed_edges);
    }

    igraph_t compress() {
        igraph_set_attribute_table(&igraph_cattribute_table);
        igraph_strvector_init(&mpi_function_map, 0);
        igraph_vector_init(&process_id_map, 0);

        VASV(&graph, "mpi_function", &mpi_function_map);
        VANV(&graph, "process_id", &process_id_map);
        std::cout << "Loaded attributes..." << std::endl;

        graph_num_nodes = igraph_vcount(&graph);
        node_to_process_map.resize(graph_num_nodes);

        for (igraph_integer_t i = 0; i < graph_num_nodes; i++) {
            mpi_function_map_encoded.push_back(
                mpi_func_encode_map.at(igraph_strvector_get(&mpi_function_map, i)));
        }

        compr_func_name.resize(num_procs);
        compr_num_funcs.resize(num_procs);
        compr_proc_id.resize(num_procs);
        compr_isend_seq.resize(num_procs);
        std::cout << "Created MPI encoding..." << std::endl;

        for (igraph_integer_t i = 0; i < graph_num_nodes; i++) {
            node_to_compressed_node_map.push_back(-1);
        }
        std::cout << "Initialized compressed node map..." << std::endl;

        populate_process_nodes();
        std::cout << "Initialized process to node map..." << std::endl;

        #pragma omp parallel for
        for (int proc = 0; proc < num_procs; proc++) {
            populate_compressed_nodes(proc);
            std::cout << "Populated for process #" << proc << std::endl;
        }

        int compressed_num_nodes = 0;
        for (auto v : compressed_node_vec)
            compressed_num_nodes += v.size();

        igraph_vector_init(&atr_num_funcs, compressed_num_nodes);
        igraph_vector_init(&atr_proc_id, compressed_num_nodes);
        igraph_strvector_init(&atr_mpi_func, compressed_num_nodes);
        std::cout << "Initializing isend_seq..." << std::endl;
        igraph_strvector_init(&atr_isend_seq, compressed_num_nodes);

        #pragma omp parallel for
        for (int p = 0; p < num_procs; p++) {
            adjust_compressed_node_labels(p);
            std::cout << "Adjusted node labels..." << std::endl;
        }

        igraph_vector_int_init(&compressed_edges, 0);

        #pragma omp parallel for
        for (int p = 0; p < num_procs; p++) {
            populate_compressed_edges(p);
            std::cout << "Created edges for process #" << p << std::endl;
        }

        compr_ev_nodes = get_compressed_node_to_ev_node_mapping(compressed_num_nodes);
        igraph_strvector_init(&atr_compr_ev_nodes, compressed_num_nodes);
        for (igraph_integer_t i = 0; i < compressed_num_nodes; i++) {
            igraph_strvector_set(&atr_compr_ev_nodes, i, compr_ev_nodes[i].c_str());
        }

        igraph_empty(&compressed_graph, compressed_num_nodes, IGRAPH_DIRECTED);
        igraph_add_edges(&compressed_graph, &compressed_edges, 0);

        SETVANV(&compressed_graph, "process_id", &atr_proc_id);
        SETVASV(&compressed_graph, "mpi_function", &atr_mpi_func);
        SETVANV(&compressed_graph, "num_main_func", &atr_num_funcs);
        SETVASV(&compressed_graph, "isend_seq", &atr_isend_seq);
        SETVASV(&compressed_graph, "ev_nodes", &atr_compr_ev_nodes);

        igraph_simplify(&compressed_graph, 1, 1, 0);
        return compressed_graph;
    }
};

const std::unordered_map<std::string, op_compressor::MPI_FUNCTION> op_compressor::mpi_func_encode_map = {
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

const std::unordered_map<op_compressor::MPI_FUNCTION, char*> op_compressor::mpi_encode_to_string_map = {
    {op_compressor::MPI_Isend,    "MPI_Isend"},
    {op_compressor::MPI_Send,     "MPI_Send"},
    {op_compressor::MPI_Recv,     "MPI_Recv"},
    {op_compressor::MPI_Irecv,    "MPI_Irecv"},
    {op_compressor::MPI_Allreduce,"MPI_Allreduce"},
    {op_compressor::MPI_Bcast,    "MPI_Bcast"},
    {op_compressor::MPI_Reduce,   "MPI_Reduce"},
    {op_compressor::MPI_Alltoall, "MPI_Alltoall"},
    {op_compressor::MPI_Alltoallv,"MPI_Alltoallv"},
    {op_compressor::MPI_Wait,     "MPI_Wait"},
    {op_compressor::MPI_Waitall,  "MPI_Waitall"},
    {op_compressor::MPI_Waitany,  "MPI_Waitany"},
    {op_compressor::MPI_Waitsome, "MPI_Waitsome"},
    {op_compressor::MPI_Barrier,  "MPI_Barrier"},
    {op_compressor::MPI_Test,     "MPI_Test"},
    {op_compressor::MPI_Init,     "MPI_Init"},
    {op_compressor::MPI_Finalize, "MPI_Finalize"}
};
