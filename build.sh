rm ev_graph_compress/build/* -rf
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build