import sys; sys.path.append("build")
import nd_compress
l = nd_compress.compress("/home/exouser/amr_dumpi_files/event_graph.graphml", 1, 64)
print(len(l))