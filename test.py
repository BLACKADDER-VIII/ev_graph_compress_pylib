import py_api_connector
l = py_api_connector.compress("/home/exouser/amr_dumpi_files/event_graph.graphml", 3, 8)
print(l.vs[0])