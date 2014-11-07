Graphalorithm
=============

Graph Algorithm is an algorithm dedicated for interrogating and training signaling networks (represented as interaction graphs) with experimental data from stimulus - response experiments. Graph Algorithm uses Floyd - Warshall algorithm and custom Breadth - First algorithm. Given an interaction graph topology (stored in a file Network Generic) and a set of experiments in each of which some nodes were perturbed (defined in a file stimuli ) and the resulting qualitative response (activated or unchanged) of some nodes were measured (defined in a file signal ) and a matrix of experimental values between the respective stimuli and signals (defined in a file phosphos), we address the basic problem of determining a sub-graph of the given network topology that can fit the measurements for a set of scenarios the best way possible.

This problem seeks to match the network topology with measurements from a single stimulus-
response experiment. More specifically,it operates on a set of scenarios and seek to train the network
structure over all scenarios by constructing pathways between the graph compounds (nodes) that
best fit the subordinate relationships (dependencies) imposed by the experimental design. In the
end, we are sure that through our algorithm we determine, if not the optimal sub-graph that
minimizes the number of inconsistencies between measurements and predictions, the subset of the
initial network where the optimum solution lies in.
