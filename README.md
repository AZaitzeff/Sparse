# Sparse
This program runs our algorithm and 6 other ranking algorithms on a graph (highest in-degree, closeness centrality, Betweeness centrality, eigenvector centrality, Katz Centrality, and PageRank)

Usage:
python sparse.py <name> <numtopnodes> <graph to load in>

The script will create a folder named <name> and put these items in it
<name>.npy - The adjacency matrix as a numpy array
<name>.json -  upload file to zaitzeff.org/graph-transitions to see the rankings by different algorithms
<name>ranks.csv - table of rankings by the 8 algorithms
<name>norChange.csv - table of how the adjacency matrix rankings changed with bigger epsilons 
<name>rwChange.csv -- table of how the random walk matrix rankings changed with bigger epsilons 

The rankings list the <numtopnodes> nodes for each algorithm and will color the <numtopnodes> top nodes in the json file.

<graph to load in> can be used in 4 ways

1. Own graph 
usage python sparse.py <name> <numtopnodes> <filename> <numberofnodes>
example of usage: python sparse.py test 6  matrix.txt 50

<filename> should be a txt file as an edge list in the following format:
for all edges between <startnode> <endnode> the txt file should have

<startnode> <endnode>
Example of the cycle graph on three nodes
0 1
1 2
2 0


2. Random 0-1 graph
usage python sparse.py <name> <numtopnodes> -r <numberofnodes> <fill chance>
example of usage: python sparse.py test 6 -r 40 .5

creates a graph of that <numberofnodes> nodes big. Between any two nodes there is a <fill chance> chance of a directed edge between them.
<fill chance> should be between 0 and 1. Does not create simple cycles.

3. Random cluster 0-1 graph

usage python sparse.py <name> <numtopnodes> -c <numberofnodes> <numberofgroups> <alpha> <beta>
example of usage: python sparse.py test 6 -c 40 4 .2 .05

Creates a graph of <numberofnodes> nodes with <numberofgroups> groups. If two nodes are in the same cluster there is a <alpha> chance of an directed edge being between them. If two nodes are in different clusters then there is a <beta> change of an directed edge being between them.
<alpha> and <beta> should between 0 and 1. <numberofgroups> should divide <numberofnodes> evenly

4. Power Law graph
usage python sparse.py <name> <numtopnodes> -p <numberofnodes> <exponent>
example of usage: python sparse.py test 6 -r 40 2.1
Creates a power law graph of <numberofnodes> nodes using <exp> as beta. 

Note for 2, 3 and 4:
Will run until creates a strongly connected graph. Not symmetric in general.
