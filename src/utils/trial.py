import networkx as nx

def create_bipartite_graph(G):
    '''Return a bipartite graph with top nodes = G.nodes and
    bottom nodes = {1,2,3...} which 1=G.nodes[1], etc.
    Each edge has capacity=1.0.

    G is a directed graph.

    '''
    B = nx.DiGraph()

    bottom_nodes = {}
    node_index = ['S']

    node_id = 1 # First node is 1
    for node in G.nodes():
        B.add_edge('S', node, capacity=1.0) # add edge to [S]ource node

        bottom_nodes[node] = node_id
        node_index.append(node)

        B.add_edge(node_id, 'T', capacity=1.0) # add edge to [T]arget node

        node_id += 1
    # print bottom_nodes;

    for node in G.nodes():
        for e in G[node].keys():
            node_id = bottom_nodes[e]
            B.add_edge(node, node_id, capacity=1.0)

    # for edge in B.edges():
    #     print edge
    return B, node_index

def run_max_flow(G):
    '''Return edges in max flow analysis.

    G is a bipartite graph created from create_bipartite_graph().

    '''
    edges = nx.max_flow_min_cost(G, 'S', 'T')
    return edges

def rebuild_edges(edges, node_index):
    '''Rebuild edges from bipartite graph to orginal directed graph.

    edges are output from run_max_flow().

    '''
    K = nx.DiGraph()
    for node in edges:
        if node == 'S': continue # ignore S node
        for e in edges[node]:
            if e == 'T': continue # ignore T node

            # print node, e

            if edges[node][e] > 0.0:
                K.add_edge(node, node_index[e])

    return K.edges()

def get_min_paths(G):
    '''Returns minimal paths including all edges.
    G is a directed graph.

    '''
    B, node_index = create_bipartite_graph(G)
    edges = rebuild_edges(run_max_flow(B), node_index)
    for e in edges:
        print e

if __name__=='__main__':
    G = nx.DiGraph()
    G.add_path(['A', 'B', 'C', 'D', 'E'])
    G.add_edge('A', 'C')
    G.add_edge('B', 'D')
    G.add_edge('C', 'E')
    paths = get_min_paths(G)
