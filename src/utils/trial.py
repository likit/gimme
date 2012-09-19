import networkx as nx

def create_bipartite_graph(G):
    '''Return a bipartite graph with top nodes = G.nodes and
    bottom nodes = {1,2,3...} which 1=G.nodes[1], etc.
    Each edge has capacity=1.0.

    G is a directed graph. Start and End nodes are not included
    in a bipartite graph.

    '''
    B = nx.DiGraph()
    g = G.copy()
    g.remove_nodes_from(['Start', 'End'])

    bottom_nodes = {}
    node_index = [] # note, nodes are not id'd in topological order

    node_id = 0
    for node in g.nodes():
        B.add_edge('S', node, capacity=1.0) # add edge to [S]ource node

        bottom_nodes[node] = node_id
        node_index.append(node)

        B.add_edge(node_id, 'T', capacity=1.0) # add edge to [T]arget node

        node_id += 1
    # print bottom_nodes;

    for node in g.nodes():
        for e in g[node].keys():
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
        if node == 'S':
            continue # ignore S node
        for e in edges[node]:
            if e == 'T':
                continue # ignore T node

            # print node, e

            if edges[node][e] > 0.0:
                K.add_edge(node, node_index[e])

    return K.edges()

def remove_maxflow_edges(mf_edges, B):
    '''Remove edges from max flow path of bipartite graph B.'''

    for node in mf_edges:
        if node == 'S': continue # ignore S node
        for e in mf_edges[node]:
            if e == 'T':
                continue # ignore T node
            if mf_edges[node][e] == 1.0:
                B[node][e]['capacity'] = 0.0 # disable the edge

def add_path(edges, paths, G):
    '''Build paths from a given edges.

    G is an original directed graph.

    '''
    K = nx.DiGraph()
    K.add_edges_from(edges)
    print '-' * 40
    for SG in nx.algorithms.weakly_connected_component_subgraphs(K):
        print '*' * 40
        dfs_path = [dp for dp in nx.algorithms.dfs_preorder_nodes(SG)]
        SG.add_path(nx.algorithms.shortest_path(G, 'Start', dfs_path[0])) 
        SG.add_path(nx.algorithms.shortest_path(G, dfs_path[-1], 'End'))
        print [p for p in nx.algorithms.dfs_preorder_nodes(SG, 'Start')]
        # for e in SG.edges():
        #     print e

def get_min_paths(G):
    '''Returns minimal paths including all edges.
    G is a directed graph.

    '''
    paths = []
    all_edges = set()
    B, node_index = create_bipartite_graph(G)

    mf_edges = run_max_flow(B)
    edges = set(rebuild_edges(mf_edges, node_index))

    while edges: # a max flow path is found
        all_edges.update(edges)
        add_path(edges, paths, G)
        remove_maxflow_edges(mf_edges, B)
        mf_edges = run_max_flow(B)
        edges = set(rebuild_edges(mf_edges, node_index))

    if len(all_edges) != len(G.edges()) - 2:
        print "Error: Some edges are missing."
        print set(G.edges()).difference(all_edges)

if __name__=='__main__':
    G = nx.DiGraph()
    G.add_path(['Start', 'A', 'B', 'C', 'D', 'E', 'End'])
    G.add_edge('A', 'C')
    G.add_edge('B', 'D')
    G.add_edge('C', 'E')
    paths = get_min_paths(G)
