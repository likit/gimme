import networkx as nx

def split(graph):
    strands = set([graph.get_edge_data(*e)['strand'] for e in graph.edges()])
    g1 = nx.DiGraph()
    g2 = nx.DiGraph()
    for e in graph.edges():
        st = graph.get_edge_data(*e)['strand']
        if st == '+': g1.add_edge(*e)
        elif st == '-': g2.add_edge(*e)
        else : pass

    return strands, g1, g2
