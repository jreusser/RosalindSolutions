# Botryococcene pathway -- draw it before you model it
# pip install matplotlib networkx

from collections import defaultdict

import networkx as nx
import matplotlib.pyplot as plt


def hierarchical_layout(G: nx.DiGraph, x_scale: float = 4.5, y_scale: float = 2.5) -> dict:
    """
    Sugiyama-style layered layout for DAGs.
    Each node is placed on the deepest layer its predecessors allow,
    then nodes within each layer are spread evenly on the x-axis.
    """
    # Longest-path layering (guarantees edges only go downward)
    layers: dict[str, int] = {}
    for node in nx.topological_sort(G):
        preds = list(G.predecessors(node))
        layers[node] = (max(layers[p] for p in preds) + 1) if preds else 0

    layer_nodes: dict[int, list] = defaultdict(list)
    for node, layer in layers.items():
        layer_nodes[layer].append(node)

    pos: dict[str, tuple[float, float]] = {}
    for layer, nodes in layer_nodes.items():
        n = len(nodes)
        for i, node in enumerate(nodes):
            pos[node] = ((i - (n - 1) / 2) * x_scale, -layer * y_scale)
    return pos


def get_layout(G: nx.DiGraph) -> dict:
    """Try Graphviz dot (best for DAGs), fall back to our hierarchical layout."""
    for fn in (
        lambda: nx.nx_agraph.graphviz_layout(G, prog="dot"),
        lambda: nx.nx_pydot.graphviz_layout(G, prog="dot"),
    ):
        try:
            return fn()
        except Exception:
            pass
    return hierarchical_layout(G)

# Each tuple is: (input molecule, output molecule, enzyme that does it)
PATHWAY = [
    ("CO2 + H2O",       "G3P + Pyruvate",   "Photosynthesis"),
    ("G3P + Pyruvate",  "DXP",              "DXS"),        # rate-limiting step 1
    ("DXP",             "MEP",              "DXR"),        # rate-limiting step 2
    ("MEP",             "IPP",              "MEP pathway"),
    ("IPP",             "DMAPP",            "IDI"),
    ("IPP + DMAPP",     "FPP",              "FPPS"),
    ("FPP + FPP",       "PSPP",             "SSL-1"),      # your gene 1
    ("PSPP",            "Squalene",         "SSL-2"),      # your gene 2 (membrane)
    ("PSPP",            "Botryococcene",    "SSL-3"),      # your gene 3 (PRODUCT)
]

def draw_pathway(pathway: list[tuple[str, str, str]]) -> None:
    G = nx.DiGraph()

    for source, target, enzyme in pathway:
        G.add_edge(source, target, label=enzyme)

    pos = get_layout(G)

    # Color the product node differently
    node_colors = [
        "gold" if node == "Botryococcene"
        else "lightblue" if node in ("SSL-1", "SSL-2", "SSL-3")
        else "lightgray"
        for node in G.nodes()
    ]

    plt.figure(figsize=(16, 10))
    nx.draw(
        G, pos,
        with_labels=True,
        node_color=node_colors,
        node_size=2000,
        font_size=8,
        arrows=True,
        arrowsize=20,
    )

    edge_labels = nx.get_edge_attributes(G, "label")
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, font_size=7)

    plt.title("Botryococcene Biosynthesis Pathway in PCC 11901")
    plt.tight_layout()
    plt.savefig("pcc_11901/botryococcene_pathway.png", dpi=150)
    plt.show()
    print("Saved: pcc_11901/botryococcene_pathway.png")

if __name__ == "__main__":
    draw_pathway(PATHWAY)