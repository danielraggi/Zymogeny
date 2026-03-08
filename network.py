"""Level 1 Saccharomyces species phylogenetic network."""

import networkx as nx


def build_species_network() -> nx.DiGraph:
    """Build the Level 1 species network: 8 species + outgroup + 5 hybrid nodes."""
    G = nx.DiGraph()

    # ── Species nodes ──────────────────────────────────────────────────

    species = [
        {
            "id": "N_castellii",
            "type": "species",
            "display_name": "Naumovozyma castellii",
            "is_hybrid": False,
            "fermentation": [],
            "wild": True,
            "geography": [],
            "notes": "Outgroup",
        },
        {
            "id": "S_arboricola",
            "type": "species",
            "display_name": "Saccharomyces arboricola",
            "is_hybrid": False,
            "fermentation": [],
            "wild": True,
            "geography": ["East Asia"],
            "notes": "East Asian oaks; wild only",
        },
        {
            "id": "S_kudriavzevii",
            "type": "species",
            "display_name": "Saccharomyces kudriavzevii",
            "is_hybrid": False,
            "fermentation": [],
            "wild": True,
            "geography": ["Europe", "Asia"],
            "notes": "European/Asian oak bark; never isolated from fermentation",
        },
        {
            "id": "S_mikatae",
            "type": "species",
            "display_name": "Saccharomyces mikatae",
            "is_hybrid": False,
            "fermentation": [],
            "wild": True,
            "geography": ["Japan"],
            "notes": "Japan; wild only",
        },
        {
            "id": "S_jurei",
            "type": "species",
            "display_name": "Saccharomyces jurei",
            "is_hybrid": False,
            "fermentation": [],
            "wild": True,
            "geography": ["Europe"],
            "notes": "European pre-Alps/UK oaks; described 2017",
        },
        {
            "id": "S_paradoxus",
            "type": "species",
            "display_name": "Saccharomyces paradoxus",
            "is_hybrid": False,
            "fermentation": ["spontaneous"],
            "wild": True,
            "geography": ["Global"],
            "notes": "Global; wild, rare spontaneous fermentation",
        },
        {
            "id": "S_uvarum",
            "type": "species",
            "display_name": "Saccharomyces uvarum",
            "is_hybrid": False,
            "fermentation": ["wine", "cider"],
            "wild": False,
            "geography": ["Global"],
            "notes": "Wine, cider, cold fermentation; partly domesticated",
        },
        {
            "id": "S_eubayanus",
            "type": "species",
            "display_name": "Saccharomyces eubayanus",
            "is_hybrid": False,
            "fermentation": [],
            "wild": True,
            "geography": ["Patagonia", "Tibet"],
            "notes": "Wild; lager parent",
        },
        {
            "id": "S_cerevisiae",
            "type": "species",
            "display_name": "Saccharomyces cerevisiae",
            "is_hybrid": False,
            "fermentation": ["ale", "wine", "bread", "sake"],
            "wild": False,
            "geography": ["Global"],
            "notes": "Primary domesticated fermentation yeast",
        },
    ]

    for sp in species:
        G.add_node(sp["id"], **sp)

    # ── Backbone divergence edges (parent → child) ─────────────────────
    # Follows the topology from the brief.

    backbone = [
        # (parent, child, divergence_mya)
        ("N_castellii", "S_arboricola", 20.0),       # outgroup → earliest branch
        ("N_castellii", "S_kudriavzevii", 20.0),     # same root split
        ("S_kudriavzevii", "S_mikatae", None),
        ("S_kudriavzevii", "S_jurei", None),
        ("S_mikatae", "S_paradoxus", None),
        ("S_paradoxus", "S_uvarum", None),
        ("S_uvarum", "S_eubayanus", None),
        ("S_eubayanus", "S_cerevisiae", None),
    ]

    # Actually, the brief shows a nested clade structure. Let me re-read:
    #   Outgroup
    #   └── sensu stricto
    #       ├── S. arboricola
    #       └──┬── S. kudriavzevii
    #          └──┬── S. mikatae
    #             ├── S. jurei
    #             └──┬── S. paradoxus
    #                └──┬── S. uvarum
    #                   └──┬── S. eubayanus
    #                      └── S. cerevisiae
    #
    # We need internal ancestor nodes to represent branching points.

    G.clear()

    # Re-add species nodes
    for sp in species:
        G.add_node(sp["id"], **sp)

    # Internal ancestor nodes (unnamed ancestral populations)
    ancestors = [
        "anc_sensu_stricto",  # root of the sensu stricto clade
        "anc_1",              # parent of S_kudriavzevii + rest
        "anc_2",              # parent of S_mikatae + (S_jurei + rest)
        "anc_3",              # parent of S_jurei + (S_paradoxus + rest)
        "anc_4",              # parent of S_paradoxus + (S_uvarum + rest)
        "anc_5",              # parent of S_uvarum + (S_eubayanus + S_cerevisiae)
        "anc_6",              # parent of S_eubayanus + S_cerevisiae
    ]
    for anc in ancestors:
        G.add_node(anc, id=anc, type="ancestor", display_name=anc,
                   is_hybrid=False, fermentation=[], wild=True,
                   geography=[], notes="Internal ancestral node")

    def _div(parent, child, mya=None):
        """Add a divergence edge."""
        G.add_edge(parent, child,
                   type="divergence",
                   admixture_fraction=None,
                   divergence_mya=mya,
                   confidence="high")

    # Outgroup → sensu stricto root (~20 Mya crown age)
    _div("N_castellii", "anc_sensu_stricto", 20.0)

    # Sensu stricto root splits
    _div("anc_sensu_stricto", "S_arboricola")
    _div("anc_sensu_stricto", "anc_1")

    _div("anc_1", "S_kudriavzevii")
    _div("anc_1", "anc_2")

    _div("anc_2", "S_mikatae")
    _div("anc_2", "anc_3")

    _div("anc_3", "S_jurei")
    _div("anc_3", "anc_4")

    _div("anc_4", "S_paradoxus")
    _div("anc_4", "anc_5")

    _div("anc_5", "S_uvarum")
    _div("anc_5", "anc_6")

    _div("anc_6", "S_eubayanus")
    _div("anc_6", "S_cerevisiae")

    # ── Hybrid / reticulation nodes ────────────────────────────────────

    hybrids = [
        {
            "id": "S_pastorianus",
            "type": "hybrid",
            "display_name": "Saccharomyces pastorianus",
            "is_hybrid": True,
            "fermentation": ["lager"],
            "wild": False,
            "geography": ["Europe", "Global"],
            "notes": "Lager beer; two independent events (Saaz/Group 1, Frohberg/Group 2)",
        },
        {
            "id": "Sc_Sk_hybrids",
            "type": "hybrid",
            "display_name": "S. cerevisiae × S. kudriavzevii hybrids",
            "is_hybrid": True,
            "fermentation": ["wine"],
            "wild": False,
            "geography": ["Central Europe"],
            "notes": "Central European wine; at least 6 independent events",
        },
        {
            "id": "Sc_Su_hybrids",
            "type": "hybrid",
            "display_name": "S. cerevisiae × S. uvarum hybrids",
            "is_hybrid": True,
            "fermentation": ["wine", "cider"],
            "wild": False,
            "geography": ["Europe"],
            "notes": "Wine and cider",
        },
        {
            "id": "Sc_Sk_Su_hybrids",
            "type": "hybrid",
            "display_name": "S. cerevisiae × S. kudriavzevii × S. uvarum triple hybrids",
            "is_hybrid": True,
            "fermentation": ["wine"],
            "wild": False,
            "geography": ["Europe"],
            "notes": "Wine; triple hybrid",
        },
        {
            "id": "S_bayanus",
            "type": "hybrid",
            "display_name": "Saccharomyces bayanus (CBS 380T)",
            "is_hybrid": True,
            "fermentation": [],
            "wild": False,
            "geography": [],
            "notes": "Historical type strain; three-way mosaic (S. uvarum ~67%, "
                     "S. eubayanus ~33%, S. cerevisiae introgression)",
        },
    ]

    for h in hybrids:
        G.add_node(h["id"], **h)

    def _hyb(parent, child, fraction, confidence="high"):
        """Add a hybridisation edge."""
        G.add_edge(parent, child,
                   type="hybridisation",
                   admixture_fraction=fraction,
                   divergence_mya=None,
                   confidence=confidence)

    def _intro(parent, child, fraction, confidence="high"):
        """Add an introgression edge."""
        G.add_edge(parent, child,
                   type="introgression",
                   admixture_fraction=fraction,
                   divergence_mya=None,
                   confidence=confidence)

    # S. pastorianus = S. cerevisiae × S. eubayanus
    _hyb("S_cerevisiae", "S_pastorianus", 0.5)
    _hyb("S_eubayanus", "S_pastorianus", 0.5)

    # S. cerevisiae × S. kudriavzevii hybrids
    _hyb("S_cerevisiae", "Sc_Sk_hybrids", 0.5)
    _hyb("S_kudriavzevii", "Sc_Sk_hybrids", 0.5)

    # S. cerevisiae × S. uvarum hybrids
    _hyb("S_cerevisiae", "Sc_Su_hybrids", 0.5, confidence="medium")
    _hyb("S_uvarum", "Sc_Su_hybrids", 0.5, confidence="medium")

    # S. cerevisiae × S. kudriavzevii × S. uvarum triple hybrids
    _hyb("S_cerevisiae", "Sc_Sk_Su_hybrids", 0.34, confidence="medium")
    _hyb("S_kudriavzevii", "Sc_Sk_Su_hybrids", 0.33, confidence="medium")
    _hyb("S_uvarum", "Sc_Sk_Su_hybrids", 0.33, confidence="medium")

    # S. bayanus = S. uvarum (~67%) + S. eubayanus (~33%) + S. cerevisiae introgression
    _hyb("S_uvarum", "S_bayanus", 0.67)
    _hyb("S_eubayanus", "S_bayanus", 0.33)
    _intro("S_cerevisiae", "S_bayanus", None)

    return G


def to_extended_newick(G: nx.DiGraph) -> str:
    """Export the network to Extended Newick (Rich Newick) format.

    Hybrid nodes (multiple parents) are labelled with #H suffixes.
    Only the first occurrence expands the subtree; subsequent references
    use the #H label alone.
    """
    roots = [n for n in G.nodes() if G.in_degree(n) == 0]
    if len(roots) != 1:
        raise ValueError(f"Expected exactly one root, found {len(roots)}: {roots}")

    root = roots[0]

    # Assign #H labels to hybrid nodes (nodes with >1 parent)
    hybrid_labels = {}
    h_counter = 1
    for node in G.nodes():
        if G.in_degree(node) > 1:
            hybrid_labels[node] = f"#H{h_counter}"
            h_counter += 1

    visited_hybrids = set()

    def _newick(node):
        children = list(G.successors(node))
        data = G.nodes[node]
        label = data.get("display_name", node).replace(" ", "_")

        h_tag = hybrid_labels.get(node, "")

        # If this is a hybrid we've already expanded, just emit the label
        if node in hybrid_labels and node in visited_hybrids:
            return f"{label}{h_tag}"

        if node in hybrid_labels:
            visited_hybrids.add(node)

        # Skip ancestor nodes in the newick label — they're structural only
        if data.get("type") == "ancestor":
            label = ""

        if not children:
            return f"{label}{h_tag}"

        child_strs = []
        for child in children:
            edge = G.edges[node, child]
            subtree = _newick(child)
            branch_len = edge.get("divergence_mya")
            if branch_len is not None:
                subtree += f":{branch_len}"
            child_strs.append(subtree)

        return f"({','.join(child_strs)}){label}{h_tag}"

    return _newick(root) + ";"


def print_summary(G: nx.DiGraph):
    """Print a summary of the network."""
    species = [n for n, d in G.nodes(data=True) if d.get("type") == "species"]
    hybrids = [n for n, d in G.nodes(data=True) if d.get("type") == "hybrid"]
    ancestors = [n for n, d in G.nodes(data=True) if d.get("type") == "ancestor"]

    print(f"Nodes: {G.number_of_nodes()} "
          f"({len(species)} species, {len(hybrids)} hybrids, {len(ancestors)} ancestors)")
    print(f"Edges: {G.number_of_edges()}")

    div_edges = [(u, v) for u, v, d in G.edges(data=True) if d["type"] == "divergence"]
    hyb_edges = [(u, v) for u, v, d in G.edges(data=True) if d["type"] == "hybridisation"]
    int_edges = [(u, v) for u, v, d in G.edges(data=True) if d["type"] == "introgression"]
    print(f"  Divergence: {len(div_edges)}, Hybridisation: {len(hyb_edges)}, "
          f"Introgression: {len(int_edges)}")

    print("\nHybrid nodes:")
    for h in hybrids:
        parents = list(G.predecessors(h))
        parent_names = [G.nodes[p].get("display_name", p) for p in parents]
        print(f"  {G.nodes[h]['display_name']}")
        for p, pname in zip(parents, parent_names):
            edge = G.edges[p, h]
            frac = edge.get("admixture_fraction")
            frac_str = f" ({frac:.0%})" if frac is not None else ""
            print(f"    ← {pname}{frac_str} [{edge['type']}]")


if __name__ == "__main__":
    G = build_species_network()
    print_summary(G)
    print(f"\nExtended Newick:\n{to_extended_newick(G)}")
