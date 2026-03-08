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


def add_population_nodes(G: nx.DiGraph) -> nx.DiGraph:
    """Add Level 2: S. cerevisiae population nodes and introgression edges."""

    populations = [
        {
            "id": "pop_beer1",
            "type": "population",
            "display_name": "Beer 1",
            "is_hybrid": False,
            "fermentation": ["ale"],
            "wild": False,
            "geography": ["Europe", "Global"],
            "notes": "British, Belgian abbey, many craft strains; largest domesticated ale cluster",
        },
        {
            "id": "pop_beer2",
            "type": "population",
            "display_name": "Beer 2",
            "is_hybrid": False,
            "fermentation": ["ale"],
            "wild": False,
            "geography": ["Belgium", "Europe"],
            "notes": "Belgian saison, wheat beer; distinct from Beer 1",
        },
        {
            "id": "pop_wine",
            "type": "population",
            "display_name": "Wine / European",
            "is_hybrid": False,
            "fermentation": ["wine", "bread"],
            "wild": False,
            "geography": ["Europe", "Global"],
            "notes": "Largest overall cluster; includes many commercial wine yeasts",
        },
        {
            "id": "pop_sake",
            "type": "population",
            "display_name": "Sake",
            "is_hybrid": False,
            "fermentation": ["sake"],
            "wild": False,
            "geography": ["Japan"],
            "notes": "",
        },
        {
            "id": "pop_west_african",
            "type": "population",
            "display_name": "West African",
            "is_hybrid": False,
            "fermentation": ["palm wine"],
            "wild": False,
            "geography": ["West Africa"],
            "notes": "Palm wine, local fermentation",
        },
        {
            "id": "pop_malaysian",
            "type": "population",
            "display_name": "Malaysian",
            "is_hybrid": False,
            "fermentation": [],
            "wild": True,
            "geography": ["Southeast Asia"],
            "notes": "Wild/fermentation",
        },
        {
            "id": "pop_north_american",
            "type": "population",
            "display_name": "North American",
            "is_hybrid": False,
            "fermentation": [],
            "wild": True,
            "geography": ["North America"],
            "notes": "Wild",
        },
        {
            "id": "pop_farmhouse",
            "type": "population",
            "display_name": "Farmhouse (European landrace)",
            "is_hybrid": False,
            "fermentation": ["farmhouse ale"],
            "wild": False,
            "geography": ["Norway", "Baltics"],
            "notes": "Kveik, Lithuanian, Latvian; mixed Beer 1 + Asian domesticated ancestry; Preiss et al. 2024",
        },
    ]

    # Wild basal lineages — separate node type
    wild_node = {
        "id": "wild_basal",
        "type": "wild",
        "display_name": "Wild S. cerevisiae (basal lineages)",
        "is_hybrid": False,
        "fermentation": [],
        "wild": True,
        "geography": ["East Asia", "Global"],
        "notes": "Extant wild isolates representing the deepest S. cerevisiae diversity; "
                 "mainly East Asian (China); paraphyletic — all domesticated lineages "
                 "emerge from within this wild diversity",
    }
    G.add_node(wild_node["id"], **wild_node)
    G.add_edge("S_cerevisiae", wild_node["id"],
               type="divergence",
               admixture_fraction=None,
               divergence_mya=None,
               confidence="high")

    chicha = {
        "id": "pop_chicha",
        "type": "population",
        "display_name": "Andean chicha",
        "is_hybrid": False,
        "fermentation": ["chicha"],
        "wild": False,
        "geography": ["Ecuador", "Peru"],
        "notes": "Chicha (maize beer); related to Mexican agave and French Guiana strains; carries STA1 diastatic gene",
    }
    populations.append(chicha)

    for pop in populations:
        G.add_node(pop["id"], **pop)
        G.add_edge("S_cerevisiae", pop["id"],
                    type="divergence",
                    admixture_fraction=None,
                    divergence_mya=None,
                    confidence="high")

    # Introgression from S. paradoxus into Neotropical S. cerevisiae populations
    # (notably enriched in Andean chicha strains)
    G.add_edge("S_paradoxus", "pop_chicha",
               type="introgression",
               admixture_fraction=None,
               divergence_mya=None,
               confidence="medium")

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
    populations = [n for n, d in G.nodes(data=True) if d.get("type") == "population"]
    wilds = [n for n, d in G.nodes(data=True) if d.get("type") == "wild"]

    print(f"Nodes: {G.number_of_nodes()} "
          f"({len(species)} species, {len(hybrids)} hybrids, "
          f"{len(populations)} populations, {len(wilds)} wild, "
          f"{len(ancestors)} ancestors)")
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


def plot_static(G: nx.DiGraph, outfile="network.png"):
    """Render the network as a static matplotlib image."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    fig, ax = plt.subplots(figsize=(16, 10))

    # Hierarchical layout: assign y by depth, spread x within each depth
    depths = {}
    def _assign_depth(node, d=0):
        if node not in depths or d < depths[node]:
            depths[node] = d
            for child in G.successors(node):
                edge = G.edges[node, child]
                if edge["type"] == "divergence":
                    _assign_depth(child, d + 1)

    roots = [n for n in G.nodes() if G.in_degree(n) == 0]
    for r in roots:
        _assign_depth(r)

    # Hybrids not reached via divergence — place them one level below deepest parent
    for node in G.nodes():
        if node not in depths:
            parent_depths = [depths[p] for p in G.predecessors(node) if p in depths]
            depths[node] = max(parent_depths) + 1 if parent_depths else 0

    # Group nodes by depth, assign x positions
    from collections import defaultdict
    by_depth = defaultdict(list)
    for node, d in depths.items():
        by_depth[d].append(node)

    pos = {}
    max_depth = max(depths.values()) if depths else 0
    for d, nodes in by_depth.items():
        n = len(nodes)
        for i, node in enumerate(nodes):
            x = (i - (n - 1) / 2) * 2.5
            y = -d * 2.0
            pos[node] = (x, y)

    # Colours and sizes by node type
    node_colors = []
    node_sizes = []
    labels = {}
    for node in G.nodes():
        data = G.nodes[node]
        ntype = data.get("type", "")
        if ntype == "species":
            node_colors.append("#4A90D9")
            node_sizes.append(800)
            # Short label: genus initial + species
            name = data.get("display_name", node)
            parts = name.split()
            labels[node] = f"{parts[0][0]}. {parts[1]}" if len(parts) >= 2 else name
        elif ntype == "hybrid":
            node_colors.append("#E74C3C")
            node_sizes.append(600)
            name = data.get("display_name", node)
            if "×" in name:
                labels[node] = name.split("hybrids")[0].strip() if "hybrids" in name else name
                # Shorten further
                labels[node] = labels[node].replace("S. cerevisiae", "Sc") \
                    .replace("S. kudriavzevii", "Sk") \
                    .replace("S. uvarum", "Su")
            else:
                parts = name.split()
                labels[node] = f"{parts[0][0]}. {parts[1]}" if len(parts) >= 2 else name
        elif ntype == "population":
            node_colors.append("#27AE60")
            node_sizes.append(500)
            labels[node] = data.get("display_name", node)
        elif ntype == "wild":
            node_colors.append("#8E44AD")
            node_sizes.append(600)
            labels[node] = data.get("display_name", node)
        else:  # ancestor
            node_colors.append("#95A5A6")
            node_sizes.append(200)
            labels[node] = ""

    # Draw edges by type
    div_edges = [(u, v) for u, v, d in G.edges(data=True) if d["type"] == "divergence"]
    hyb_edges = [(u, v) for u, v, d in G.edges(data=True) if d["type"] == "hybridisation"]
    int_edges = [(u, v) for u, v, d in G.edges(data=True) if d["type"] == "introgression"]

    nx.draw_networkx_edges(G, pos, edgelist=div_edges, ax=ax,
                           edge_color="#2C3E50", width=2.0, arrows=True,
                           arrowsize=15, connectionstyle="arc3,rad=0.0")
    nx.draw_networkx_edges(G, pos, edgelist=hyb_edges, ax=ax,
                           edge_color="#E74C3C", width=1.5, style="dashed",
                           arrows=True, arrowsize=12, connectionstyle="arc3,rad=0.1")
    nx.draw_networkx_edges(G, pos, edgelist=int_edges, ax=ax,
                           edge_color="#F39C12", width=1.2, style="dotted",
                           arrows=True, arrowsize=10, connectionstyle="arc3,rad=0.15")

    # Draw nodes
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors,
                           node_size=node_sizes, edgecolors="white", linewidths=1.5)

    # Labels
    nx.draw_networkx_labels(G, pos, labels, ax=ax, font_size=7,
                            font_weight="bold")

    # Edge labels for admixture fractions on hybrid edges
    hyb_edge_labels = {}
    for u, v, d in G.edges(data=True):
        if d["type"] in ("hybridisation", "introgression"):
            frac = d.get("admixture_fraction")
            if frac is not None:
                hyb_edge_labels[(u, v)] = f"{frac:.0%}"
    nx.draw_networkx_edge_labels(G, pos, hyb_edge_labels, ax=ax,
                                 font_size=6, font_color="#E74C3C")

    # Legend
    legend_items = [
        mpatches.Patch(color="#4A90D9", label="Species"),
        mpatches.Patch(color="#E74C3C", label="Hybrid"),
        mpatches.Patch(color="#27AE60", label="Population"),
        mpatches.Patch(color="#8E44AD", label="Wild"),
        mpatches.Patch(color="#95A5A6", label="Ancestor"),
        plt.Line2D([0], [0], color="#2C3E50", lw=2, label="Divergence"),
        plt.Line2D([0], [0], color="#E74C3C", lw=1.5, ls="--", label="Hybridisation"),
        plt.Line2D([0], [0], color="#F39C12", lw=1.2, ls=":", label="Introgression"),
    ]
    ax.legend(handles=legend_items, loc="lower left", fontsize=8)

    ax.set_title("Saccharomyces Phylogenetic Network", fontsize=14, fontweight="bold")
    ax.axis("off")
    plt.tight_layout()
    plt.savefig(outfile, dpi=150, bbox_inches="tight")
    print(f"Static plot saved to {outfile}")
    plt.close()


def plot_interactive(G: nx.DiGraph, outfile="network.html"):
    """Render the network as an interactive HTML file using pyvis."""
    from pyvis.network import Network

    net = Network(height="800px", width="100%", directed=True, notebook=False)
    net.barnes_hut(gravity=-8000, central_gravity=0.3, spring_length=200)

    # Add nodes
    for node, data in G.nodes(data=True):
        ntype = data.get("type", "")
        name = data.get("display_name", node)
        notes = data.get("notes", "")

        if ntype == "species":
            color = "#4A90D9"
            size = 25
            title = f"<b>{name}</b><br>{notes}"
            ferm = data.get("fermentation", [])
            geo = data.get("geography", [])
            if ferm:
                title += f"<br>Fermentation: {', '.join(ferm)}"
            if geo:
                title += f"<br>Geography: {', '.join(geo)}"
            parts = name.split()
            label = f"{parts[0][0]}. {parts[1]}" if len(parts) >= 2 else name
        elif ntype == "hybrid":
            color = "#E74C3C"
            size = 20
            title = f"<b>{name}</b><br>{notes}"
            ferm = data.get("fermentation", [])
            if ferm:
                title += f"<br>Fermentation: {', '.join(ferm)}"
            if "×" in name:
                label = name.replace("S. cerevisiae", "Sc") \
                    .replace("S. kudriavzevii", "Sk") \
                    .replace("S. uvarum", "Su") \
                    .replace(" hybrids", "") \
                    .replace(" triple", "")
            else:
                parts = name.split()
                label = f"{parts[0][0]}. {parts[1]}" if len(parts) >= 2 else name
        elif ntype == "population":
            color = "#27AE60"
            size = 18
            label = name
            title = f"<b>{name}</b><br>{notes}"
            ferm = data.get("fermentation", [])
            geo = data.get("geography", [])
            if ferm:
                title += f"<br>Fermentation: {', '.join(ferm)}"
            if geo:
                title += f"<br>Geography: {', '.join(geo)}"
        elif ntype == "wild":
            color = "#8E44AD"
            size = 22
            label = name
            title = f"<b>{name}</b><br>{notes}"
            geo = data.get("geography", [])
            if geo:
                title += f"<br>Geography: {', '.join(geo)}"
        else:  # ancestor
            color = "#95A5A6"
            size = 8
            label = ""
            title = data.get("display_name", node)

        net.add_node(node, label=label, title=title, color=color,
                     size=size, font={"size": 10})

    # Add edges
    for u, v, data in G.edges(data=True):
        etype = data.get("type", "")
        frac = data.get("admixture_fraction")

        if etype == "divergence":
            color = "#2C3E50"
            width = 2.5
            dashes = False
            title = "Divergence"
        elif etype == "hybridisation":
            color = "#E74C3C"
            width = 2.0
            dashes = True
            frac_str = f" ({frac:.0%})" if frac is not None else ""
            title = f"Hybridisation{frac_str}"
        else:  # introgression
            color = "#F39C12"
            width = 1.5
            dashes = True
            title = "Introgression"

        net.add_edge(u, v, color=color, width=width, dashes=dashes,
                     title=title, arrows="to")

    net.save_graph(outfile)
    print(f"Interactive plot saved to {outfile}")


if __name__ == "__main__":
    G = build_species_network()
    add_population_nodes(G)
    print_summary(G)
    print(f"\nExtended Newick:\n{to_extended_newick(G)}")
    plot_static(G)
    plot_interactive(G)
