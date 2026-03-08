"""
Saccharomyces phylogenetic network — graph construction.

Builds a NetworkX DiGraph covering:
- Level 1: species backbone + interspecific hybrids
- Level 2: S. cerevisiae population nodes
- Level 3: representative commercial strains mapped to populations
"""

import networkx as nx


def build_network() -> nx.DiGraph:
    G = nx.DiGraph()

    # ── Level 1: species backbone ────────────────────────────────────
    species = [
        ("outgroup", {"display_name": "N. castellii", "type": "species",
                       "is_hybrid": False, "fermentation": [], "wild": True}),
        ("S_arboricola", {"display_name": "S. arboricola", "type": "species",
                          "is_hybrid": False, "fermentation": [], "wild": True}),
        ("S_kudriavzevii", {"display_name": "S. kudriavzevii", "type": "species",
                            "is_hybrid": False, "fermentation": [], "wild": True}),
        ("S_mikatae", {"display_name": "S. mikatae", "type": "species",
                       "is_hybrid": False, "fermentation": [], "wild": True}),
        ("S_jurei", {"display_name": "S. jurei", "type": "species",
                     "is_hybrid": False, "fermentation": [], "wild": True}),
        ("S_paradoxus", {"display_name": "S. paradoxus", "type": "species",
                         "is_hybrid": False, "fermentation": [], "wild": True}),
        ("S_uvarum", {"display_name": "S. uvarum", "type": "species",
                      "is_hybrid": False, "fermentation": ["wine", "cider"], "wild": False}),
        ("S_eubayanus", {"display_name": "S. eubayanus", "type": "species",
                         "is_hybrid": False, "fermentation": [], "wild": True}),
        ("S_cerevisiae", {"display_name": "S. cerevisiae", "type": "species",
                          "is_hybrid": False, "fermentation": ["ale", "wine", "sake", "bread"],
                          "wild": False}),
    ]
    G.add_nodes_from(species)

    # Backbone divergence edges (parent → child in phylogenetic sense)
    backbone = [
        ("outgroup", "S_arboricola", {"type": "divergence", "divergence_mya": 20}),
        ("outgroup", "S_kudriavzevii", {"type": "divergence", "divergence_mya": 20}),
        ("S_kudriavzevii", "S_mikatae", {"type": "divergence"}),
        ("S_mikatae", "S_jurei", {"type": "divergence"}),
        ("S_mikatae", "S_paradoxus", {"type": "divergence"}),
        ("S_paradoxus", "S_uvarum", {"type": "divergence"}),
        ("S_uvarum", "S_eubayanus", {"type": "divergence"}),
        ("S_uvarum", "S_cerevisiae", {"type": "divergence"}),
    ]
    G.add_edges_from(backbone)

    # ── Interspecific hybrids ────────────────────────────────────────
    hybrids = [
        ("S_pastorianus", {"display_name": "S. pastorianus", "type": "hybrid",
                           "is_hybrid": True, "fermentation": ["lager"]}),
        ("Sc_Sk_hybrid", {"display_name": "S. cer. × S. kud.", "type": "hybrid",
                          "is_hybrid": True, "fermentation": ["wine"]}),
        ("Sc_Su_hybrid", {"display_name": "S. cer. × S. uva.", "type": "hybrid",
                          "is_hybrid": True, "fermentation": ["wine", "cider"]}),
        ("triple_hybrid", {"display_name": "S. cer. × kud. × uva.", "type": "hybrid",
                           "is_hybrid": True, "fermentation": ["wine"]}),
        ("S_bayanus", {"display_name": "S. bayanus CBS 380", "type": "hybrid",
                       "is_hybrid": True, "fermentation": ["wine"]}),
    ]
    G.add_nodes_from(hybrids)

    hybrid_edges = [
        ("S_cerevisiae", "S_pastorianus", {"type": "hybridisation", "admixture_fraction": 0.5}),
        ("S_eubayanus", "S_pastorianus", {"type": "hybridisation", "admixture_fraction": 0.5}),
        ("S_cerevisiae", "Sc_Sk_hybrid", {"type": "hybridisation", "admixture_fraction": 0.5}),
        ("S_kudriavzevii", "Sc_Sk_hybrid", {"type": "hybridisation", "admixture_fraction": 0.5}),
        ("S_cerevisiae", "Sc_Su_hybrid", {"type": "hybridisation", "admixture_fraction": 0.5}),
        ("S_uvarum", "Sc_Su_hybrid", {"type": "hybridisation", "admixture_fraction": 0.5}),
        ("S_cerevisiae", "triple_hybrid", {"type": "hybridisation", "admixture_fraction": 0.34}),
        ("S_kudriavzevii", "triple_hybrid", {"type": "hybridisation", "admixture_fraction": 0.33}),
        ("S_uvarum", "triple_hybrid", {"type": "hybridisation", "admixture_fraction": 0.33}),
        ("S_uvarum", "S_bayanus", {"type": "hybridisation", "admixture_fraction": 0.67}),
        ("S_eubayanus", "S_bayanus", {"type": "hybridisation", "admixture_fraction": 0.33}),
    ]
    G.add_edges_from(hybrid_edges)

    # ── Level 2: S. cerevisiae populations ───────────────────────────
    populations = [
        ("beer1", {"display_name": "Beer 1", "type": "population",
                   "is_hybrid": False, "fermentation": ["ale"]}),
        ("beer2", {"display_name": "Beer 2", "type": "population",
                   "is_hybrid": False, "fermentation": ["wheat beer", "saison"]}),
        ("wine_european", {"display_name": "Wine / European", "type": "population",
                           "is_hybrid": False, "fermentation": ["wine", "bread"]}),
        ("sake", {"display_name": "Sake", "type": "population",
                  "is_hybrid": False, "fermentation": ["sake"]}),
        ("west_african", {"display_name": "West African", "type": "population",
                          "is_hybrid": False, "fermentation": ["palm wine"]}),
        ("malaysian", {"display_name": "Malaysian", "type": "population",
                       "is_hybrid": False, "fermentation": [], "wild": True}),
        ("north_american", {"display_name": "North American", "type": "population",
                            "is_hybrid": False, "fermentation": [], "wild": True}),
        ("farmhouse", {"display_name": "Farmhouse", "type": "population",
                       "is_hybrid": False, "fermentation": ["farmhouse ale"]}),
        ("wild_root", {"display_name": "Wild (root)", "type": "population",
                       "is_hybrid": False, "fermentation": [], "wild": True}),
        ("andean_chicha", {"display_name": "Andean chicha", "type": "population",
                           "is_hybrid": False, "fermentation": ["chicha"]}),
    ]
    G.add_nodes_from(populations)

    for pop_id, _ in populations:
        G.add_edge("S_cerevisiae", pop_id, type="divergence")

    # Farmhouse has mixed Beer 1 + Asian ancestry
    G.add_edge("beer1", "farmhouse", type="introgression", admixture_fraction=0.3,
               confidence="medium")

    # S. paradoxus introgression into Neotropical populations
    G.add_edge("S_paradoxus", "andean_chicha", type="introgression",
               admixture_fraction=0.1, confidence="medium")

    # ── Level 3: representative commercial strains ───────────────────
    strains = {
        "beer1": [
            ("WLP001", "WLP001 Cal. Ale"),
            ("WY1056", "WY1056 American Ale"),
            ("US05", "Fermentis US-05"),
            ("WLP002", "WLP002 English Ale"),
            ("WY1968", "WY1968 London ESB"),
            ("OYL004", "Omega OYL-004 West Coast Ale I"),
            ("S04", "Fermentis S-04"),
        ],
        "beer2": [
            ("WLP400", "WLP400 Belgian Wit"),
            ("WY3944", "WY3944 Belgian Wit"),
            ("WLP565", "WLP565 Belgian Saison I"),
            ("WY3724", "WY3724 Belgian Saison"),
            ("OYL026", "Omega OYL-026 French Saison"),
        ],
        "wine_european": [
            ("EC1118", "Lalvin EC-1118"),
            ("D47", "Lalvin ICV-D47"),
            ("WLP715", "WLP715 Champagne"),
            ("QA23", "Lalvin QA23"),
            ("BM45", "Lalvin BM45"),
        ],
        "sake": [
            ("WLP705", "WLP705 Sake"),
            ("WY4134", "WY4134 Sake #9"),
        ],
        "farmhouse": [
            ("OYL033", "Omega OYL-033 Jovaru"),
            ("OYL036", "Omega OYL-036 Simonaitis"),
            ("OYL057", "Omega OYL-057 HotHead"),
            ("OYL061", "Omega OYL-061 Voss Kveik"),
            ("OYL071", "Omega OYL-071 Lutra"),
            ("OYL090", "Omega OYL-090 Espe Kveik"),
            ("ESCLIT", "Escarpment Ecobrau Lithuanian"),
        ],
        "west_african": [
            ("WLP300", "WLP300 Hefeweizen"),
            ("WY3068", "WY3068 Weihenstephan Weizen"),
        ],
    }
    for pop_id, strain_list in strains.items():
        for strain_id, display_name in strain_list:
            G.add_node(strain_id, display_name=display_name, type="strain",
                       is_hybrid=False)
            G.add_edge(pop_id, strain_id, type="divergence")

    return G


if __name__ == "__main__":
    G = build_network()
    print(f"Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")
    for node, data in G.nodes(data=True):
        print(f"  {node}: {data.get('display_name', node)} ({data.get('type', '?')})")
