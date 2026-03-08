"""Saccharomyces phylogenetic network (species + population levels)."""

import json

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
            "strain_count": None,
            "primary_substrate": [],
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
            "strain_count": 17,
            "primary_substrate": ["oak bark"],
            "notes": "East Asian oaks; wild only",
        },
        {
            "id": "S_kudriavzevii",
            "type": "species",
            "display_name": "Saccharomyces kudriavzevii",
            "is_hybrid": False,
            "fermentation": [],
            "wild": True,
            "geography": ["Europe", "East Asia"],
            "strain_count": 27,
            "primary_substrate": ["oak bark", "decayed leaf"],
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
            "strain_count": 11,
            "primary_substrate": ["decayed leaf", "soil"],
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
            "strain_count": 6,
            "primary_substrate": ["oak bark"],
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
            "strain_count": 197,
            "primary_substrate": ["oak bark", "soil", "insect"],
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
            "strain_count": 36,
            "primary_substrate": ["grape must", "cider", "oak bark"],
            "notes": "Wine, cider, cold fermentation; partly domesticated",
        },
        {
            "id": "S_eubayanus",
            "type": "species",
            "display_name": "Saccharomyces eubayanus",
            "is_hybrid": False,
            "fermentation": [],
            "wild": True,
            "geography": ["Patagonia", "Tibet", "North America"],
            "strain_count": 152,
            "primary_substrate": ["Nothofagus bark", "soil"],
            "notes": "Wild; lager parent; hotspot of diversity in Patagonia",
        },
        {
            "id": "S_cerevisiae",
            "type": "species",
            "display_name": "Saccharomyces cerevisiae",
            "is_hybrid": False,
            "fermentation": [],
            "wild": False,
            "geography": ["Global"],
            "strain_count": 1282,
            "primary_substrate": [],
            "notes": "Species node; fermentation associations are on population-level nodes below",
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
    """Add Level 2: S. cerevisiae population subgraph with hierarchical structure.

    Topology based on Gallone et al. (2016), Peter et al. (2018),
    Liti et al. (2009), and Preiss et al. (2024). Beer 1 originated from
    admixture between Wine/European and Sake/Asian ancestors; Beer 2 is
    a separate, more divergent domestication. Wine and Beer 1 share a
    common European domesticated ancestor.
    """

    def _add(node_dict):
        G.add_node(node_dict["id"], **node_dict)

    def _div(parent, child, confidence="high"):
        G.add_edge(parent, child, type="divergence",
                   admixture_fraction=None, divergence_mya=None,
                   confidence=confidence)

    def _intro(parent, child, fraction=None, confidence="medium"):
        G.add_edge(parent, child, type="introgression",
                   admixture_fraction=fraction, divergence_mya=None,
                   confidence=confidence)

    # ── Wild basal lineages ──────────────────────────────────────────
    _add({
        "id": "wild_basal",
        "type": "wild",
        "display_name": "Wild S. cerevisiae (basal lineages)",
        "is_hybrid": False, "fermentation": [], "wild": True,
        "geography": ["East Asia", "Global"],
        "notes": "Extant wild isolates representing the deepest S. cerevisiae "
                 "diversity; mainly East Asian (China); paraphyletic — all "
                 "domesticated lineages emerge from within this wild diversity",
    })
    _div("S_cerevisiae", "wild_basal")

    # ── Independent wild/early-diverging populations ─────────────────
    for pop in [
        {"id": "pop_malaysian", "type": "population",
         "display_name": "Malaysian", "is_hybrid": False,
         "fermentation": [], "wild": True,
         "geography": ["Southeast Asia"],
         "notes": "Wild/fermentation; early-diverging lineage"},
        {"id": "pop_north_american", "type": "population",
         "display_name": "North American", "is_hybrid": False,
         "fermentation": [], "wild": True,
         "geography": ["North America"],
         "notes": "Wild"},
        {"id": "pop_west_african", "type": "population",
         "display_name": "West African", "is_hybrid": False,
         "fermentation": ["palm wine"], "wild": False,
         "geography": ["West Africa"],
         "notes": "Palm wine, local fermentation; independent domestication"},
    ]:
        _add(pop)
        _div("S_cerevisiae", pop["id"])

    # ── Asian domesticated ancestor ──────────────────────────────────
    _add({"id": "anc_asian_domestic", "type": "ancestor",
          "display_name": "anc_asian_domestic", "is_hybrid": False,
          "fermentation": [], "wild": False, "geography": ["East Asia"],
          "notes": "Ancestral Asian domesticated lineage"})
    _div("S_cerevisiae", "anc_asian_domestic")

    _add({"id": "pop_sake", "type": "population",
          "display_name": "Sake", "is_hybrid": False,
          "fermentation": ["sake"], "wild": False,
          "geography": ["Japan"],
          "notes": "Independent Asian domestication"})
    _div("anc_asian_domestic", "pop_sake")

    # ── European domesticated ancestor (Wine + Beer 1 share this) ────
    _add({"id": "anc_european_domestic", "type": "ancestor",
          "display_name": "anc_european_domestic", "is_hybrid": False,
          "fermentation": [], "wild": False, "geography": ["Europe"],
          "notes": "Common ancestor of Wine/European and Beer 1 lineages"})
    _div("S_cerevisiae", "anc_european_domestic")

    # Wine / European
    _add({"id": "pop_wine", "type": "population",
          "display_name": "Wine / European", "is_hybrid": False,
          "fermentation": ["wine", "bread"], "wild": False,
          "geography": ["Europe", "Global"],
          "notes": "Largest overall cluster; includes commercial wine yeasts "
                   "and bread strains"})
    _div("anc_european_domestic", "pop_wine")

    # Beer 1 — polyploid admixed origin (European + Asian ancestry)
    _add({"id": "pop_beer1", "type": "population",
          "display_name": "Beer 1", "is_hybrid": False,
          "fermentation": ["ale"], "wild": False,
          "geography": ["Europe", "Global"],
          "notes": "Largest domesticated ale cluster; polyploid admixed origin "
                   "from European wine and Asian domesticated ancestors; "
                   "Gallone et al. 2016"})
    _div("anc_european_domestic", "pop_beer1")
    # Asian admixture into Beer 1
    _intro("anc_asian_domestic", "pop_beer1", confidence="high")

    # ── Beer 1 subclusters ───────────────────────────────────────────
    # Branching order from Gallone et al. 2016/2019 time-calibrated
    # phylogeny of Beer 1:
    #   pop_beer1
    #   ├── beer1_wheat (basal; wheat beer + lager S. cerevisiae parent)
    #   └── beer1_ale_core
    #       ├── beer1_belgian (continental European)
    #       └── beer1_british_us
    #           ├── beer1_british
    #           └── beer1_american (diverged from British during colonisation)

    # German wheat beer — basal in Beer 1 (sister to all other ale strains;
    # the S. cerevisiae parent of lager yeasts branches from here too)
    _add({"id": "beer1_wheat", "type": "population",
          "display_name": "Beer 1 — wheat beer", "is_hybrid": False,
          "fermentation": ["ale"], "wild": False,
          "geography": ["Germany"],
          "notes": "German Hefeweizen strains; basal in Beer 1; the "
                   "S. cerevisiae parent of lager hybrids branches from "
                   "near this subclade; POF+ (4-vinylguaiacol)"})
    _div("pop_beer1", "beer1_wheat")

    # Core ale clade (after wheat beer splits off)
    _add({"id": "beer1_ale_core", "type": "ancestor",
          "display_name": "beer1_ale_core", "is_hybrid": False,
          "fermentation": [], "wild": False,
          "geography": ["Europe"],
          "notes": "Ancestor of continental European and British/US ale strains"})
    _div("pop_beer1", "beer1_ale_core")

    # Belgian / continental European subclade
    _add({"id": "beer1_belgian", "type": "population",
          "display_name": "Beer 1 — Belgian / German ale", "is_hybrid": False,
          "fermentation": ["ale"], "wild": False,
          "geography": ["Belgium", "Germany"],
          "notes": "Belgian abbey, Trappist, Kölsch, Altbier strains"})
    _div("beer1_ale_core", "beer1_belgian")

    # British + US ancestor
    _add({"id": "beer1_british_us", "type": "ancestor",
          "display_name": "beer1_british_us", "is_hybrid": False,
          "fermentation": [], "wild": False,
          "geography": ["UK"],
          "notes": "Common ancestor of British and American ale strains"})
    _div("beer1_ale_core", "beer1_british_us")

    # British ale subclade — contains several distinct subgroups
    _add({"id": "beer1_british", "type": "population",
          "display_name": "Beer 1 — British ale", "is_hybrid": False,
          "fermentation": ["ale"], "wild": False,
          "geography": ["UK", "Ireland"],
          "notes": "British ale strains; multiple subgroups including "
                   "Whitbread B, Whitbread II, Brakspear/Burton"})
    _div("beer1_british_us", "beer1_british")

    # Whitbread B subgroup
    _add({"id": "brit_whitbread_b", "type": "population",
          "display_name": "Whitbread B family", "is_hybrid": False,
          "fermentation": ["ale"], "wild": False,
          "geography": ["UK"],
          "notes": "Perhaps the most important group of British industrial "
                   "yeasts; includes Fullers, Gale's/Hales lineage"})
    _div("beer1_british", "brit_whitbread_b")

    # Whitbread II subgroup (distinct from Whitbread B)
    _add({"id": "brit_whitbread_ii", "type": "population",
          "display_name": "Whitbread II family", "is_hybrid": False,
          "fermentation": ["ale"], "wild": False,
          "geography": ["UK"],
          "notes": "Separate Whitbread lineage; includes Boddington's; "
                   "WLP017 (Vault strain)"})
    _div("beer1_british", "brit_whitbread_ii")

    # Brakspear / Burton subgroup
    _add({"id": "brit_brakspear", "type": "population",
          "display_name": "Brakspear / Burton family", "is_hybrid": False,
          "fermentation": ["ale"], "wild": False,
          "geography": ["UK"],
          "notes": "Brakspear (via Marston's) lineage; derived from Mann "
                   "yeast widely used across SE England / Greene King"})
    _div("beer1_british", "brit_brakspear")

    # Bedford / London subgroup (S-04, WLP006, WLP013)
    _add({"id": "brit_bedford", "type": "population",
          "display_name": "Bedford / London family", "is_hybrid": False,
          "fermentation": ["ale"], "wild": False,
          "geography": ["UK"],
          "notes": "WLP006 Bedford (Charles Wells), WLP013 London, S-04; "
                   "distinct from Whitbread B despite common assumption"})
    _div("beer1_british", "brit_bedford")

    # American ale subclade — two distinct Chico sub-families
    _add({"id": "beer1_american", "type": "population",
          "display_name": "Beer 1 — American ale", "is_hybrid": False,
          "fermentation": ["ale"], "wild": False,
          "geography": ["North America", "Global"],
          "notes": "Chico family; diverged from British subclade during "
                   "colonisation era; Ballantine → Siebel → Sierra Nevada; "
                   "two distinct sub-families (WLP001 vs WY1056)"})
    _div("beer1_british_us", "beer1_american")

    # WLP001 Chico sub-family
    _add({"id": "chico_wlp001", "type": "population",
          "display_name": "Chico — WLP001 sub-family", "is_hybrid": False,
          "fermentation": ["ale"], "wild": False,
          "geography": ["North America"],
          "notes": "WLP001-centred group; chromosome VIII recombination "
                   "restoring wild-type BAT1; includes Omega, Escarpment, "
                   "Gigayeast Chico variants"})
    _div("beer1_american", "chico_wlp001")

    # WY1056 Chico sub-family
    _add({"id": "chico_wy1056", "type": "population",
          "display_name": "Chico — WY1056 sub-family", "is_hybrid": False,
          "fermentation": ["ale"], "wild": False,
          "geography": ["North America"],
          "notes": "WY1056-centred group; retains BAT1 A234D mutation; "
                   "includes US-05 and Imperial A07 Flagship"})
    _div("beer1_american", "chico_wy1056")

    # ── Specific strains ─────────────────────────────────────────────
    # evidence: "sequenced" = whole-genome placed on tree
    #           "equivalent" = known same strain as a sequenced one
    #           "inferred" = placed by phenotype, history, or interdelta PCR

    strains = [
        # ── British — Whitbread B ────────────────────────────────────
        {"id": "WLP002", "parent": "brit_whitbread_b",
         "display_name": "WLP002 English Ale", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["UK"],
         "notes": "Fullers; Whitbread B family", "confidence": "high"},
        {"id": "WLP007", "parent": "brit_whitbread_b",
         "display_name": "WLP007 Dry English Ale", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["UK"],
         "notes": "Whitbread B family", "confidence": "high"},
        {"id": "WY1968", "parent": "brit_whitbread_b",
         "display_name": "WY1968 London ESB", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["UK"],
         "notes": "Supposedly Fullers; near WLP002 but not identical; "
                  "Fay et al. 2019", "confidence": "high"},
        {"id": "WY1332", "parent": "brit_whitbread_b",
         "display_name": "WY1332 Northwest Ale", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["North America"],
         "notes": "Hales of Seattle via Gale's (England); Whitbread B "
                  "lineage; near WLP041 Pacific (Redhook)",
         "confidence": "high"},
        {"id": "escarpment_vermont", "parent": "brit_whitbread_b",
         "display_name": "Escarpment Vermont Ale", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["North America"],
         "notes": "Classic NEIPA strain; closely related to WY1968; "
                  "Whitbread B group", "confidence": "high"},

        # ── British — Whitbread II ───────────────────────────────────
        {"id": "WY1098", "parent": "brit_whitbread_ii",
         "display_name": "WY1098 British Ale", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["UK"],
         "notes": "Whitbread origin; NOT equivalent to WLP007 despite "
                  "internet charts", "confidence": "high"},
        {"id": "WY1318", "parent": "brit_whitbread_ii",
         "display_name": "WY1318 London Ale III", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["UK"],
         "notes": "Boddington's origin; NEIPA/biotransformation strain; "
                  "whole-genome places it here despite interdelta PCR "
                  "suggesting German group", "confidence": "high"},

        # ── British — Brakspear / Burton ─────────────────────────────
        {"id": "WLP023", "parent": "brit_brakspear",
         "display_name": "WLP023 Burton Ale", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["UK"],
         "notes": "Brakspear / Marston's origin", "confidence": "high"},
        {"id": "WY1275", "parent": "brit_brakspear",
         "display_name": "WY1275 Thames Valley", "evidence": "equivalent",
         "fermentation": ["ale"], "geography": ["UK"],
         "notes": "Brakspear via Marston's; same origin as WLP023; "
                  "confirmed by Wyeast via mrmalty", "confidence": "high"},

        # ── British — Bedford / London ───────────────────────────────
        {"id": "WLP013", "parent": "brit_bedford",
         "display_name": "WLP013 London Ale", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["UK"],
         "notes": "Worthington White Shield origin", "confidence": "high"},
        {"id": "S04", "parent": "brit_bedford",
         "display_name": "Fermentis S-04", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["UK"],
         "notes": "Near WLP006 Bedford and WLP013; NOT Whitbread B "
                  "despite common assumption; coded CFG in Gallone",
         "confidence": "high"},

        # ── British — other (less certain sub-position) ──────────────
        {"id": "WLP004", "parent": "beer1_british",
         "display_name": "WLP004 Irish Stout", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["Ireland"],
         "notes": "Guinness strain", "confidence": "medium"},
        {"id": "WLP028", "parent": "beer1_british",
         "display_name": "WLP028 Edinburgh Ale", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["UK"],
         "notes": "McEwan's origin; distant from WY1728 despite shared "
                  "attribution; mosaic strain", "confidence": "medium"},
        {"id": "WY1028", "parent": "beer1_british",
         "display_name": "WY1028 London Ale", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["UK"],
         "notes": "Worthington White Shield; unexpectedly close to WY1728",
         "confidence": "medium"},
        {"id": "WY1728", "parent": "beer1_british",
         "display_name": "WY1728 Scottish Ale", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["UK"],
         "notes": "McEwan's; close to WY1028 not WLP028; near WLP011 "
                  "European and WLP072 French", "confidence": "medium"},
        {"id": "nottingham", "parent": "beer1_british",
         "display_name": "Lallemand Nottingham", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["UK"],
         "notes": "Beer 1 British outlier; near WLP039 East Midlands; "
                  "from Boots multi-strain culture",
         "confidence": "medium"},

        # ── American — WLP001 sub-family ─────────────────────────────
        {"id": "WLP001", "parent": "chico_wlp001",
         "display_name": "WLP001 California Ale", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["North America"],
         "notes": "Sierra Nevada 'Chico' strain", "confidence": "high"},
        {"id": "WLP090", "parent": "chico_wlp001",
         "display_name": "WLP090 San Diego Super", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["North America"],
         "notes": "Clean American ale; closely related to WLP001",
         "confidence": "high"},
        {"id": "WY1792", "parent": "chico_wlp001",
         "display_name": "WY1792 (Fat Tire / VSS)", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["North America"],
         "notes": "New Belgium strain; WLP001 sub-family",
         "confidence": "high"},
        {"id": "escarpment_cali", "parent": "chico_wlp001",
         "display_name": "Escarpment Cali Ale", "evidence": "equivalent",
         "fermentation": ["ale"], "geography": ["North America"],
         "notes": "WLP001 sub-family", "confidence": "high"},
        {"id": "OYL004", "parent": "chico_wlp001",
         "display_name": "Omega OYL-004", "evidence": "equivalent",
         "fermentation": ["ale"], "geography": ["North America"],
         "notes": "Omega Chico equivalent; WLP001 sub-family",
         "confidence": "high"},

        # ── American — WY1056 sub-family ─────────────────────────────
        {"id": "WY1056", "parent": "chico_wy1056",
         "display_name": "WY1056 American Ale", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["North America"],
         "notes": "Chico strain; NOT identical to WLP001 — retains "
                  "BAT1 A234D mutation", "confidence": "high"},
        {"id": "US05", "parent": "chico_wy1056",
         "display_name": "Fermentis US-05", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["North America"],
         "notes": "Dry Chico; WY1056 sub-family; coded CFD in Gallone",
         "confidence": "high"},
        {"id": "WY1764", "parent": "chico_wy1056",
         "display_name": "WY1764 Pacman", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["North America"],
         "notes": "Rogue's yeast; Chico derivative; = Imperial A18 "
                  "Joystick", "confidence": "high"},
        {"id": "imperial_A07", "parent": "chico_wy1056",
         "display_name": "Imperial A07 Flagship", "evidence": "equivalent",
         "fermentation": ["ale"], "geography": ["North America"],
         "notes": "WY1056 sub-family", "confidence": "high"},
        {"id": "WY1272", "parent": "beer1_american",
         "display_name": "WY1272 American Ale II", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["North America"],
         "notes": "Near Chico/US group; exact sub-family uncertain",
         "confidence": "medium"},

        # ── Belgian / German ale ─────────────────────────────────────
        {"id": "WLP530", "parent": "beer1_belgian",
         "display_name": "WLP530 Abbey Ale", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["Belgium"],
         "notes": "Westmalle origin", "confidence": "high"},
        {"id": "WY3787", "parent": "beer1_belgian",
         "display_name": "WY3787 Trappist HG", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["Belgium"],
         "notes": "Westmalle origin; near WLP530 but not identical; "
                  "Hittinger lab (yHAB43)", "confidence": "high"},
        {"id": "WLP400", "parent": "beer1_belgian",
         "display_name": "WLP400 Belgian Wit", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["Belgium"],
         "notes": "Hoegaarden origin", "confidence": "high"},
        {"id": "WLP550", "parent": "beer1_belgian",
         "display_name": "WLP550 Belgian Ale", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["Belgium"],
         "notes": "Achouffe origin", "confidence": "high"},
        {"id": "WLP500", "parent": "beer1_belgian",
         "display_name": "WLP500 Monastery Ale", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["Belgium"],
         "notes": "Chimay origin", "confidence": "high"},
        {"id": "WLP003", "parent": "beer1_belgian",
         "display_name": "WLP003 German Ale II", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["Germany"],
         "notes": "Near WY1007", "confidence": "high"},
        {"id": "WY1007", "parent": "beer1_belgian",
         "display_name": "WY1007 German Ale", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["Germany"],
         "notes": "Close to WLP003 German II; Fay et al. 2019",
         "confidence": "high"},
        {"id": "K97", "parent": "beer1_belgian",
         "display_name": "Fermentis K-97", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["Germany"],
         "notes": "German ale; near WY1007/WLP036 group",
         "confidence": "medium"},
        {"id": "WLP029", "parent": "beer1_belgian",
         "display_name": "WLP029 German Ale / Kölsch", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["Germany"],
         "notes": "Kölsch-style; in WLP002/007 vicinity per White Labs "
                  "catalogue", "confidence": "medium"},
        {"id": "WY2565", "parent": "beer1_belgian",
         "display_name": "WY2565 Kölsch", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["Germany"],
         "notes": "Close to WLP800 Pilsner; Fay et al. 2019",
         "confidence": "high"},

        # ── German wheat beer (basal Beer 1) ─────────────────────────
        {"id": "WLP300", "parent": "beer1_wheat",
         "display_name": "WLP300 Hefeweizen", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["Germany"],
         "notes": "German wheat beer; POF+ (4-vinylguaiacol)",
         "confidence": "high"},
        {"id": "WY3068", "parent": "beer1_wheat",
         "display_name": "WY3068 Weihenstephan Weizen", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["Germany"],
         "notes": "Weihenstephan 68; POF+; Hittinger lab",
         "confidence": "high"},
    ]

    for s in strains:
        parent = s.pop("parent")
        conf = s.pop("confidence", "high")
        s["type"] = "strain"
        s["is_hybrid"] = False
        s["wild"] = False
        _add(s)
        _div(parent, s["id"], confidence=conf)

    # ── Beer 2 — separate domestication, STA1+ diastatic ─────────────
    _add({"id": "pop_beer2", "type": "population",
          "display_name": "Beer 2", "is_hybrid": False,
          "fermentation": ["ale"], "wild": False,
          "geography": ["Belgium", "Europe"],
          "notes": "Belgian saison, wheat beer; separate domestication from "
                   "Beer 1; STA1 diastatic gene prevalent; lacks geographic "
                   "structure; Gallone et al. 2016"})
    _div("S_cerevisiae", "pop_beer2")

    beer2_strains = [
        # Saison group
        {"id": "WLP565", "parent": "pop_beer2",
         "display_name": "WLP565 Belgian Saison I", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["Belgium"],
         "notes": "Dupont origin; classic saison; WY3724 equivalent; STA1+",
         "confidence": "high"},
        {"id": "WY3724", "parent": "pop_beer2",
         "display_name": "WY3724 Belgian Saison", "evidence": "equivalent",
         "fermentation": ["ale"], "geography": ["Belgium"],
         "notes": "Dupont; = WLP565; STA1+", "confidence": "high"},
        {"id": "WLP566", "parent": "pop_beer2",
         "display_name": "WLP566 Belgian Saison II", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["Belgium"],
         "notes": "Second saison strain; POF- despite saison label "
                  "(PAD1/FDC1 nonsense mutations)", "confidence": "high"},
        {"id": "belle_saison", "parent": "pop_beer2",
         "display_name": "Lallemand Belle Saison", "evidence": "equivalent",
         "fermentation": ["ale"], "geography": ["Belgium"],
         "notes": "= Fermentis BE-134; STA1+ diastaticus",
         "confidence": "medium"},

        # Duvel group
        {"id": "WLP570", "parent": "pop_beer2",
         "display_name": "WLP570 Belgian Golden Ale", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["Belgium"],
         "notes": "Duvel origin; supposedly from McEwan's",
         "confidence": "high"},
        {"id": "WY1388", "parent": "pop_beer2",
         "display_name": "WY1388 Belgian Strong Ale", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["Belgium"],
         "notes": "Duvel; STA1+; clusters near WLP570",
         "confidence": "high"},
        {"id": "WB06", "parent": "pop_beer2",
         "display_name": "Fermentis WB-06", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["Belgium"],
         "notes": "Marketed as wheat beer but genetically Beer 2 / Duvel "
                  "family; STA1+; near WLP570/WY1388",
         "confidence": "high"},
        {"id": "WLP644", "parent": "pop_beer2",
         "display_name": "WLP644 Sacch. 'Trois'", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["Belgium"],
         "notes": "Originally mislabeled as Brettanomyces; actually "
                  "S. cerevisiae; from 3 Fonteinen; near Duvel strains; "
                  "POF-, STA1+", "confidence": "high"},

        # Surprise British strains in Beer 2
        {"id": "WLP026", "parent": "pop_beer2",
         "display_name": "WLP026 Premium Bitter", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["UK"],
         "notes": "Marston's; British origin but falls in Beer 2 among "
                  "Belgian saison strains", "confidence": "high"},
        {"id": "WLP037", "parent": "pop_beer2",
         "display_name": "WLP037 Yorkshire Square", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["UK"],
         "notes": "Sam Smith's; British origin but Beer 2",
         "confidence": "high"},
        {"id": "WLP038", "parent": "pop_beer2",
         "display_name": "WLP038 Manchester Ale", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["UK"],
         "notes": "Manchester; British origin but Beer 2",
         "confidence": "high"},
    ]

    for s in beer2_strains:
        parent = s.pop("parent")
        conf = s.pop("confidence", "high")
        s["type"] = "strain"
        s["is_hybrid"] = False
        s["wild"] = False
        _add(s)
        _div(parent, s["id"], confidence=conf)

    # ── Mixed clade (bottle refermentation, bread) ───────────────────
    _add({"id": "pop_mixed", "type": "population",
          "display_name": "Mixed", "is_hybrid": False,
          "fermentation": ["ale", "bread"], "wild": False,
          "geography": ["Belgium", "Europe"],
          "notes": "Atypical beer yeasts for bottle refermentation of strong "
                   "Belgian ales; contains all bread strains; mosaic ancestry; "
                   "Gallone et al. 2016"})
    _div("S_cerevisiae", "pop_mixed")

    mixed_strains = [
        {"id": "S33", "parent": "pop_mixed",
         "display_name": "Fermentis S-33", "evidence": "sequenced",
         "fermentation": ["ale", "bread"], "geography": ["UK"],
         "notes": "EDME origin; near Lallemand Windsor and bread yeasts",
         "confidence": "high"},
        {"id": "windsor", "parent": "pop_mixed",
         "display_name": "Lallemand Windsor", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["UK"],
         "notes": "From Boots multi-strain culture; very close to S-33; "
                  "near bread yeasts", "confidence": "high"},
        {"id": "T58", "parent": "pop_mixed",
         "display_name": "Fermentis T-58", "evidence": "sequenced",
         "fermentation": ["ale"], "geography": ["Belgium"],
         "notes": "Belgian strain; Mixed group, not Beer 1 or Beer 2",
         "confidence": "medium"},
    ]

    for s in mixed_strains:
        parent = s.pop("parent")
        conf = s.pop("confidence", "high")
        s["type"] = "strain"
        s["is_hybrid"] = False
        s["wild"] = False
        _add(s)
        _div(parent, s["id"], confidence=conf)

    # ── Farmhouse / kveik ────────────────────────────────────────────
    _add({"id": "pop_farmhouse", "type": "population",
          "display_name": "Farmhouse (European landrace)", "is_hybrid": False,
          "fermentation": ["farmhouse ale"], "wild": False,
          "geography": ["Norway", "Baltics"],
          "notes": "Kveik, Lithuanian, Latvian; mixed Beer 1 + Asian "
                   "domesticated ancestry; sister group to Beer 1; "
                   "majority tetraploid; Preiss et al. 2018, 2024"})
    _div("S_cerevisiae", "pop_farmhouse")
    _intro("pop_beer1", "pop_farmhouse", confidence="high")
    _intro("anc_asian_domestic", "pop_farmhouse", confidence="high")

    # Kveik sub-population
    _add({"id": "kveik_norwegian", "type": "population",
          "display_name": "Norwegian kveik", "is_hybrid": False,
          "fermentation": ["farmhouse ale"], "wild": False,
          "geography": ["Norway"],
          "notes": "Genetically distinct group; POF-; non-diastatic; "
                   "high thermotolerance; Preiss et al. 2018"})
    _div("pop_farmhouse", "kveik_norwegian")

    # Baltic sub-population
    _add({"id": "kveik_baltic", "type": "population",
          "display_name": "Baltic landrace", "is_hybrid": False,
          "fermentation": ["farmhouse ale"], "wild": False,
          "geography": ["Lithuania", "Latvia"],
          "notes": "Lithuanian and Latvian landrace strains; separate "
                   "sub-population from Norwegian kveik; Preiss et al. 2024"})
    _div("pop_farmhouse", "kveik_baltic")

    kveik_strains = [
        # ── Norwegian kveik ──────────────────────────────────────────
        {"id": "OYL061", "parent": "kveik_norwegian",
         "display_name": "Omega OYL-061 Voss Kveik", "evidence": "sequenced",
         "fermentation": ["farmhouse ale"], "geography": ["Norway"],
         "notes": "From Sigmund Gjernes (Voss); thermotolerant; "
                  "Preiss et al. 2018 sequenced",
         "confidence": "high"},
        {"id": "lallemand_voss", "parent": "kveik_norwegian",
         "display_name": "Lallemand Voss Kveik", "evidence": "equivalent",
         "fermentation": ["farmhouse ale"], "geography": ["Norway"],
         "notes": "Dry format of Voss kveik; = OYL-061 origin",
         "confidence": "high"},
        {"id": "OYL091", "parent": "kveik_norwegian",
         "display_name": "Omega OYL-091 Hornindal", "evidence": "sequenced",
         "fermentation": ["farmhouse ale"], "geography": ["Norway"],
         "notes": "Blend from Terje Raftevold (Hornindal); stone fruit "
                  "and pineapple at high temps; Preiss et al. 2018",
         "confidence": "high"},
        {"id": "OYL071", "parent": "kveik_norwegian",
         "display_name": "Omega OYL-071 Lutra", "evidence": "sequenced",
         "fermentation": ["farmhouse ale"], "geography": ["Norway"],
         "notes": "Single isolate from Hornindal blend; very clean; "
                  "pseudo-lager capability",
         "confidence": "high"},
        {"id": "OYL057", "parent": "kveik_norwegian",
         "display_name": "Omega OYL-057 HotHead", "evidence": "sequenced",
         "fermentation": ["farmhouse ale"], "geography": ["Norway"],
         "notes": "From Stranda kveik; collected by Lars Garshol; "
                  "Preiss et al. 2018",
         "confidence": "high"},
        {"id": "OYL090", "parent": "kveik_norwegian",
         "display_name": "Omega OYL-090 Espe Kveik", "evidence": "sequenced",
         "fermentation": ["farmhouse ale"], "geography": ["Norway"],
         "notes": "From Rivenes farm (Espe); tropical/citrus; "
                  "collected by Lars Garshol",
         "confidence": "high"},
        {"id": "imperial_B48", "parent": "kveik_norwegian",
         "display_name": "Imperial B48 Triple Double", "evidence": "equivalent",
         "fermentation": ["farmhouse ale"], "geography": ["Norway"],
         "notes": "Stranda kveik origin; = OYL-057 HotHead source culture",
         "confidence": "high"},
        {"id": "escarpment_krispy", "parent": "kveik_norwegian",
         "display_name": "Escarpment Krispy Kveik", "evidence": "sequenced",
         "fermentation": ["farmhouse ale"], "geography": ["Norway"],
         "notes": "Clean kveik isolate; Escarpment Labs; sequenced in "
                  "Preiss et al. 2018",
         "confidence": "high"},

        # ── Baltic landrace ──────────────────────────────────────────
        {"id": "OYL033", "parent": "kveik_baltic",
         "display_name": "Omega OYL-033 Jovaru", "evidence": "sequenced",
         "fermentation": ["farmhouse ale"], "geography": ["Lithuania"],
         "notes": "From Jovaru Alus brewery (Jovarai, Lithuania); "
                  "traditional Lithuanian farmhouse ale; POF+; "
                  "Preiss et al. 2018",
         "confidence": "high"},
        {"id": "OYL036", "parent": "kveik_baltic",
         "display_name": "Omega OYL-036 Simonaitis", "evidence": "sequenced",
         "fermentation": ["farmhouse ale"], "geography": ["Lithuania"],
         "notes": "Lithuanian farmhouse strain; from Ramunas Cizas "
                  "brewery; collected by Lars Garshol; "
                  "Preiss et al. 2018",
         "confidence": "high"},
        {"id": "escarpment_ecobrau", "parent": "kveik_baltic",
         "display_name": "Escarpment Ecobrau Lithuanian", "evidence": "sequenced",
         "fermentation": ["farmhouse ale"], "geography": ["Lithuania"],
         "notes": "Lithuanian farmhouse isolate; Escarpment Labs; "
                  "Preiss et al. 2018",
         "confidence": "high"},
    ]

    for s in kveik_strains:
        parent = s.pop("parent")
        conf = s.pop("confidence", "high")
        s["type"] = "strain"
        s["is_hybrid"] = False
        s["wild"] = False
        _add(s)
        _div(parent, s["id"], confidence=conf)

    # ── Andean chicha ────────────────────────────────────────────────
    _add({"id": "pop_chicha", "type": "population",
          "display_name": "Andean chicha", "is_hybrid": False,
          "fermentation": ["chicha"], "wild": False,
          "geography": ["Ecuador", "Peru"],
          "notes": "Chicha (maize beer); related to Mexican agave and French "
                   "Guiana strains; carries STA1 diastatic gene"})
    _div("S_cerevisiae", "pop_chicha")

    # Introgression from S. paradoxus into Neotropical populations
    _intro("S_paradoxus", "pop_chicha")

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


def to_graphml(G: nx.DiGraph, outfile="network.graphml"):
    """Export the network to GraphML format.

    List-valued node attributes are serialised as JSON strings
    since GraphML only supports scalar types.
    """
    G_export = G.copy()
    for node in G_export.nodes():
        for key, val in G_export.nodes[node].items():
            if isinstance(val, list):
                G_export.nodes[node][key] = json.dumps(val)
            elif val is None:
                G_export.nodes[node][key] = ""
    for u, v in G_export.edges():
        for key, val in G_export.edges[u, v].items():
            if val is None:
                G_export.edges[u, v][key] = ""
            elif isinstance(val, list):
                G_export.edges[u, v][key] = json.dumps(val)
    nx.write_graphml(G_export, outfile)
    print(f"GraphML exported to {outfile}")


def print_summary(G: nx.DiGraph):
    """Print a summary of the network."""
    by_type = {}
    for n, d in G.nodes(data=True):
        t = d.get("type", "unknown")
        by_type.setdefault(t, []).append(n)

    counts = ", ".join(f"{len(v)} {k}" for k, v in by_type.items())
    print(f"Nodes: {G.number_of_nodes()} ({counts})")
    print(f"Edges: {G.number_of_edges()}")

    div_edges = [(u, v) for u, v, d in G.edges(data=True) if d["type"] == "divergence"]
    hyb_edges = [(u, v) for u, v, d in G.edges(data=True) if d["type"] == "hybridisation"]
    int_edges = [(u, v) for u, v, d in G.edges(data=True) if d["type"] == "introgression"]
    print(f"  Divergence: {len(div_edges)}, Hybridisation: {len(hyb_edges)}, "
          f"Introgression: {len(int_edges)}")

    print("\nSpecies strain counts (Peris et al. 2023):")
    for s in by_type.get("species", []):
        data = G.nodes[s]
        count = data.get("strain_count")
        if count is not None:
            substrate = data.get("primary_substrate", [])
            sub_str = f"  substrate: {', '.join(substrate)}" if substrate else ""
            print(f"  {data['display_name']}: {count} strains{sub_str}")

    print("\nHybrid nodes:")
    for h in by_type.get("hybrid", []):
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
        elif ntype == "strain":
            node_colors.append("#F39C12")
            node_sizes.append(350)
            labels[node] = data.get("display_name", node)
        elif ntype == "wild":
            node_colors.append("#8E44AD")
            node_sizes.append(600)
            labels[node] = data.get("display_name", node)
        else:  # ancestor
            node_colors.append("#95A5A6")
            node_sizes.append(200)
            labels[node] = ""

    # Draw edges by type and confidence
    div_high = [(u, v) for u, v, d in G.edges(data=True)
                if d["type"] == "divergence" and d.get("confidence") == "high"]
    div_med = [(u, v) for u, v, d in G.edges(data=True)
               if d["type"] == "divergence" and d.get("confidence") == "medium"]
    div_low = [(u, v) for u, v, d in G.edges(data=True)
               if d["type"] == "divergence" and d.get("confidence") == "low"]
    hyb_edges = [(u, v) for u, v, d in G.edges(data=True) if d["type"] == "hybridisation"]
    int_edges = [(u, v) for u, v, d in G.edges(data=True) if d["type"] == "introgression"]

    nx.draw_networkx_edges(G, pos, edgelist=div_high, ax=ax,
                           edge_color="#2C3E50", width=2.0, arrows=True,
                           arrowsize=15, connectionstyle="arc3,rad=0.0")
    nx.draw_networkx_edges(G, pos, edgelist=div_med, ax=ax,
                           edge_color="#2C3E50", width=1.5, style="dashed",
                           arrows=True, arrowsize=12, connectionstyle="arc3,rad=0.0")
    nx.draw_networkx_edges(G, pos, edgelist=div_low, ax=ax,
                           edge_color="#2C3E50", width=1.0, style="dotted",
                           arrows=True, arrowsize=10, connectionstyle="arc3,rad=0.0")
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
        mpatches.Patch(color="#F39C12", label="Strain"),
        mpatches.Patch(color="#8E44AD", label="Wild"),
        mpatches.Patch(color="#95A5A6", label="Ancestor"),
        plt.Line2D([0], [0], color="#2C3E50", lw=2, label="Divergence (high)"),
        plt.Line2D([0], [0], color="#2C3E50", lw=1.5, ls="--", label="Divergence (medium)"),
        plt.Line2D([0], [0], color="#2C3E50", lw=1.0, ls=":", label="Divergence (low)"),
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
            count = data.get("strain_count")
            substrate = data.get("primary_substrate", [])
            if count is not None:
                title += f"<br>Strains: {count}"
            if substrate:
                title += f"<br>Substrate: {', '.join(substrate)}"
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
        elif ntype == "strain":
            color = "#F39C12"
            size = 14
            label = name
            evidence = data.get("evidence", "")
            title = f"<b>{name}</b><br>{notes}"
            if evidence:
                title += f"<br>Evidence: {evidence}"
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

        conf = data.get("confidence", "high")
        if etype == "divergence":
            color = "#2C3E50"
            if conf == "high":
                width = 2.5
                dashes = False
            elif conf == "medium":
                width = 1.8
                dashes = [10, 5]
            else:  # low
                width = 1.2
                dashes = [3, 5]
            title = f"Divergence ({conf})"
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

    newick = to_extended_newick(G)
    print(f"\nExtended Newick:\n{newick}")
    with open("network.nwk", "w") as f:
        f.write(newick + "\n")
    print("Extended Newick saved to network.nwk")

    to_graphml(G)
    plot_static(G)
    plot_interactive(G)
