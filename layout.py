"""
Sugiyama-style hierarchical DAG layout with edge-type weighting.

Key improvements over naive Sugiyama:
1. Edge weights: divergence edges (primary relationships) are weighted 3×
   more than hybridisation/introgression, so primary parent-child edges
   stay short and straight while secondary edges are allowed to be longer.
2. Sibling compaction: children of the same parent are grouped
   contiguously and centred around the parent.
3. Iterative coordinate refinement: 12 down+up passes with weighted
   median placement.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Any

import networkx as nx

# Edge type → weight for barycenter and coordinate calculations.
# Higher weight = stronger pull (edge kept shorter/straighter).
EDGE_WEIGHT = {
    "divergence": 3.0,
    "hybridisation": 1.0,
    "introgression": 0.5,
}


def _edge_weight(G: nx.DiGraph, u: str, v: str) -> float:
    """Get the layout weight for an edge."""
    data = G.edges.get((u, v), {})
    return EDGE_WEIGHT.get(data.get("type", "divergence"), 1.0)


# ─────────────────────────────────────────────────────────────────────
# Phase 1: Layer assignment
# ─────────────────────────────────────────────────────────────────────

def assign_layers(G: nx.DiGraph) -> dict[str, int]:
    """
    Assign each node to a layer using longest-path-from-source.

    Only divergence and hybridisation edges drive the layer assignment.
    Introgression edges are treated as secondary — they don't push a
    node to a later layer.  This prevents e.g. Farmhouse (which has an
    introgression from Beer 1) from being pushed a layer later than
    the other population nodes.
    """
    layer = {}
    for n in nx.topological_sort(G):
        # Only consider "structural" parents (divergence/hybridisation)
        structural_preds = [
            p for p in G.predecessors(n)
            if G.edges[p, n].get("type", "divergence") != "introgression"
        ]
        if not structural_preds:
            # Fall back to all predecessors if there are no structural ones
            # but still no predecessors at all → root node
            all_preds = [p for p in G.predecessors(n) if p in layer]
            if not all_preds:
                layer[n] = 0
            else:
                layer[n] = max(layer[p] for p in all_preds) + 1
        else:
            layer[n] = max(layer[p] for p in structural_preds) + 1
    return layer


# ─────────────────────────────────────────────────────────────────────
# Phase 2: Dummy nodes for long edges
# ─────────────────────────────────────────────────────────────────────

def insert_dummy_nodes(G: nx.DiGraph, layers: dict[str, int]):
    """Insert dummy nodes so every edge connects adjacent layers."""
    G2 = G.copy()
    layers2 = dict(layers)
    dummy_id = 0

    for u, v, data in list(G.edges(data=True)):
        span = layers[v] - layers[u]
        if span > 1:
            G2.remove_edge(u, v)
            prev = u
            for i in range(1, span):
                d = f"__dummy_{dummy_id}"
                dummy_id += 1
                G2.add_node(d, display_name="", type="dummy", is_dummy=True)
                layers2[d] = layers[u] + i
                G2.add_edge(prev, d, **data)
                prev = d
            G2.add_edge(prev, v, **data)

    return G2, layers2


# ─────────────────────────────────────────────────────────────────────
# Phase 3: Crossing minimisation
# ─────────────────────────────────────────────────────────────────────

def _weighted_barycenter(G: nx.DiGraph, node: str,
                         layer_order: dict[str, float],
                         direction: str) -> float:
    """
    Weighted barycenter: each neighbour's position is multiplied by the
    edge weight, so divergence edges pull harder than introgression.
    """
    if direction == "down":
        neighbours = [(n, _edge_weight(G, n, node))
                      for n in G.predecessors(node) if n in layer_order]
    else:
        neighbours = [(n, _edge_weight(G, node, n))
                      for n in G.successors(node) if n in layer_order]
    if not neighbours:
        return layer_order.get(node, 0.0)
    total_w = sum(w for _, w in neighbours)
    if total_w == 0:
        return layer_order.get(node, 0.0)
    return sum(layer_order[n] * w for n, w in neighbours) / total_w


def _count_crossings(G: nx.DiGraph, layer_a: list[str],
                     layer_b: list[str]) -> int:
    """Count edge crossings between two adjacent layers."""
    pos_a = {n: i for i, n in enumerate(layer_a)}
    pos_b = {n: i for i, n in enumerate(layer_b)}
    edges = []
    for u in layer_a:
        for v in G.successors(u):
            if v in pos_b:
                edges.append((pos_a[u], pos_b[v]))
    crossings = 0
    for i in range(len(edges)):
        for j in range(i + 1, len(edges)):
            if (edges[i][0] - edges[j][0]) * (edges[i][1] - edges[j][1]) < 0:
                crossings += 1
    return crossings


def _total_crossings(G, lbr):
    ranks = sorted(lbr.keys())
    return sum(_count_crossings(G, lbr[ranks[i]], lbr[ranks[i+1]])
               for i in range(len(ranks) - 1))


def _swap_reduce(G: nx.DiGraph, fixed: list[str], free: list[str],
                 fixed_is_upper: bool) -> list[str]:
    """
    Greedily swap adjacent pairs in `free` to reduce crossings with
    `fixed`.  `fixed_is_upper=True` means fixed is layer_a, free is
    layer_b; otherwise reversed.
    """
    result = list(free)
    improved = True
    while improved:
        improved = False
        for i in range(len(result) - 1):
            if fixed_is_upper:
                cur = _count_crossings(G, fixed, result)
            else:
                cur = _count_crossings(G, result, fixed)
            swapped = list(result)
            swapped[i], swapped[i + 1] = swapped[i + 1], swapped[i]
            if fixed_is_upper:
                new = _count_crossings(G, fixed, swapped)
            else:
                new = _count_crossings(G, swapped, fixed)
            if new < cur:
                result = swapped
                improved = True
    return result


def minimise_crossings(G: nx.DiGraph, layers_by_rank: dict[int, list[str]],
                       sweeps: int = 24) -> dict[int, list[str]]:
    """Weighted barycenter + adjacent-swap, alternating up/down sweeps."""
    ranks = sorted(layers_by_rank.keys())
    best = {r: list(layers_by_rank[r]) for r in ranks}
    best_xings = _total_crossings(G, best)
    current = {r: list(best[r]) for r in ranks}

    for sweep in range(sweeps):
        if sweep % 2 == 0:
            for ri in range(1, len(ranks)):
                r = ranks[ri]
                order = {n: float(i) for i, n in enumerate(current[ranks[ri-1]])}
                bc = {n: _weighted_barycenter(G, n, order, "down")
                      for n in current[r]}
                current[r] = sorted(current[r], key=lambda n: bc[n])
                current[r] = _swap_reduce(G, current[ranks[ri-1]], current[r],
                                          fixed_is_upper=True)
        else:
            for ri in range(len(ranks) - 2, -1, -1):
                r = ranks[ri]
                order = {n: float(i) for i, n in enumerate(current[ranks[ri+1]])}
                bc = {n: _weighted_barycenter(G, n, order, "up")
                      for n in current[r]}
                current[r] = sorted(current[r], key=lambda n: bc[n])
                current[r] = _swap_reduce(G, current[ranks[ri+1]], current[r],
                                          fixed_is_upper=False)

        total = _total_crossings(G, current)
        if total < best_xings:
            best = {r: list(current[r]) for r in ranks}
            best_xings = total

    return best


# ─────────────────────────────────────────────────────────────────────
# Phase 4: Coordinate assignment
# ─────────────────────────────────────────────────────────────────────

def _weighted_median(positions_weights: list[tuple[float, float]]) -> float:
    """Weighted median of (position, weight) pairs."""
    if not positions_weights:
        return 0.0
    pw = sorted(positions_weights)
    total_w = sum(w for _, w in pw)
    if total_w == 0:
        return pw[len(pw) // 2][0]
    cum = 0.0
    for pos, w in pw:
        cum += w
        if cum >= total_w / 2:
            return pos
    return pw[-1][0]


def _subtree_size(G: nx.DiGraph, node: str, cache: dict[str, int]) -> int:
    """Count of all descendants (including self)."""
    if node in cache:
        return cache[node]
    children = list(G.successors(node))
    size = 1 + sum(_subtree_size(G, c, cache) for c in children)
    cache[node] = size
    return size


def assign_coordinates(G: nx.DiGraph, layers_by_rank: dict[int, list[str]],
                       node_spacing: float = 180,
                       layer_spacing: float = 220) -> dict[str, tuple[float, float]]:
    """
    Assign (x, y) coordinates using subtree-aware spacing and
    iterative weighted-median refinement.

    Nodes with larger subtrees get more space around them; leaf nodes
    are packed tighter.  This prevents a node with 7 children from
    being allocated the same space as a childless node.
    """
    ranks = sorted(layers_by_rank.keys())
    pos: dict[str, tuple[float, float]] = {}

    # Compute subtree sizes for proportional spacing
    st_cache: dict[str, int] = {}
    for n in G.nodes:
        _subtree_size(G, n, st_cache)

    # Initial placement: proportional spacing based on subtree size
    for r in ranks:
        nodes = layers_by_rank[r]
        # Each node gets space proportional to max(1, subtree_size - 1)
        weights = [max(1, st_cache.get(n, 1) - 1) for n in nodes]
        total_weight = sum(weights)
        # Total span for this layer
        total_span = total_weight * node_spacing * 0.7
        # Compute cumulative positions
        cum = 0.0
        positions = []
        for w in weights:
            positions.append(cum + w / 2)
            cum += w
        # Centre and scale
        mid = cum / 2
        for i, n in enumerate(nodes):
            y = (positions[i] - mid) * node_spacing * 0.7
            pos[n] = (r * layer_spacing, y)

    # Iterative refinement: 12 down+up passes
    for iteration in range(12):
        for ri in range(1, len(ranks)):
            _refine_layer_weighted(G, layers_by_rank[ranks[ri]], pos,
                                   ranks[ri] * layer_spacing, node_spacing,
                                   "down")
        for ri in range(len(ranks) - 2, -1, -1):
            _refine_layer_weighted(G, layers_by_rank[ranks[ri]], pos,
                                   ranks[ri] * layer_spacing, node_spacing,
                                   "up")

    return pos


def _refine_layer_weighted(G: nx.DiGraph, nodes: list[str],
                           pos: dict[str, tuple[float, float]],
                           x: float, min_spacing: float, direction: str):
    """
    Refine Y positions for one layer using weighted median of
    connected neighbours.
    """
    desired = {}
    for n in nodes:
        pw = []  # (position, weight) pairs
        if direction == "down":
            for p in G.predecessors(n):
                if p in pos:
                    pw.append((pos[p][1], _edge_weight(G, p, n)))
        else:
            for c in G.successors(n):
                if c in pos:
                    pw.append((pos[c][1], _edge_weight(G, n, c)))
            for p in G.predecessors(n):
                if p in pos:
                    pw.append((pos[p][1], _edge_weight(G, p, n)))

        if pw:
            desired[n] = _weighted_median(pw)
        else:
            desired[n] = pos[n][1]

    # Place in order, resolving overlaps
    ordered = sorted(nodes, key=lambda n: desired.get(n, pos[n][1]))
    _compact_place(ordered, desired, pos, x, min_spacing)


def _compact_place(ordered: list[str], desired: dict[str, float],
                   pos: dict[str, tuple[float, float]],
                   x: float, min_spacing: float):
    """Place nodes at their desired Y, compacting to resolve overlaps."""
    if not ordered:
        return

    # Forward pass: enforce minimum spacing
    placed = []
    for n in ordered:
        want = desired.get(n, 0)
        if placed:
            want = max(want, placed[-1] + min_spacing)
        placed.append(want)

    # Backward pass: pull nodes up toward their desired positions
    for i in range(len(ordered) - 2, -1, -1):
        upper = placed[i + 1] - min_spacing
        ideal = desired.get(ordered[i], placed[i])
        new_y = min(upper, max(ideal, placed[i]))
        if i > 0:
            new_y = max(new_y, placed[i - 1] + min_spacing)
        placed[i] = new_y

    # Apply
    for i, n in enumerate(ordered):
        pos[n] = (x, placed[i])


# ─────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────

def hierarchical_layout(G: nx.DiGraph,
                        node_spacing: float = 180,
                        layer_spacing: float = 220,
                        sweeps: int = 24) -> dict[str, dict[str, Any]]:
    """Full Sugiyama layout pipeline with edge-type weighting."""
    # 1. Layer assignment
    layers = assign_layers(G)

    # 2. Insert dummy nodes
    G2, layers2 = insert_dummy_nodes(G, layers)

    # 3. Build layer lists with sibling-aware initial ordering
    layers_by_rank: dict[int, list[str]] = defaultdict(list)
    for n, r in layers2.items():
        layers_by_rank[r].append(n)

    ranks = sorted(layers_by_rank.keys())
    for ri in range(1, len(ranks)):
        r = ranks[ri]
        parent_pos = {n: i for i, n in enumerate(layers_by_rank[ranks[ri-1]])}

        def sort_key(node, _pp=parent_pos):
            parents = [(p, _edge_weight(G2, p, node))
                       for p in G2.predecessors(node) if p in _pp]
            if not parents:
                return (float('inf'),)
            total_w = sum(w for _, w in parents)
            if total_w == 0:
                return (float('inf'),)
            wavg = sum(_pp[p] * w for p, w in parents) / total_w
            return (wavg,)

        layers_by_rank[r] = sorted(layers_by_rank[r], key=sort_key)

    # 4. Crossing minimisation
    layers_by_rank = minimise_crossings(G2, layers_by_rank, sweeps=sweeps)

    # 5. Coordinate assignment
    pos = assign_coordinates(G2, layers_by_rank,
                             node_spacing=node_spacing,
                             layer_spacing=layer_spacing)

    # Separate real vs dummy
    real_pos, dummy_pos = {}, {}
    for n, (x, y) in pos.items():
        target = dummy_pos if G2.nodes[n].get("is_dummy") else real_pos
        target[n] = {"x": x, "y": y}

    # Build edge routes through dummy waypoints
    edge_routes = {}
    for u, v, data in G.edges(data=True):
        route = []
        if u in real_pos:
            route.append(real_pos[u])
        for d in _find_dummy_chain(G2, u, v, layers2):
            if d in dummy_pos:
                route.append(dummy_pos[d])
        if v in real_pos:
            route.append(real_pos[v])
        edge_routes[(u, v)] = {"points": route, "data": data}

    return {"nodes": real_pos, "edges": edge_routes, "dummy_positions": dummy_pos}


def _find_dummy_chain(G2: nx.DiGraph, orig_u: str, orig_v: str,
                      layers: dict[str, int]) -> list[str]:
    """Find the chain of dummy nodes inserted between orig_u and orig_v."""
    span = layers[orig_v] - layers[orig_u]
    if span <= 1:
        return []
    chain = []
    current = orig_u
    for _ in range(span - 1):
        for succ in G2.successors(current):
            if G2.nodes[succ].get("is_dummy") and layers[succ] == layers[current] + 1:
                chain.append(succ)
                current = succ
                break
        else:
            break
    return chain
