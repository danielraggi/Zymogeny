"""
Sugiyama-style hierarchical DAG layout.

Addresses three common problems with naive hierarchical layouts:
1. Children of a single node getting spread across a layer
   → Priority grouping: siblings are placed contiguously
2. Avoidable edge crossings
   → Barycenter ordering with multiple up/down sweeps
3. Same-level nodes clumped regardless of parent position
   → Coordinate assignment uses median of connected neighbours,
     not a uniform grid

Algorithm phases:
  1. Layer assignment (longest-path from roots)
  2. Dummy node insertion for long edges
  3. Crossing minimisation (barycenter + adjacent swap)
  4. Coordinate assignment (priority-based placement)
"""

from __future__ import annotations

import math
from collections import defaultdict
from typing import Any

import networkx as nx


# ─────────────────────────────────────────────────────────────────────
# Phase 1: Layer assignment
# ─────────────────────────────────────────────────────────────────────

def assign_layers(G: nx.DiGraph) -> dict[str, int]:
    """Assign each node to a layer using longest-path-from-source."""
    layer = {}
    topo = list(nx.topological_sort(G))
    for n in topo:
        preds = list(G.predecessors(n))
        if not preds:
            layer[n] = 0
        else:
            layer[n] = max(layer[p] for p in preds) + 1
    return layer


# ─────────────────────────────────────────────────────────────────────
# Phase 2: Dummy nodes for long edges
# ─────────────────────────────────────────────────────────────────────

def insert_dummy_nodes(G: nx.DiGraph, layers: dict[str, int]):
    """
    For edges spanning more than one layer, insert dummy nodes so every
    edge connects adjacent layers.  Returns (G', layers') with dummies.
    """
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

def _barycenter(G: nx.DiGraph, node: str, layer_order: dict[str, int],
                direction: str) -> float:
    """Compute barycenter of a node's neighbours in the adjacent layer."""
    if direction == "down":
        neighbours = [n for n in G.predecessors(node) if n in layer_order]
    else:
        neighbours = [n for n in G.successors(node) if n in layer_order]
    if not neighbours:
        return layer_order.get(node, 0)
    return sum(layer_order[n] for n in neighbours) / len(neighbours)


def _count_crossings(G: nx.DiGraph, layer_a: list[str], layer_b: list[str]) -> int:
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
            a1, b1 = edges[i]
            a2, b2 = edges[j]
            if (a1 - a2) * (b1 - b2) < 0:
                crossings += 1
    return crossings


def _adjacent_swap(G: nx.DiGraph, fixed_layer: list[str],
                   free_layer: list[str]) -> list[str]:
    """Greedily swap adjacent pairs in free_layer to reduce crossings."""
    improved = True
    while improved:
        improved = False
        for i in range(len(free_layer) - 1):
            cur = _count_crossings(G, fixed_layer, free_layer)
            swapped = list(free_layer)
            swapped[i], swapped[i + 1] = swapped[i + 1], swapped[i]
            if _count_crossings(G, fixed_layer, swapped) < cur:
                free_layer = swapped
                improved = True
    return free_layer


def minimise_crossings(G: nx.DiGraph, layers_by_rank: dict[int, list[str]],
                       sweeps: int = 24) -> dict[int, list[str]]:
    """
    Barycenter heuristic with adjacent-swap refinement.
    Runs alternating down-sweeps and up-sweeps.
    """
    ranks = sorted(layers_by_rank.keys())
    best = {r: list(layers_by_rank[r]) for r in ranks}
    best_crossings = sum(
        _count_crossings(G, best[ranks[i]], best[ranks[i + 1]])
        for i in range(len(ranks) - 1)
    )

    current = {r: list(best[r]) for r in ranks}

    for sweep in range(sweeps):
        if sweep % 2 == 0:
            # Down sweep: fix upper layers, reorder lower
            for ri in range(1, len(ranks)):
                r = ranks[ri]
                order_map = {n: i for i, n in enumerate(current[ranks[ri - 1]])}
                bc = {}
                for n in current[r]:
                    bc[n] = _barycenter(G, n, order_map, "down")
                current[r] = sorted(current[r], key=lambda n: bc[n])
                current[r] = _adjacent_swap(G, current[ranks[ri - 1]], current[r])
        else:
            # Up sweep: fix lower layers, reorder upper
            for ri in range(len(ranks) - 2, -1, -1):
                r = ranks[ri]
                order_map = {n: i for i, n in enumerate(current[ranks[ri + 1]])}
                bc = {}
                for n in current[r]:
                    bc[n] = _barycenter(G, n, order_map, "up")
                current[r] = sorted(current[r], key=lambda n: bc[n])
                # Swap adjacent nodes in current[r] to reduce crossings
                # with the layer below (ranks[ri+1])
                current[r] = _adjacent_swap_free(
                    G, current[ranks[ri + 1]], current[r])

        total = sum(
            _count_crossings(G, current[ranks[i]], current[ranks[i + 1]])
            for i in range(len(ranks) - 1)
        )
        if total < best_crossings:
            best = {r: list(current[r]) for r in ranks}
            best_crossings = total

    return best


def _adjacent_swap_free(G: nx.DiGraph, lower_layer: list[str],
                        upper_layer: list[str]) -> list[str]:
    """Swap adjacent nodes in upper_layer to reduce crossings with lower_layer."""
    improved = True
    result = list(upper_layer)
    while improved:
        improved = False
        for i in range(len(result) - 1):
            cur = _count_crossings(G, result, lower_layer)
            swapped = list(result)
            swapped[i], swapped[i + 1] = swapped[i + 1], swapped[i]
            if _count_crossings(G, swapped, lower_layer) < cur:
                result = swapped
                improved = True
    return result


# ─────────────────────────────────────────────────────────────────────
# Phase 4: Coordinate assignment
# ─────────────────────────────────────────────────────────────────────

def assign_coordinates(G: nx.DiGraph, layers_by_rank: dict[int, list[str]],
                       node_spacing: float = 180,
                       layer_spacing: float = 220) -> dict[str, tuple[float, float]]:
    """
    Assign (x, y) coordinates.  X follows the layer (left-to-right),
    Y is the position within a layer.

    Uses a priority method: place each node at the median position of its
    parents (for a down pass) then shift to resolve overlaps, keeping
    siblings grouped together.
    """
    ranks = sorted(layers_by_rank.keys())
    pos: dict[str, tuple[float, float]] = {}

    # Initial placement: centre each layer
    for r in ranks:
        nodes = layers_by_rank[r]
        for i, n in enumerate(nodes):
            x = r * layer_spacing
            y = i * node_spacing
            pos[n] = (x, y)

    # Down pass: pull nodes toward parent median
    for ri in range(1, len(ranks)):
        r = ranks[ri]
        nodes = layers_by_rank[r]
        desired = {}
        for n in nodes:
            parents = [p for p in G.predecessors(n) if p in pos]
            if parents:
                parent_ys = sorted(pos[p][1] for p in parents)
                median_y = parent_ys[len(parent_ys) // 2]
                desired[n] = median_y
            else:
                desired[n] = pos[n][1]

        # Sort by desired position, then assign with minimum spacing
        ordered = sorted(nodes, key=lambda n: desired[n])
        # Re-assign positions respecting minimum spacing
        _compact_place(ordered, desired, pos, r * layer_spacing, node_spacing)

    # Up pass: pull nodes toward child median (refine)
    for ri in range(len(ranks) - 2, -1, -1):
        r = ranks[ri]
        nodes = layers_by_rank[r]
        desired = {}
        for n in nodes:
            children = [c for c in G.successors(n) if c in pos]
            parents = [p for p in G.predecessors(n) if p in pos]
            connected = children + parents
            if connected:
                ys = sorted(pos[c][1] for c in connected)
                median_y = ys[len(ys) // 2]
                desired[n] = median_y
            else:
                desired[n] = pos[n][1]

        ordered = sorted(nodes, key=lambda n: desired[n])
        _compact_place(ordered, desired, pos, r * layer_spacing, node_spacing)

    # Final centering pass: shift each layer so its centre aligns better
    # with the overall graph centre
    all_ys = [pos[n][1] for n in pos]
    if all_ys:
        global_mid = (min(all_ys) + max(all_ys)) / 2
        for r in ranks:
            nodes = layers_by_rank[r]
            if not nodes:
                continue
            layer_ys = [pos[n][1] for n in nodes]
            layer_mid = (min(layer_ys) + max(layer_ys)) / 2
            # Pull gently toward global centre (30% weight)
            shift = (global_mid - layer_mid) * 0.3
            for n in nodes:
                x, y = pos[n]
                pos[n] = (x, y + shift)

    return pos


def _compact_place(ordered: list[str], desired: dict[str, float],
                   pos: dict[str, tuple[float, float]],
                   x: float, min_spacing: float):
    """Place nodes at their desired Y, compacting to resolve overlaps."""
    if not ordered:
        return
    # First pass: place greedily downward
    placed_y = []
    for n in ordered:
        want = desired[n]
        if placed_y:
            want = max(want, placed_y[-1] + min_spacing)
        placed_y.append(want)
        pos[n] = (x, want)

    # Second pass: shift upward where possible to reduce deviation
    for i in range(len(ordered) - 2, -1, -1):
        n = ordered[i]
        upper_bound = placed_y[i + 1] - min_spacing
        ideal = desired[n]
        new_y = min(upper_bound, max(ideal, placed_y[i]))
        if i > 0:
            new_y = max(new_y, placed_y[i - 1] + min_spacing)
        placed_y[i] = new_y
        pos[n] = (x, new_y)


# ─────────────────────────────────────────────────────────────────────
# Public API
# ─────────────────────────────────────────────────────────────────────

def hierarchical_layout(G: nx.DiGraph,
                        node_spacing: float = 180,
                        layer_spacing: float = 220,
                        sweeps: int = 24) -> dict[str, dict[str, Any]]:
    """
    Full Sugiyama layout pipeline.  Returns dict mapping node id → {x, y}.
    Dummy nodes are removed from the final output but their positions
    are returned separately for edge routing.
    """
    # 1. Layer assignment
    layers = assign_layers(G)

    # 2. Insert dummy nodes
    G2, layers2 = insert_dummy_nodes(G, layers)

    # 3. Build layer lists
    layers_by_rank: dict[int, list[str]] = defaultdict(list)
    for n, r in layers2.items():
        layers_by_rank[r].append(n)

    # Initial ordering: sort by parent barycenter within each layer
    ranks = sorted(layers_by_rank.keys())
    for ri in range(1, len(ranks)):
        r = ranks[ri]
        parent_rank = ranks[ri - 1]
        parent_pos = {n: i for i, n in enumerate(layers_by_rank[parent_rank])}
        def sort_key(node):
            parents = [p for p in G2.predecessors(node) if p in parent_pos]
            if not parents:
                return 0
            return sum(parent_pos[p] for p in parents) / len(parents)
        layers_by_rank[r] = sorted(layers_by_rank[r], key=sort_key)

    # 4. Crossing minimisation
    layers_by_rank = minimise_crossings(G2, layers_by_rank, sweeps=sweeps)

    # 5. Coordinate assignment
    pos = assign_coordinates(G2, layers_by_rank,
                             node_spacing=node_spacing,
                             layer_spacing=layer_spacing)

    # Separate real nodes from dummies
    real_pos = {}
    dummy_pos = {}
    for n, (x, y) in pos.items():
        if G2.nodes[n].get("is_dummy"):
            dummy_pos[n] = {"x": x, "y": y}
        else:
            real_pos[n] = {"x": x, "y": y}

    # Build edge routes (through dummy waypoints)
    edge_routes = {}
    for u, v, data in G.edges(data=True):
        route = [real_pos[u]] if u in real_pos else []
        # Find dummy chain from u to v in G2
        chain = _find_dummy_chain(G2, u, v, layers2)
        for d in chain:
            if d in dummy_pos:
                route.append(dummy_pos[d])
        if v in real_pos:
            route.append(real_pos[v])
        edge_routes[(u, v)] = {"points": route, "data": data}

    return {"nodes": real_pos, "edges": edge_routes, "dummy_positions": dummy_pos}


def _find_dummy_chain(G2: nx.DiGraph, orig_u: str, orig_v: str,
                      layers: dict[str, int]) -> list[str]:
    """Find the chain of dummy nodes between orig_u and orig_v."""
    span = layers[orig_v] - layers[orig_u]
    if span <= 1:
        return []
    chain = []
    current = orig_u
    for _ in range(span - 1):
        found = False
        for succ in G2.successors(current):
            if G2.nodes[succ].get("is_dummy") and layers[succ] == layers[current] + 1:
                chain.append(succ)
                current = succ
                found = True
                break
        if not found:
            break
    return chain
