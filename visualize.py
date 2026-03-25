#!/usr/bin/env python3
"""
Generate an interactive HTML visualization of the phylogenetic network.

Usage:
    python visualize.py            # writes phylogeny.html
    python visualize.py -o out.html
"""

import argparse
import json

from graph import build_network
from layout import hierarchical_layout

# ─────────────────────────────────────────────────────────────────────
# Colour palette by node type
# ─────────────────────────────────────────────────────────────────────
TYPE_COLOURS = {
    "species": "#2e8b57",       # sea green
    "hybrid": "#e67e22",        # orange
    "population": "#3498db",    # blue
    "strain": "#f39c12",        # amber
}


def build_vis_data(G, layout_result):
    """Convert graph + layout into JSON-serialisable dicts for the HTML."""
    nodes = []
    for nid, coord in layout_result["nodes"].items():
        data = G.nodes[nid]
        nodes.append({
            "id": nid,
            "x": coord["x"],
            "y": coord["y"],
            "label": data.get("display_name", nid),
            "type": data.get("type", "unknown"),
            "color": TYPE_COLOURS.get(data.get("type"), "#999"),
            "is_hybrid": data.get("is_hybrid", False),
            "fermentation": data.get("fermentation", []),
            "wild": data.get("wild", False),
        })

    edges = []
    for (u, v), info in layout_result["edges"].items():
        edge_data = info.get("data", {})
        edges.append({
            "source": u,
            "target": v,
            "points": info["points"],
            "type": edge_data.get("type", "divergence"),
            "admixture": edge_data.get("admixture_fraction"),
        })

    return {"nodes": nodes, "edges": edges}


HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Saccharomyces Phylogenetic Network</title>
<style>
* { margin: 0; padding: 0; box-sizing: border-box; }
body { background: #1a1a2e; font-family: 'Segoe UI', system-ui, sans-serif; overflow: hidden; }
#canvas { display: block; cursor: grab; }
#canvas.grabbing { cursor: grabbing; }

/* ─── Search bar ─── */
#search-container {
  position: fixed; top: 16px; left: 50%; transform: translateX(-50%);
  z-index: 10;
}
#search {
  width: 320px; padding: 10px 16px 10px 38px;
  background: rgba(30,30,50,0.94); color: #eee;
  border: 1px solid rgba(255,255,255,0.12); border-radius: 24px;
  font-size: 14px; outline: none;
  transition: border-color 0.2s, box-shadow 0.2s;
}
#search:focus {
  border-color: rgba(52,152,219,0.6);
  box-shadow: 0 0 12px rgba(52,152,219,0.2);
}
#search::placeholder { color: #666; }
#search-icon {
  position: absolute; left: 14px; top: 50%; transform: translateY(-50%);
  color: #666; font-size: 14px; pointer-events: none;
}
#search-results {
  position: absolute; top: 44px; left: 0; right: 0;
  background: rgba(25,25,45,0.97); border-radius: 12px;
  border: 1px solid rgba(255,255,255,0.08);
  max-height: 260px; overflow-y: auto; display: none;
}
.search-item {
  padding: 8px 16px; color: #ccc; cursor: pointer;
  font-size: 13px; border-bottom: 1px solid rgba(255,255,255,0.04);
}
.search-item:hover { background: rgba(52,152,219,0.15); color: #fff; }
.search-item .type-badge {
  display: inline-block; padding: 1px 7px; border-radius: 8px;
  font-size: 10px; margin-left: 8px; color: #fff; opacity: 0.85;
}

/* ─── Filter bar ─── */
#filter-bar {
  position: fixed; top: 60px; left: 50%; transform: translateX(-50%);
  z-index: 9; display: flex; gap: 6px; flex-wrap: wrap; justify-content: center;
}
.filter-chip {
  padding: 4px 12px; border-radius: 14px; font-size: 11px;
  background: rgba(30,30,50,0.9); color: #aaa; cursor: pointer;
  border: 1px solid rgba(255,255,255,0.08);
  transition: all 0.2s; user-select: none;
}
.filter-chip.active { color: #fff; border-color: rgba(255,255,255,0.25); }
.filter-chip:hover { background: rgba(50,50,80,0.9); }

/* ─── Legend ─── */
#legend {
  position: fixed; top: 16px; right: 16px;
  background: rgba(30,30,50,0.92); border-radius: 10px; padding: 14px 18px;
  color: #ccc; font-size: 13px; line-height: 1.9;
  border: 1px solid rgba(255,255,255,0.08);
  z-index: 10;
}
.legend-row {
  cursor: pointer; padding: 1px 4px; border-radius: 4px;
  transition: opacity 0.2s;
}
.legend-row:hover { background: rgba(255,255,255,0.05); }
.legend-row.dimmed { opacity: 0.35; }
.legend-dot {
  display: inline-block; width: 12px; height: 12px; border-radius: 50%;
  margin-right: 6px; vertical-align: middle;
}

/* ─── Tooltip ─── */
#tooltip {
  position: fixed; display: none; pointer-events: none;
  background: rgba(20,20,40,0.96); color: #eee; padding: 10px 14px;
  border-radius: 8px; font-size: 12px; line-height: 1.6;
  border: 1px solid rgba(255,255,255,0.1);
  max-width: 300px; z-index: 20;
  backdrop-filter: blur(8px);
}
.tooltip-label { font-size: 14px; font-weight: 600; margin-bottom: 4px; }
.tooltip-row { color: #aaa; }
.tooltip-row span { color: #ddd; }
.tooltip-tag {
  display: inline-block; padding: 1px 7px; border-radius: 8px;
  font-size: 10px; margin: 1px 2px; background: rgba(255,255,255,0.08);
  color: #ccc;
}

/* ─── Controls ─── */
#controls {
  position: fixed; bottom: 16px; left: 16px;
  color: #666; font-size: 11px; z-index: 10;
}

/* ─── Stats badge ─── */
#stats {
  position: fixed; bottom: 16px; right: 16px;
  color: #555; font-size: 11px; z-index: 10;
}

/* Scrollbar */
#search-results::-webkit-scrollbar { width: 5px; }
#search-results::-webkit-scrollbar-thumb { background: rgba(255,255,255,0.1); border-radius: 3px; }
</style>
</head>
<body>

<div id="search-container">
  <span id="search-icon">&#128269;</span>
  <input id="search" type="text" placeholder="Search strains, species, populations..." autocomplete="off">
  <div id="search-results"></div>
</div>

<div id="filter-bar"></div>

<canvas id="canvas"></canvas>

<div id="legend">
  <strong style="color:#fff">Node types</strong> <span style="font-size:10px;color:#666">(click to toggle)</span><br>
  <div class="legend-row" data-type="species"><span class="legend-dot" style="background:#2e8b57"></span>Species</div>
  <div class="legend-row" data-type="hybrid"><span class="legend-dot" style="background:#e67e22"></span>Hybrid</div>
  <div class="legend-row" data-type="population"><span class="legend-dot" style="background:#3498db"></span>Population</div>
  <div class="legend-row" data-type="strain"><span class="legend-dot" style="background:#f39c12"></span>Strain</div>
  <br>
  <strong style="color:#fff">Edge types</strong><br>
  <span style="color:#667">━━</span> Divergence<br>
  <span style="color:#e67e22">╌╌</span> Hybridisation<br>
  <span style="color:#9b59b6">┈┈</span> Introgression
</div>

<div id="tooltip"></div>
<div id="controls">Scroll to zoom · Drag to pan · Hover for details · Click node to focus</div>
<div id="stats"></div>

<script>
const DATA = __DATA__;

const canvas = document.getElementById('canvas');
const ctx = canvas.getContext('2d');
const tooltip = document.getElementById('tooltip');
const searchInput = document.getElementById('search');
const searchResults = document.getElementById('search-results');
const filterBar = document.getElementById('filter-bar');
const statsEl = document.getElementById('stats');

let W, H, dpr;
let camX = 0, camY = 0, zoom = 1;
let dragging = false, dragStartX, dragStartY, camStartX, camStartY;
let hoveredNode = null;
let focusedNode = null;       // click-to-focus
let searchHighlight = null;   // search highlight
let hiddenTypes = new Set();  // legend toggle
let activeFermentation = null; // fermentation filter

// Build edge lookup: node id → connected edges
const edgesByNode = {};
for (const n of DATA.nodes) edgesByNode[n.id] = [];
for (const e of DATA.edges) {
  if (edgesByNode[e.source]) edgesByNode[e.source].push(e);
  if (edgesByNode[e.target]) edgesByNode[e.target].push(e);
}

// Build adjacency: node id → set of connected node ids
const adjacency = {};
for (const n of DATA.nodes) adjacency[n.id] = new Set();
for (const e of DATA.edges) {
  if (adjacency[e.source]) adjacency[e.source].add(e.target);
  if (adjacency[e.target]) adjacency[e.target].add(e.source);
}

// ─── Collect all fermentation types for filter chips ───
const allFermentations = new Set();
for (const n of DATA.nodes) {
  for (const f of (n.fermentation || [])) allFermentations.add(f);
}
const fermList = [...allFermentations].sort();
for (const f of fermList) {
  const chip = document.createElement('div');
  chip.className = 'filter-chip';
  chip.textContent = f;
  chip.addEventListener('click', () => {
    if (activeFermentation === f) {
      activeFermentation = null;
      chip.classList.remove('active');
    } else {
      activeFermentation = f;
      document.querySelectorAll('.filter-chip').forEach(c => c.classList.remove('active'));
      chip.classList.add('active');
    }
    focusedNode = null;
    draw();
  });
  filterBar.appendChild(chip);
}

// ─── Stats ───
const typeCounts = {};
for (const n of DATA.nodes) typeCounts[n.type] = (typeCounts[n.type] || 0) + 1;
statsEl.textContent = `${DATA.nodes.length} nodes · ${DATA.edges.length} edges`;

// ─── Visibility helpers ───
function isNodeVisible(n) {
  if (hiddenTypes.has(n.type)) return false;
  return true;
}

function isNodeHighlighted(n) {
  // If fermentation filter is active, check membership
  if (activeFermentation) {
    return (n.fermentation || []).includes(activeFermentation);
  }
  // If a node is focused, only it and neighbors are highlighted
  if (focusedNode) {
    return n.id === focusedNode.id || adjacency[focusedNode.id].has(n.id);
  }
  // Search highlight
  if (searchHighlight) {
    return n.id === searchHighlight.id;
  }
  return true; // all highlighted when no filter active
}

function isEdgeHighlighted(e) {
  if (activeFermentation) {
    const sNode = DATA.nodes.find(n => n.id === e.source);
    const tNode = DATA.nodes.find(n => n.id === e.target);
    const sMatch = sNode && (sNode.fermentation || []).includes(activeFermentation);
    const tMatch = tNode && (tNode.fermentation || []).includes(activeFermentation);
    return sMatch || tMatch;
  }
  if (focusedNode) {
    return e.source === focusedNode.id || e.target === focusedNode.id;
  }
  return true;
}

// ─── Resize ───
function resize() {
  dpr = window.devicePixelRatio || 1;
  W = window.innerWidth; H = window.innerHeight;
  canvas.width = W * dpr; canvas.height = H * dpr;
  canvas.style.width = W + 'px'; canvas.style.height = H + 'px';
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  draw();
}
window.addEventListener('resize', resize);

function centreView() {
  if (!DATA.nodes.length) return;
  let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
  for (const n of DATA.nodes) {
    minX = Math.min(minX, n.x); maxX = Math.max(maxX, n.x);
    minY = Math.min(minY, n.y); maxY = Math.max(maxY, n.y);
  }
  const gw = maxX - minX + 200, gh = maxY - minY + 200;
  zoom = Math.min(W / gw, H / gh, 1.5);
  camX = (minX + maxX) / 2 - W / (2 * zoom);
  camY = (minY + maxY) / 2 - H / (2 * zoom);
}

function screenToWorld(sx, sy) {
  return { x: sx / zoom + camX, y: sy / zoom + camY };
}
function worldToScreen(wx, wy) {
  return { x: (wx - camX) * zoom, y: (wy - camY) * zoom };
}

// Smooth pan-to-node animation
function panToNode(n, cb) {
  const targetCamX = n.x - W / (2 * zoom);
  const targetCamY = n.y - H / (2 * zoom);
  const startX = camX, startY = camY;
  const duration = 400;
  const start = performance.now();
  function step(t) {
    const elapsed = t - start;
    const p = Math.min(elapsed / duration, 1);
    const ease = 1 - Math.pow(1 - p, 3); // ease-out cubic
    camX = startX + (targetCamX - startX) * ease;
    camY = startY + (targetCamY - startY) * ease;
    draw();
    if (p < 1) requestAnimationFrame(step);
    else if (cb) cb();
  }
  requestAnimationFrame(step);
}

// ─── Drawing ───────────────────────────────────────────────────
const EDGE_COLORS = {
  divergence: 'rgba(100,120,140,__A__)',
  hybridisation: 'rgba(230,126,34,__A__)',
  introgression: 'rgba(155,89,182,__A__)',
};

function edgeColor(type, alpha) {
  const base = {
    divergence: [100,120,140],
    hybridisation: [230,126,34],
    introgression: [155,89,182],
  };
  const c = base[type] || base.divergence;
  return `rgba(${c[0]},${c[1]},${c[2]},${alpha})`;
}

function drawEdge(e) {
  const pts = e.points;
  if (pts.length < 2) return;

  // Check if source/target types are hidden
  const sNode = DATA.nodes.find(n => n.id === e.source);
  const tNode = DATA.nodes.find(n => n.id === e.target);
  if (sNode && hiddenTypes.has(sNode.type)) return;
  if (tNode && hiddenTypes.has(tNode.type)) return;

  const highlighted = isEdgeHighlighted(e);
  const alpha = highlighted ? 0.55 : 0.08;

  ctx.beginPath();
  ctx.strokeStyle = edgeColor(e.type, alpha);
  ctx.lineWidth = (e.type === 'divergence' ? 1.5 : 2) / zoom;

  if (e.type === 'introgression') {
    ctx.setLineDash([6 / zoom, 4 / zoom]);
  } else if (e.type === 'hybridisation') {
    ctx.setLineDash([8 / zoom, 3 / zoom]);
  } else {
    ctx.setLineDash([]);
  }

  const s0 = worldToScreen(pts[0].x, pts[0].y);
  ctx.moveTo(s0.x, s0.y);

  if (pts.length === 2) {
    const s1 = worldToScreen(pts[1].x, pts[1].y);
    const mx = (s0.x + s1.x) / 2;
    ctx.bezierCurveTo(mx, s0.y, mx, s1.y, s1.x, s1.y);
  } else {
    for (let i = 1; i < pts.length; i++) {
      const si = worldToScreen(pts[i].x, pts[i].y);
      const prev = worldToScreen(pts[i-1].x, pts[i-1].y);
      const mx = (prev.x + si.x) / 2;
      ctx.bezierCurveTo(mx, prev.y, mx, si.y, si.x, si.y);
    }
  }
  ctx.stroke();
  ctx.setLineDash([]);

  // Arrowhead
  if (highlighted) {
    const last = worldToScreen(pts[pts.length - 1].x, pts[pts.length - 1].y);
    const prev = worldToScreen(pts[pts.length - 2].x, pts[pts.length - 2].y);
    const angle = Math.atan2(last.y - prev.y, last.x - prev.x);
    const aw = 8;
    ctx.beginPath();
    ctx.fillStyle = edgeColor(e.type, alpha);
    ctx.moveTo(last.x, last.y);
    ctx.lineTo(last.x - aw * Math.cos(angle - 0.4), last.y - aw * Math.sin(angle - 0.4));
    ctx.lineTo(last.x - aw * Math.cos(angle + 0.4), last.y - aw * Math.sin(angle + 0.4));
    ctx.fill();
  }
}

function drawNode(n) {
  if (!isNodeVisible(n)) return;

  const s = worldToScreen(n.x, n.y);
  const r = (n.type === 'species' ? 14 : n.type === 'population' ? 11 : 8);
  const highlighted = isNodeHighlighted(n);
  const isSearchTarget = searchHighlight && n.id === searchHighlight.id;
  const isFocused = focusedNode && n.id === focusedNode.id;
  const nodeAlpha = highlighted ? 1.0 : 0.12;

  // Shadow for species/hybrid nodes
  if (highlighted && (n.type === 'species' || n.type === 'hybrid')) {
    ctx.save();
    ctx.shadowColor = n.color;
    ctx.shadowBlur = 12;
    ctx.beginPath();
    ctx.arc(s.x, s.y, r, 0, Math.PI * 2);
    ctx.fillStyle = 'transparent';
    ctx.fill();
    ctx.restore();
  }

  // Search pulse ring
  if (isSearchTarget) {
    ctx.beginPath();
    ctx.arc(s.x, s.y, r + 10, 0, Math.PI * 2);
    ctx.strokeStyle = '#3498db';
    ctx.lineWidth = 2;
    ctx.stroke();
  }

  // Focus ring
  if (isFocused) {
    ctx.beginPath();
    ctx.arc(s.x, s.y, r + 8, 0, Math.PI * 2);
    ctx.strokeStyle = 'rgba(255,255,255,0.4)';
    ctx.lineWidth = 2;
    ctx.stroke();
  }

  // Glow for hovered
  if (hoveredNode === n) {
    ctx.beginPath();
    ctx.arc(s.x, s.y, r + 6, 0, Math.PI * 2);
    ctx.fillStyle = n.color + '33';
    ctx.fill();
  }

  // Node circle
  ctx.globalAlpha = nodeAlpha;
  ctx.beginPath();
  ctx.arc(s.x, s.y, r, 0, Math.PI * 2);
  ctx.fillStyle = n.color;
  ctx.fill();
  ctx.strokeStyle = 'rgba(255,255,255,0.25)';
  ctx.lineWidth = 1;
  ctx.stroke();

  // Label
  const fontSize = n.type === 'strain' ? 10 : 12;
  ctx.font = (n.type === 'species' ? 'bold ' : '') + `${fontSize}px 'Segoe UI', system-ui, sans-serif`;
  ctx.textAlign = 'center';
  ctx.fillStyle = highlighted ? '#ddd' : '#555';
  ctx.fillText(n.label, s.x, s.y + r + fontSize + 3);
  ctx.globalAlpha = 1.0;
}

function draw() {
  ctx.clearRect(0, 0, W, H);
  for (const e of DATA.edges) drawEdge(e);
  for (const n of DATA.nodes) drawNode(n);
}

// ─── Interaction ───────────────────────────────────────────────
canvas.addEventListener('wheel', (e) => {
  e.preventDefault();
  const factor = e.deltaY > 0 ? 0.9 : 1.1;
  const wx = e.offsetX / zoom + camX;
  const wy = e.offsetY / zoom + camY;
  zoom *= factor;
  zoom = Math.max(0.1, Math.min(zoom, 8));
  camX = wx - e.offsetX / zoom;
  camY = wy - e.offsetY / zoom;
  draw();
}, { passive: false });

let dragMoved = false;
canvas.addEventListener('mousedown', (e) => {
  dragging = true;
  dragMoved = false;
  dragStartX = e.clientX; dragStartY = e.clientY;
  camStartX = camX; camStartY = camY;
  canvas.classList.add('grabbing');
});

canvas.addEventListener('mousemove', (e) => {
  if (dragging) {
    const dx = e.clientX - dragStartX, dy = e.clientY - dragStartY;
    if (Math.abs(dx) > 3 || Math.abs(dy) > 3) dragMoved = true;
    camX = camStartX - dx / zoom;
    camY = camStartY - dy / zoom;
    draw();
  }
  // Hover detection
  const w = screenToWorld(e.offsetX, e.offsetY);
  let found = null;
  for (const n of DATA.nodes) {
    if (!isNodeVisible(n)) continue;
    const dx = n.x - w.x, dy = n.y - w.y;
    if (dx * dx + dy * dy < 400) { found = n; break; }
  }
  if (found !== hoveredNode) {
    hoveredNode = found;
    draw();
    if (found) {
      tooltip.style.display = 'block';
      updateTooltip(found, e.clientX, e.clientY);
    } else {
      tooltip.style.display = 'none';
    }
  } else if (found) {
    tooltip.style.left = (e.clientX + 14) + 'px';
    tooltip.style.top = (e.clientY - 10) + 'px';
  }
});

canvas.addEventListener('mouseup', (e) => {
  dragging = false;
  canvas.classList.remove('grabbing');
  // Click-to-focus (only if we didn't drag)
  if (!dragMoved) {
    const w = screenToWorld(e.offsetX, e.offsetY);
    let clicked = null;
    for (const n of DATA.nodes) {
      if (!isNodeVisible(n)) continue;
      const dx = n.x - w.x, dy = n.y - w.y;
      if (dx * dx + dy * dy < 400) { clicked = n; break; }
    }
    if (clicked) {
      focusedNode = (focusedNode === clicked) ? null : clicked;
      searchHighlight = null;
    } else {
      focusedNode = null;
    }
    draw();
  }
});

function updateTooltip(n, mx, my) {
  let html = `<div class="tooltip-label" style="color:${n.color}">${n.label}</div>`;
  const typeLabel = n.type.charAt(0).toUpperCase() + n.type.slice(1);
  html += `<div class="tooltip-row">Type: <span>${typeLabel}</span></div>`;
  if (n.is_hybrid) html += `<div class="tooltip-row">Hybrid: <span>Yes</span></div>`;
  if (n.wild) html += `<div class="tooltip-row">Origin: <span>Wild</span></div>`;
  else if (n.type === 'species' || n.type === 'population')
    html += `<div class="tooltip-row">Origin: <span>Domesticated</span></div>`;
  if (n.fermentation && n.fermentation.length > 0) {
    html += `<div class="tooltip-row" style="margin-top:4px">Fermentation:</div>`;
    html += n.fermentation.map(f => `<span class="tooltip-tag">${f}</span>`).join('');
  }
  // Show number of connections
  const connCount = (adjacency[n.id] || new Set()).size;
  html += `<div class="tooltip-row" style="margin-top:4px;color:#666">${connCount} connection${connCount !== 1 ? 's' : ''}</div>`;
  tooltip.innerHTML = html;
  tooltip.style.left = (mx + 14) + 'px';
  tooltip.style.top = (my - 10) + 'px';
}

// ─── Legend toggle ───
document.querySelectorAll('.legend-row').forEach(row => {
  row.addEventListener('click', () => {
    const type = row.dataset.type;
    if (hiddenTypes.has(type)) {
      hiddenTypes.delete(type);
      row.classList.remove('dimmed');
    } else {
      hiddenTypes.add(type);
      row.classList.add('dimmed');
    }
    draw();
  });
});

// ─── Search ───
searchInput.addEventListener('input', () => {
  const q = searchInput.value.trim().toLowerCase();
  if (!q) { searchResults.style.display = 'none'; return; }
  const matches = DATA.nodes
    .filter(n => n.label.toLowerCase().includes(q) || n.id.toLowerCase().includes(q))
    .slice(0, 12);
  if (matches.length === 0) { searchResults.style.display = 'none'; return; }
  const typeColors = { species: '#2e8b57', hybrid: '#e67e22', population: '#3498db', strain: '#f39c12' };
  searchResults.innerHTML = matches.map(n =>
    `<div class="search-item" data-id="${n.id}">${n.label}` +
    `<span class="type-badge" style="background:${typeColors[n.type] || '#666'}">${n.type}</span></div>`
  ).join('');
  searchResults.style.display = 'block';
  // Attach click handlers
  searchResults.querySelectorAll('.search-item').forEach(item => {
    item.addEventListener('click', () => {
      const node = DATA.nodes.find(n => n.id === item.dataset.id);
      if (node) {
        searchHighlight = node;
        focusedNode = null;
        panToNode(node);
        searchResults.style.display = 'none';
        searchInput.value = node.label;
        searchInput.blur();
      }
    });
  });
});

searchInput.addEventListener('keydown', (e) => {
  if (e.key === 'Escape') {
    searchInput.value = '';
    searchResults.style.display = 'none';
    searchHighlight = null;
    draw();
  }
  if (e.key === 'Enter') {
    const first = searchResults.querySelector('.search-item');
    if (first) first.click();
  }
});

// Close search results on outside click
document.addEventListener('click', (e) => {
  if (!e.target.closest('#search-container')) {
    searchResults.style.display = 'none';
  }
});

// Keyboard shortcut: / to focus search
document.addEventListener('keydown', (e) => {
  if (e.key === '/' && document.activeElement !== searchInput) {
    e.preventDefault();
    searchInput.focus();
  }
  if (e.key === 'Escape') {
    focusedNode = null;
    searchHighlight = null;
    activeFermentation = null;
    document.querySelectorAll('.filter-chip').forEach(c => c.classList.remove('active'));
    draw();
  }
});

// Init – resize first so W/H are defined for centreView
resize();
centreView();
draw();
</script>
</body>
</html>
"""


def generate_html(output_path: str = "phylogeny.html"):
    G = build_network()
    result = hierarchical_layout(G, node_spacing=160, layer_spacing=260, sweeps=30)
    vis_data = build_vis_data(G, result)
    html = HTML_TEMPLATE.replace("__DATA__", json.dumps(vis_data))
    with open(output_path, "w") as f:
        f.write(html)
    print(f"Wrote {output_path} ({len(vis_data['nodes'])} nodes, {len(vis_data['edges'])} edges)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--output", default="phylogeny.html")
    args = parser.parse_args()
    generate_html(args.output)
