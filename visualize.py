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
#canvas { display: block; }
#legend {
  position: fixed; top: 16px; right: 16px;
  background: rgba(30,30,50,0.92); border-radius: 8px; padding: 14px 18px;
  color: #ccc; font-size: 13px; line-height: 1.8;
  border: 1px solid rgba(255,255,255,0.08);
}
.legend-dot {
  display: inline-block; width: 12px; height: 12px; border-radius: 50%;
  margin-right: 6px; vertical-align: middle;
}
#tooltip {
  position: fixed; display: none; pointer-events: none;
  background: rgba(20,20,40,0.95); color: #eee; padding: 8px 12px;
  border-radius: 6px; font-size: 12px; border: 1px solid rgba(255,255,255,0.1);
  max-width: 260px;
}
#controls {
  position: fixed; bottom: 16px; left: 16px;
  color: #888; font-size: 12px;
}
</style>
</head>
<body>
<canvas id="canvas"></canvas>
<div id="legend">
  <strong style="color:#fff">Node types</strong><br>
  <span class="legend-dot" style="background:#2e8b57"></span>Species<br>
  <span class="legend-dot" style="background:#e67e22"></span>Hybrid<br>
  <span class="legend-dot" style="background:#3498db"></span>Population<br>
  <span class="legend-dot" style="background:#f39c12"></span>Strain<br>
  <br>
  <strong style="color:#fff">Edge types</strong><br>
  <span style="color:#667">━━</span> Divergence<br>
  <span style="color:#e67e22">╌╌</span> Hybridisation<br>
  <span style="color:#9b59b6">┈┈</span> Introgression
</div>
<div id="tooltip"></div>
<div id="controls">Scroll to zoom · Drag to pan · Hover nodes for details</div>

<script>
const DATA = __DATA__;

const canvas = document.getElementById('canvas');
const ctx = canvas.getContext('2d');
const tooltip = document.getElementById('tooltip');

let W, H, dpr;
let camX = 0, camY = 0, zoom = 1;
let dragging = false, dragStartX, dragStartY, camStartX, camStartY;
let hoveredNode = null;

function resize() {
  dpr = window.devicePixelRatio || 1;
  W = window.innerWidth; H = window.innerHeight;
  canvas.width = W * dpr; canvas.height = H * dpr;
  canvas.style.width = W + 'px'; canvas.style.height = H + 'px';
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  draw();
}
window.addEventListener('resize', resize);

// Centre the graph initially
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

// ─── Drawing ───────────────────────────────────────────────────
const EDGE_COLORS = {
  divergence: 'rgba(100,120,140,0.45)',
  hybridisation: 'rgba(230,126,34,0.55)',
  introgression: 'rgba(155,89,182,0.55)',
};

function drawEdge(e) {
  const pts = e.points;
  if (pts.length < 2) return;

  ctx.beginPath();
  ctx.strokeStyle = EDGE_COLORS[e.type] || EDGE_COLORS.divergence;
  ctx.lineWidth = (e.type === 'divergence' ? 1.5 : 2) / zoom;

  if (e.type === 'introgression') {
    ctx.setLineDash([6 / zoom, 4 / zoom]);
  } else if (e.type === 'hybridisation') {
    ctx.setLineDash([8 / zoom, 3 / zoom]);
  } else {
    ctx.setLineDash([]);
  }

  // Draw smooth curve through waypoints
  const s0 = worldToScreen(pts[0].x, pts[0].y);
  ctx.moveTo(s0.x, s0.y);

  if (pts.length === 2) {
    const s1 = worldToScreen(pts[1].x, pts[1].y);
    // Bezier with horizontal tangents for hierarchical look
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

  // Arrowhead at target
  const last = worldToScreen(pts[pts.length - 1].x, pts[pts.length - 1].y);
  const prev = worldToScreen(pts[pts.length - 2].x, pts[pts.length - 2].y);
  const angle = Math.atan2(last.y - prev.y, last.x - prev.x);
  const aw = 8, ah = 5;
  ctx.beginPath();
  ctx.fillStyle = ctx.strokeStyle;
  ctx.moveTo(last.x, last.y);
  ctx.lineTo(last.x - aw * Math.cos(angle - 0.4), last.y - aw * Math.sin(angle - 0.4));
  ctx.lineTo(last.x - aw * Math.cos(angle + 0.4), last.y - aw * Math.sin(angle + 0.4));
  ctx.fill();
}

function drawNode(n) {
  const s = worldToScreen(n.x, n.y);
  const r = (n.type === 'species' ? 14 : n.type === 'population' ? 11 : 8);
  const sr = r;

  // Glow for hovered
  if (hoveredNode === n) {
    ctx.beginPath();
    ctx.arc(s.x, s.y, sr + 6, 0, Math.PI * 2);
    ctx.fillStyle = n.color + '33';
    ctx.fill();
  }

  ctx.beginPath();
  ctx.arc(s.x, s.y, sr, 0, Math.PI * 2);
  ctx.fillStyle = n.color;
  ctx.fill();
  ctx.strokeStyle = 'rgba(255,255,255,0.2)';
  ctx.lineWidth = 1;
  ctx.stroke();

  // Label
  ctx.fillStyle = '#ccc';
  const fontSize = Math.max(10, 11);
  ctx.font = `${fontSize}px 'Segoe UI', system-ui, sans-serif`;
  ctx.textAlign = 'center';
  ctx.fillText(n.label, s.x, s.y + sr + fontSize + 2);
}

function draw() {
  ctx.clearRect(0, 0, W, H);
  // Edges first
  for (const e of DATA.edges) drawEdge(e);
  // Nodes on top
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

canvas.addEventListener('mousedown', (e) => {
  dragging = true;
  dragStartX = e.clientX; dragStartY = e.clientY;
  camStartX = camX; camStartY = camY;
  canvas.style.cursor = 'grabbing';
});
canvas.addEventListener('mousemove', (e) => {
  if (dragging) {
    camX = camStartX - (e.clientX - dragStartX) / zoom;
    camY = camStartY - (e.clientY - dragStartY) / zoom;
    draw();
  }
  // Hover detection
  const w = screenToWorld(e.offsetX, e.offsetY);
  let found = null;
  for (const n of DATA.nodes) {
    const dx = n.x - w.x, dy = n.y - w.y;
    if (dx * dx + dy * dy < 400) { found = n; break; }
  }
  if (found !== hoveredNode) {
    hoveredNode = found;
    draw();
    if (found) {
      tooltip.style.display = 'block';
      tooltip.style.left = (e.clientX + 14) + 'px';
      tooltip.style.top = (e.clientY - 10) + 'px';
      tooltip.innerHTML = `<strong>${found.label}</strong><br>Type: ${found.type}` +
        (found.is_hybrid ? '<br>Hybrid: yes' : '');
    } else {
      tooltip.style.display = 'none';
    }
  } else if (found) {
    tooltip.style.left = (e.clientX + 14) + 'px';
    tooltip.style.top = (e.clientY - 10) + 'px';
  }
});
canvas.addEventListener('mouseup', () => {
  dragging = false;
  canvas.style.cursor = 'default';
});

// Init
centreView();
resize();
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
