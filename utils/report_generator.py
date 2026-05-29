"""HTML report generator for target prediction results.

Produces self-contained HTML files — no external CSS/JS dependencies.
All styles are inline; images are embedded as inline SVG.

Each report contains:
  1. Compound header (name, SMILES, date)
  2. Molecule structure with FG highlights (inline SVG)
  3. Detected functional groups (colour-coded chips)
  4. Predicted target classes (ranked table + score bars)
  5. Key binding residues (horizontal score bars)
  6. Scoring methodology (plain-language explanation for clients)
  7. Data sources & version footer

Batch usage also generates an index.html summary page.
"""

from __future__ import annotations

import html as _html_module
from datetime import date
from math import log
from pathlib import Path
from typing import Any

import pandas as pd
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

_VERSION = "1.0"
_BIOLIP_EVENTS = "46,114"
_N_FGS_TOTAL = 32   # 31 SMARTS + Steroid Python

# ── Color palette for FG chips / highlights ────────────────────────────────────
_PALETTE_HEX: list[str] = [
    "#e74c3c", "#3498db", "#27ae60", "#f39c12", "#9b59b6",
    "#1abc9c", "#e91e63", "#8bc34a", "#ff9800", "#5c6bc0",
    "#795548", "#00bcd4", "#cddc39", "#8d6e63", "#f06292",
]


def _hex_to_rgb01(h: str) -> tuple[float, float, float]:
    """Convert #rrggbb to (r, g, b) in [0, 1]."""
    h = h.lstrip("#")
    return (int(h[0:2], 16) / 255, int(h[2:4], 16) / 255, int(h[4:6], 16) / 255)


# ── Molecule SVG with FG highlights ───────────────────────────────────────────

def _mol_svg(
    smiles: str,
    fgs_detected: list[str],
    width: int = 480,
    height: int = 300,
) -> str:
    """Return an inline SVG string of the molecule with FG atoms highlighted."""
    from constants.fg_smarts import FG_SMARTS

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (
            '<p style="color:#c0392b;padding:1rem;">'
            "⚠ Cannot parse SMILES — structure not available.</p>"
        )

    atom_colors: dict[int, tuple[float, float, float]] = {}
    highlight_atoms: list[int] = []
    atom_radii: dict[int, float] = {}

    for fg_idx, fg_name in enumerate(fgs_detected):
        color = _hex_to_rgb01(_PALETTE_HEX[fg_idx % len(_PALETTE_HEX)])

        if fg_name in FG_SMARTS:
            pat = Chem.MolFromSmarts(FG_SMARTS[fg_name])
            if pat:
                for match in mol.GetSubstructMatches(pat):
                    for atom_idx in match:
                        atom_colors.setdefault(atom_idx, color)
                        highlight_atoms.append(atom_idx)
                        atom_radii[atom_idx] = 0.35
        elif fg_name == "Steroid":
            ri = mol.GetRingInfo()
            for ring in ri.AtomRings():
                for atom_idx in ring:
                    atom_colors.setdefault(atom_idx, color)
                    highlight_atoms.append(atom_idx)
                    atom_radii[atom_idx] = 0.35

    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.drawOptions().addStereoAnnotation = True
    try:
        drawer.DrawMolecule(
            mol,
            highlightAtoms=list(set(highlight_atoms)),
            highlightAtomColors=atom_colors,
            highlightBonds=[],
            highlightAtomRadii=atom_radii,
        )
    except Exception:
        drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return svg.replace("<?xml version='1.0' encoding='utf-8'?>\n", "")


# ── CSS / HTML building blocks ─────────────────────────────────────────────────

_CSS = """
* { box-sizing: border-box; margin: 0; padding: 0; }
body {
  font-family: 'Segoe UI', Arial, sans-serif;
  background: #f0f2f5;
  color: #2c3e50;
  padding: 24px;
}
.report {
  max-width: 860px;
  margin: 0 auto;
  background: #fff;
  border-radius: 10px;
  box-shadow: 0 2px 16px rgba(0,0,0,.10);
  overflow: hidden;
}

/* Header */
.rpt-header {
  background: linear-gradient(135deg, #1e3a5f 0%, #2980b9 100%);
  color: #fff;
  padding: 28px 32px 22px;
}
.rpt-header h1 { font-size: 1.45rem; font-weight: 700; margin-bottom: 6px; }
.rpt-header .meta { font-size: .82rem; opacity: .80; line-height: 1.7; }

/* Sections */
.section { padding: 24px 32px; border-bottom: 1px solid #eaecef; }
.section:last-child { border-bottom: none; }
.section h2 {
  font-size: 1.05rem;
  font-weight: 600;
  color: #1e3a5f;
  margin-bottom: 14px;
  padding-bottom: 6px;
  border-bottom: 2px solid #2980b9;
  display: inline-block;
}

/* Molecule */
.mol-wrap { text-align: center; background: #fafbfc; border-radius: 8px;
            padding: 12px; border: 1px solid #e0e4ea; }

/* FG chips */
.fg-chips { display: flex; flex-wrap: wrap; gap: 8px; margin-top: 4px; }
.chip {
  display: inline-flex; align-items: center; gap: 6px;
  padding: 4px 12px; border-radius: 20px;
  font-size: .82rem; font-weight: 600; color: #fff;
  white-space: nowrap;
}
.chip-dot {
  width: 8px; height: 8px; border-radius: 50%;
  background: rgba(255,255,255,.55);
  flex-shrink: 0;
}

/* Target table */
table { width: 100%; border-collapse: collapse; font-size: .88rem; }
th {
  background: #f0f4f8; color: #1e3a5f;
  padding: 9px 12px; text-align: left; font-weight: 600;
  border-bottom: 2px solid #d0d7e2;
}
td { padding: 8px 12px; border-bottom: 1px solid #eaecef; vertical-align: middle; }
tr:hover td { background: #f7f9fc; }
.rank { color: #888; font-size: .78rem; font-weight: 600; text-align: center; }
.score-val { font-weight: 700; color: #2980b9; }
.votes-val { color: #888; font-size: .82rem; }
.evidence { color: #555; font-size: .80rem; font-style: italic; }

/* Score bar (shared for targets + residues) */
.bar-wrap { display: flex; align-items: center; gap: 10px; }
.bar-bg {
  flex: 1; height: 14px; background: #eaecef;
  border-radius: 7px; overflow: hidden; min-width: 80px;
}
.bar-fill {
  height: 100%; border-radius: 7px;
  background: linear-gradient(90deg, #2980b9, #1abc9c);
}
.bar-num { font-size: .78rem; color: #555; min-width: 38px; text-align: right; }

/* Residue rows */
.res-row { display: flex; align-items: center; gap: 12px;
           padding: 5px 0; border-bottom: 1px solid #f0f2f5; }
.res-name {
  font-family: monospace; font-size: .88rem; font-weight: 700;
  color: #1e3a5f; min-width: 40px;
}
.res-label { font-size: .78rem; color: #888; min-width: 36px; }

/* Methodology box */
.method-box {
  background: #f4f8fd;
  border-left: 4px solid #2980b9;
  border-radius: 0 8px 8px 0;
  padding: 18px 22px;
  font-size: .88rem;
  line-height: 1.7;
}
.method-box h3 {
  color: #1e3a5f; font-size: .95rem; font-weight: 700; margin-bottom: 10px;
}
.method-box h4 {
  color: #2980b9; font-size: .85rem; font-weight: 700;
  margin: 14px 0 5px;
}
.method-box p { margin-bottom: 8px; color: #444; }
.method-box .formula {
  font-family: monospace; background: #e8f0fb;
  padding: 6px 12px; border-radius: 4px;
  display: inline-block; margin: 4px 0; font-size: .84rem;
  color: #1e3a5f;
}
.method-box table { font-size: .82rem; margin-top: 6px; }
.method-box td { padding: 4px 10px; border-bottom: 1px solid #dde4ee; }
.method-box td:first-child { font-weight: 600; color: #1e3a5f; }

/* Footer */
.rpt-footer {
  background: #f4f6f9;
  padding: 16px 32px;
  font-size: .78rem;
  color: #888;
  line-height: 1.7;
}
.rpt-footer a { color: #2980b9; text-decoration: none; }
"""


def _chip_html(fg_name: str, idx: int) -> str:
    """Return a colored chip HTML element for a functional group."""
    color = _PALETTE_HEX[idx % len(_PALETTE_HEX)]
    name  = _html_module.escape(fg_name)
    return (
        f'<span class="chip" style="background:{color};">'
        f'<span class="chip-dot"></span>{name}</span>'
    )


def _score_bar_html(value: float, max_value: float, color: str = "") -> str:
    """Return a CSS score bar as an HTML string."""
    pct = min(100, round(value / max_value * 100, 1)) if max_value > 0 else 0
    style = f'width:{pct}%;' + (f'background:{color};' if color else "")
    return (
        f'<div class="bar-wrap">'
        f'<div class="bar-bg"><div class="bar-fill" style="{style}"></div></div>'
        f'<span class="bar-num">{value:.2f}</span>'
        f'</div>'
    )


def _methodology_html(fgs_detected: list[str], fg_db: dict) -> str:
    """Return the methodology explanation section as HTML."""
    # Compute IDF for example targets to show in the explanation table
    from math import log as _log
    tc_count: dict[str, int] = {}
    for entry in fg_db.values():
        for tc in entry.get("known_target_classes", []):
            tc_count[tc] = tc_count.get(tc, 0) + 1

    n_total = _N_FGS_TOTAL
    # Pick 4 illustrative targets: 2 generic, 2 specific
    illustrative = sorted(tc_count.items(), key=lambda x: x[1])
    # Bottom 2 (most specific) + top 2 (most generic)
    examples = illustrative[:2] + illustrative[-2:]

    rows = ""
    for tc, count in sorted(examples, key=lambda x: x[1], reverse=True):
        idf = round(_log(n_total / count), 2)
        rows += (
            f"<tr><td>{_html_module.escape(tc)}</td>"
            f"<td style='text-align:center'>{count}</td>"
            f"<td style='text-align:center'>{idf:.2f}</td></tr>"
        )

    detected_str = ", ".join(
        f"<strong>{_html_module.escape(fg)}</strong>" for fg in fgs_detected
    ) if fgs_detected else "<em>(none detected)</em>"

    return f"""
<div class="method-box">
  <h3>🔬 How are targets ranked?</h3>
  <p>
    The scoring pipeline has two complementary stages.
    For this compound, the following functional groups were identified:<br>
    {detected_str}
  </p>

  <h4>Stage 1 — Structural database matching (residue scoring)</h4>
  <p>
    Each functional group is matched against <strong>{_BIOLIP_EVENTS} protein–ligand
    binding events</strong> from the BioLiP 2.0 structural database.
    For every amino acid residue type, we count how often it appears in binding
    sites that contain each of your compound's functional groups.
  </p>
  <p>
    <strong>Fairness correction (z-score normalisation):</strong>
    Some functional groups are extremely common — Hydroxyl appears in thousands
    of ligands while Endoperoxide appears in fewer than five.
    Without correction, common groups would always dominate the final score.
    We normalise each group's contribution to a common statistical scale
    (mean&nbsp;= 0, std&nbsp;= 1) so that every functional group contributes
    equally regardless of how often it appears in the database.
  </p>

  <h4>Stage 2 — Target class voting (IDF weighting)</h4>
  <p>
    Each detected functional group "votes" for the target protein classes listed
    in our curated pharmacological annotation database.
    Not all votes carry equal weight:
  </p>
  <p>
    <span class="formula">
      Score = votes &times; log( {n_total} &divide; number&nbsp;of&nbsp;FGs&nbsp;annotated&nbsp;for&nbsp;this&nbsp;target )
    </span>
  </p>
  <p>
    A target associated with <em>many</em> different functional groups (e.g. "kinase")
    is a <em>generic</em> signal — it tells us little about what makes your compound
    unique.  A target associated with <em>only one or two</em> functional groups
    (e.g. "VKORC1" for Coumarin, or "antimalarial target" for Endoperoxide) is a
    <em>specific, high-confidence</em> signal and receives a higher score.
  </p>
  <table>
    <tr>
      <th>Target class</th>
      <th style="text-align:center">FGs annotated</th>
      <th style="text-align:center">Weight (IDF)</th>
    </tr>
    {rows}
  </table>
  <p style="margin-top:10px;color:#666;font-size:.80rem;">
    Higher weight = more specific prediction. Lower weight = generic label
    (still relevant, but less discriminating).
  </p>
</div>
"""


# ── Main HTML report function ──────────────────────────────────────────────────

def generate_html_report(
    name: str,
    smiles: str,
    pred: dict,
    fg_db: dict,
    top_residues: int = 10,
) -> str:
    """Generate a self-contained HTML report for one compound prediction.

    Args:
        name:         Compound display name.
        smiles:       Input SMILES string.
        pred:         Result dict from target_predictor.predict().
        fg_db:        functional_groups dict from fg_database.json.
        top_residues: How many residues to show in the bar chart.

    Returns:
        HTML string (self-contained, embeds all CSS and SVG inline).
    """
    today     = date.today().strftime("%Y-%m-%d")
    fgs       = pred.get("fgs_detected", [])
    tc_df: pd.DataFrame = pred.get("target_class_votes", pd.DataFrame())
    rs_df: pd.DataFrame = pred.get("residue_scores",     pd.DataFrame())
    warning   = pred.get("warning") or ""

    esc_name  = _html_module.escape(name)
    esc_smiles = _html_module.escape(smiles)

    # ── Structure SVG ──────────────────────────────────────────────────────────
    mol_svg_html = _mol_svg(smiles, fgs)

    # ── FG chips ───────────────────────────────────────────────────────────────
    if fgs:
        chips_html = "".join(_chip_html(fg, i) for i, fg in enumerate(fgs))
    else:
        chips_html = '<span style="color:#c0392b;">⚠ No functional groups detected</span>'

    # ── Warning block ──────────────────────────────────────────────────────────
    warn_block = ""
    if warning:
        warn_block = (
            f'<div style="margin:0 32px 16px;padding:12px 16px;'
            f'background:#fff3cd;border-left:4px solid #f39c12;'
            f'border-radius:0 6px 6px 0;font-size:.86rem;color:#856404;">'
            f'⚠ {_html_module.escape(warning)}</div>'
        )

    # ── Target class table ────────────────────────────────────────────────────
    if not tc_df.empty:
        max_score = tc_df["score"].max()
        tc_rows   = ""
        for rank, (_, row) in enumerate(tc_df.iterrows(), 1):
            bar = _score_bar_html(row["score"], max_score)
            tc_rows += (
                f"<tr>"
                f"<td class='rank'>#{rank}</td>"
                f"<td><strong>{_html_module.escape(str(row['target_class']))}</strong></td>"
                f"<td>{bar}</td>"
                f"<td class='votes-val'>{int(row['votes'])}</td>"
                f"<td class='evidence'>{_html_module.escape(str(row['evidence_fgs']))}</td>"
                f"</tr>"
            )
        tc_html = f"""
<table>
  <tr>
    <th style="width:36px">#</th>
    <th>Target class</th>
    <th>Score (IDF-weighted)</th>
    <th style="width:56px">Votes</th>
    <th>Evidence FGs</th>
  </tr>
  {tc_rows}
</table>"""
    else:
        tc_html = '<p style="color:#888;">No target class annotations found for detected FGs.</p>'

    # ── Residue bar chart ─────────────────────────────────────────────────────
    if not rs_df.empty:
        top_rs    = rs_df.head(top_residues)
        max_rscore = top_rs["score"].max()
        res_rows  = ""
        for aa, row in top_rs.iterrows():
            bar = _score_bar_html(row["score"], max_rscore)
            res_rows += (
                f'<div class="res-row">'
                f'<span class="res-name">{_html_module.escape(str(aa))}</span>'
                f'<span class="res-label">{_html_module.escape(str(row["residue_name"]))}</span>'
                f'{bar}'
                f'</div>'
            )
        rs_html = res_rows
    else:
        rs_html = '<p style="color:#888;">No residue scores available.</p>'

    # ── Methodology ────────────────────────────────────────────────────────────
    method_html = _methodology_html(fgs, fg_db)

    # ── Full HTML ─────────────────────────────────────────────────────────────
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Target Prediction — {esc_name}</title>
  <style>{_CSS}</style>
</head>
<body>
<div class="report">

  <!-- Header -->
  <div class="rpt-header">
    <h1>🎯 Target Prediction Report</h1>
    <div class="meta">
      <strong>Compound:</strong> {esc_name} &nbsp;|&nbsp;
      <strong>Date:</strong> {today} &nbsp;|&nbsp;
      <strong>Tool:</strong> chem_target v{_VERSION}<br>
      <strong>SMILES:</strong> <span style="font-family:monospace;font-size:.80rem;">{esc_smiles}</span>
    </div>
  </div>

  {warn_block}

  <!-- Structure -->
  <div class="section">
    <h2>Molecular Structure</h2>
    <div class="mol-wrap">{mol_svg_html}</div>
  </div>

  <!-- Functional groups -->
  <div class="section">
    <h2>Detected Functional Groups ({len(fgs)})</h2>
    <div class="fg-chips">{chips_html}</div>
  </div>

  <!-- Target classes -->
  <div class="section">
    <h2>Predicted Target Classes</h2>
    {tc_html}
  </div>

  <!-- Residues -->
  <div class="section">
    <h2>Key Binding Residues</h2>
    <p style="font-size:.82rem;color:#888;margin-bottom:12px;">
      Z-score-normalised co-occurrence with detected FGs in BioLiP 2.0.
      Higher score = stronger structural evidence for binding site involvement.
    </p>
    {rs_html}
  </div>

  <!-- Methodology -->
  <div class="section">
    <h2>Scoring Methodology</h2>
    {method_html}
  </div>

  <!-- Footer -->
  <div class="rpt-footer">
    <strong>Data sources:</strong>
    BioLiP 2.0 ({_BIOLIP_EVENTS} binding events) ·
    RCSB CCD (SMILES for 6,002 ligands) ·
    PubChem / ChEBI (FG metadata) ·
    RCSB PDB (3D poses)<br>
    <strong>chem_target v{_VERSION}</strong> ·
    Generated {today} ·
    For research use only — not a substitute for experimental validation.
  </div>

</div>
</body>
</html>"""


# ── Batch index page ───────────────────────────────────────────────────────────

def generate_index_html(
    results: list[dict[str, Any]],
    title: str = "Target Prediction Batch Report",
) -> str:
    """Generate a summary index HTML linking to individual compound reports.

    Args:
        results: List of dicts with keys:
                   name, smiles, html_filename, n_fgs, fgs,
                   top_target, top_score, warning
        title:   Page title.

    Returns:
        HTML string for the index page.
    """
    today = date.today().strftime("%Y-%m-%d")
    n_ok  = sum(1 for r in results if not r.get("warning"))
    n_err = len(results) - n_ok

    rows = ""
    for i, r in enumerate(results, 1):
        link     = _html_module.escape(r.get("html_filename", "#"))
        name     = _html_module.escape(r.get("name", ""))
        target   = _html_module.escape(str(r.get("top_target", "—")))
        score    = r.get("top_score", "")
        score_s  = f"{score:.2f}" if isinstance(score, float) else "—"
        fgs      = _html_module.escape(str(r.get("fgs", "—")))
        n_fgs    = r.get("n_fgs", 0)
        warn     = r.get("warning", "")
        warn_td  = (
            f'<td style="color:#c0392b;font-size:.80rem;">'
            f'{_html_module.escape(str(warn)[:60])}</td>'
            if warn else
            '<td style="color:#27ae60;">✓</td>'
        )

        rows += (
            f"<tr>"
            f"<td class='rank'>{i}</td>"
            f"<td><a href='{link}' style='color:#2980b9;font-weight:600;'>{name}</a></td>"
            f"<td>{n_fgs}</td>"
            f"<td style='font-size:.80rem;color:#555;'>{fgs}</td>"
            f"<td><strong>{target}</strong></td>"
            f"<td class='score-val'>{score_s}</td>"
            f"{warn_td}"
            f"</tr>"
        )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>{_html_module.escape(title)}</title>
  <style>
    {_CSS}
    .index-wrap {{ max-width: 1000px; margin: 0 auto; background:#fff;
                  border-radius:10px; box-shadow:0 2px 16px rgba(0,0,0,.10);
                  overflow:hidden; }}
    .stat-bar {{ display:flex; gap:24px; padding:16px 32px;
                background:#f4f8fd; border-bottom:1px solid #dde4ee; }}
    .stat {{ font-size:.85rem; }}
    .stat strong {{ font-size:1.2rem; color:#1e3a5f; display:block; }}
  </style>
</head>
<body>
<div class="index-wrap">

  <div class="rpt-header">
    <h1>📋 Batch Prediction Report</h1>
    <div class="meta">
      <strong>Date:</strong> {today} &nbsp;|&nbsp;
      <strong>Tool:</strong> chem_target v{_VERSION} &nbsp;|&nbsp;
      <strong>Total compounds:</strong> {len(results)}
    </div>
  </div>

  <div class="stat-bar">
    <div class="stat"><strong>{len(results)}</strong>Total</div>
    <div class="stat"><strong style="color:#27ae60">{n_ok}</strong>Predicted</div>
    <div class="stat"><strong style="color:#c0392b">{n_err}</strong>Failed</div>
  </div>

  <div class="section">
    <h2>Compounds</h2>
    <table>
      <tr>
        <th style="width:36px">#</th>
        <th>Compound</th>
        <th style="width:50px">FGs</th>
        <th>Detected functional groups</th>
        <th>Top predicted target</th>
        <th style="width:72px">Score</th>
        <th style="width:80px">Status</th>
      </tr>
      {rows}
    </table>
  </div>

  <div class="rpt-footer">
    <strong>chem_target v{_VERSION}</strong> ·
    Generated {today} ·
    For research use only — not a substitute for experimental validation.
  </div>

</div>
</body>
</html>"""
