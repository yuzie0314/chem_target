"""Molecule visualization — SVG output with color-coded functional group highlighting
and a human-readable legend.
"""

import os
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from constants.fg_smarts import FG_SMARTS

# ── Color palette (r, g, b) in 0-1 range, distinct and print-friendly ─────────
_PALETTE: list[tuple[float, float, float]] = [
    (0.92, 0.27, 0.27),   # red
    (0.20, 0.60, 0.86),   # blue
    (0.24, 0.70, 0.44),   # green
    (0.96, 0.60, 0.13),   # orange
    (0.66, 0.40, 0.80),   # purple
    (0.13, 0.70, 0.67),   # teal
    (0.93, 0.40, 0.65),   # pink
    (0.55, 0.76, 0.24),   # lime
    (0.87, 0.72, 0.20),   # gold
    (0.40, 0.56, 0.80),   # steel blue
    (0.80, 0.47, 0.20),   # brown
    (0.47, 0.80, 0.73),   # cyan
    (0.78, 0.78, 0.30),   # olive
    (0.65, 0.30, 0.47),   # burgundy
    (0.30, 0.47, 0.35),   # dark green
]

# Legend layout constants
_MOL_HEIGHT   = 400   # molecule drawing area height (px)
_MOL_WIDTH    = 600
_ROW_HEIGHT   = 24    # height of each legend row (px)
_SWATCH_W     = 16    # color swatch width
_SWATCH_H     = 14    # color swatch height
_LEGEND_PAD   = 18    # padding above legend
_LEGEND_BOT   = 16    # extra padding below last legend row


def _rgb_to_hex(r: float, g: float, b: float) -> str:
    """Convert 0-1 float RGB to CSS hex string."""
    return "#{:02x}{:02x}{:02x}".format(
        int(r * 255), int(g * 255), int(b * 255)
    )


def _build_legend_svg(
    fg_labels: list[str],
    colors: list[tuple[float, float, float]],
    y_offset: int,
    width: int,
) -> str:
    """Return SVG markup for a color-coded legend block.

    Args:
        fg_labels: Ordered list of functional group names to show.
        colors:    Matching list of (r, g, b) tuples.
        y_offset:  Top y-coordinate to start drawing from.
        width:     Total SVG width (for centering).
    """
    lines = []
    left = 16   # left margin

    for i, (label, color) in enumerate(zip(fg_labels, colors)):
        y = y_offset + i * _ROW_HEIGHT
        hex_color = _rgb_to_hex(*color)

        # Color swatch
        lines.append(
            f'<rect x="{left}" y="{y}" width="{_SWATCH_W}" height="{_SWATCH_H}" '
            f'fill="{hex_color}" rx="2" />'
        )
        # Label text
        lines.append(
            f'<text x="{left + _SWATCH_W + 6}" y="{y + _SWATCH_H - 1}" '
            f'font-size="12" font-family="sans-serif" fill="#222">'
            f'{label}</text>'
        )

    return "\n".join(lines)


def draw_compounds(
    compounds: dict[str, str],
    output_dir: str = "output/images",
) -> None:
    """Draw each compound as an SVG with color-coded FG highlighting and a legend.

    Each functional group that matches the molecule gets a unique color.
    Atoms belonging to multiple FGs keep the color of the last matched FG
    (overlap is rare for the SMARTS in FG_SMARTS).

    Args:
        compounds:  {compound_name: smiles}
        output_dir: Directory to write SVG files into. Created if absent.
    """
    os.makedirs(output_dir, exist_ok=True)

    for comp_name, smiles in compounds.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"  ! Skipped {comp_name}: invalid SMILES")
            continue

        # ── Match each FG and assign colors ───────────────────────────────
        atom_colors: dict[int, tuple[float, float, float]] = {}
        bond_colors: dict[int, tuple[float, float, float]] = {}
        matched_fgs: list[tuple[str, tuple[float, float, float]]] = []

        color_idx = 0
        for fg_label, smarts in FG_SMARTS.items():
            pattern = Chem.MolFromSmarts(smarts)
            if pattern is None:
                continue
            matches = mol.GetSubstructMatches(pattern)
            if not matches:
                continue

            color = _PALETTE[color_idx % len(_PALETTE)]
            newly_colored = 0

            for match in matches:
                for atom_idx in match:
                    if atom_idx not in atom_colors:
                        atom_colors[atom_idx] = color
                        newly_colored += 1
                for i in range(len(match)):
                    for j in range(i + 1, len(match)):
                        bond = mol.GetBondBetweenAtoms(match[i], match[j])
                        if bond and bond.GetIdx() not in bond_colors:
                            bond_colors[bond.GetIdx()] = color

            # Only add to legend if this FG actually colored new atoms.
            # Prevents duplicate-SMARTS entries (e.g. Aromatic ring == Benzene)
            # and sub-pattern entries (e.g. [OH] inside C(=O)[OH]) from appearing
            # in the legend with a color that doesn't show on the molecule.
            if newly_colored > 0:
                color_idx += 1
                matched_fgs.append((fg_label, color))

        # ── Draw molecule ─────────────────────────────────────────────────
        drawer = rdMolDraw2D.MolDraw2DSVG(_MOL_WIDTH, _MOL_HEIGHT)
        drawer.drawOptions().addAtomIndices = False
        drawer.DrawMolecule(
            mol,
            highlightAtoms=list(atom_colors.keys()),
            highlightAtomColors=atom_colors,
            highlightBonds=list(bond_colors.keys()),
            highlightBondColors=bond_colors,
        )
        drawer.FinishDrawing()
        mol_svg = drawer.GetDrawingText()

        # ── Append legend into the SVG ────────────────────────────────────
        if matched_fgs:
            n_rows  = len(matched_fgs)
            legend_h = _LEGEND_PAD + n_rows * _ROW_HEIGHT + _LEGEND_BOT
            total_h  = _MOL_HEIGHT + legend_h

            fg_labels = [label for label, _ in matched_fgs]
            fg_colors = [color for _, color in matched_fgs]

            legend_svg = _build_legend_svg(fg_labels, fg_colors, _MOL_HEIGHT + _LEGEND_PAD, _MOL_WIDTH)

            # Update height AND viewBox so the legend isn't clipped
            import re as _re
            mol_svg = _re.sub(
                r"(width='600px'\s+height=')(\d+)(px')",
                f"\\g<1>{total_h}\\g<3>",
                mol_svg,
            )
            mol_svg = _re.sub(
                r"(width=\"600px\"\s+height=\")(\d+)(px\")",
                f"\\g<1>{total_h}\\g<3>",
                mol_svg,
            )
            mol_svg = _re.sub(
                r"viewBox='0 0 \d+ \d+'",
                f"viewBox='0 0 {_MOL_WIDTH} {total_h}'",
                mol_svg,
            )
            mol_svg = _re.sub(
                r'viewBox="0 0 \d+ \d+"',
                f'viewBox="0 0 {_MOL_WIDTH} {total_h}"',
                mol_svg,
            )
            mol_svg = mol_svg.replace("</svg>", f"{legend_svg}\n</svg>")

        svg_path = os.path.join(output_dir, f"{comp_name}_fg.svg")
        with open(svg_path, "w", encoding="utf-8") as f:
            f.write(mol_svg)

        print(f"  Saved: {svg_path}")
