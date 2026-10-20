# -*- coding: utf-8 -*-
"""
tcs_figstyle.py — unified figure style for all TCS manuscript figures.

Design goals:
- One font family everywhere (Arial, falling back to Liberation Sans / DejaVu Sans).
- No baked-in figure numbers or suptitles: numbering lives in the captions.
- Consistent axes: outward ticks, 1 pt spines, white background.
- Colorblind-safe Okabe-Ito palette available as OI.
- Exports: SVG (fonts as paths) + 300 dpi PNG.

Usage at the top of a figure script, right after `import matplotlib.pyplot as plt`:
    from tcs_figstyle import apply_style, OI, panel_label
    apply_style()
"""
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt


def apply_style(base_size=12):
    plt.rcParams.update({
        # fonts
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans'],
        'font.size': base_size,
        'mathtext.fontset': 'custom',
        'mathtext.rm': 'Arial',
        'mathtext.it': 'Arial:italic',
        'mathtext.bf': 'Arial:bold',
        'svg.fonttype': 'path',
        # axes
        'axes.linewidth': 1.0,
        'axes.labelsize': base_size + 2,
        'axes.titlesize': base_size + 2,
        'axes.facecolor': 'white',
        'figure.facecolor': 'white',
        # ticks
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'xtick.labelsize': base_size,
        'ytick.labelsize': base_size,
        'xtick.major.width': 1.0,
        'ytick.major.width': 1.0,
        'xtick.minor.width': 0.8,
        'ytick.minor.width': 0.8,
        # legend
        'legend.fontsize': base_size - 1,
        'legend.framealpha': 0.9,
        'legend.edgecolor': '0.8',
        # saving
        'savefig.dpi': 300,
        'savefig.facecolor': 'white',
    })


# Okabe-Ito colorblind-safe palette
OI = {
    'blue': '#0072B2',
    'orange': '#E69F00',
    'green': '#009E73',
    'vermillion': '#D55E00',
    'sky': '#56B4E9',
    'purple': '#CC79A7',
    'yellow': '#F0E442',
    'black': '#000000',
}


def panel_label(ax, label, x=-0.14, y=1.03, fontsize=16):
    """Nature-style bold lowercase panel letter at the top-left corner."""
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=fontsize, fontweight='bold', va='bottom', ha='left')


def save_fig(fig, stem, dpi=300):
    """Save both SVG (fonts as paths) and PNG."""
    fig.savefig(stem + '.svg', bbox_inches='tight', facecolor='white')
    fig.savefig(stem + '.png', dpi=dpi, bbox_inches='tight', facecolor='white')
