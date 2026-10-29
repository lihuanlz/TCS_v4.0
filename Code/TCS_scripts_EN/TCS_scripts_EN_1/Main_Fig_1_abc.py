#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fig 1 (panels a+b+c) — SINGLE-FILE GENERATOR, Nature Methods double-column spec.

Run:
    python3 Main_Fig_1_abc.py

Output (saved next to this script, wherever it lives):
    Fig1_总图_abc.png / Fig1_总图_abc.pdf   (183.0 x ~189.3 mm, 600 dpi)

Self-contained: panels a, b, c are drawn in memory (no intermediate files,
no /tmp dependency) and placed at exact millimetre positions.

IMPORTANT — keep the ./fonts folder next to this script.  It carries the
Liberation Sans/Serif font files (Arial/Times-metric, OFL-licensed); the
script registers them at runtime, so the figure renders IDENTICALLY on any
machine.  Without it, matplotlib falls back to wider fonts (DejaVu) and
fixed-position texts may overflow.  Layout:
  page 183 mm wide (max 183 x 247 mm); all in-panel text stays within the
  5-7 pt band at final placement:
    a  placed 182.6 mm wide (native; min font 6.5 pt)
    b  placed 117.0 mm wide (6 pt -> 5.0 pt, at the floor)
    c  placed  61.5 mm wide (portrait re-render, 14 pt base -> 7.0 pt)
  b and c are top-aligned on one baseline; panel labels 8 pt bold,
  sans-serif (Liberation Sans = Arial metric); each letter carries a short
  bold title whose wording is taken from the manuscript's Fig. 1 caption.
Panel d is NOT here — it ships separately as an Extended-Data sheet via
Fig1_总图_拼版_脚本.py (four panels at the 5-pt floor would need ~285 mm,
exceeding one page).
"""
import io
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from matplotlib import font_manager as _fm

# Register the font files shipped in ./fonts next to this script, so the
# output is identical on every machine (no dependence on the OS font scan).
_FONT_DIRS = [os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fonts'),
              os.path.join(os.getcwd(), 'fonts')]
_REGISTERED = 0
for _d in _FONT_DIRS:
    if os.path.isdir(_d):
        for _f in sorted(os.listdir(_d)):
            if _f.lower().endswith(('.ttf', '.otf')):
                try:
                    _fm.fontManager.addfont(os.path.join(_d, _f))
                    _REGISTERED += 1
                except Exception:
                    pass
        if _REGISTERED:
            print(f'[fonts] registered {_REGISTERED} font files from {_d}')
        break
if not _REGISTERED:
    print('[fonts] WARNING: no ./fonts folder found next to the script or in the'
          ' current directory; falling back to system fonts (text may overflow)')

def _fig_to_image(fig, dpi=600):
    """Render a finished figure to an in-memory PIL image (tight, 600 dpi)."""
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=dpi, bbox_inches='tight',
                pad_inches=0.02)
    plt.close(fig)
    buf.seek(0)
    return Image.open(buf).copy()


# =============================================================== panel a
def render_a():
    """
    Fig 1a (v3.1): "Epitope-anchored quantification in FOUR steps",
    aligned with the four-layer stochastic model of SI S2b.1:
      (1) Molecular sampling: W ~ Poisson(M)
      (2) Finite-site binding: exact binding identity
          <B> = <(M-B)(Omega-B)>/KV (equation 1)
      (3) Partition allocation: Omega sites over N partitions of
          beta = Omega/N sites; hypergeometric allocation of C complexes
      (4) Observation: analog intensity vs digital counting
          Z ~ Binomial(N, Ppos); absolute count M = kappa*Omega*xi (eq. 4)
    Output: Fig1a svg (fonts as paths) + 600dpi png
    Canvas: 183 mm (7.2 in) double-column width, Liberation Sans.
    """
    import os
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import (Circle, Ellipse, FancyBboxPatch,
                                    FancyArrowPatch, Polygon, Rectangle)
    import matplotlib.gridspec as gridspec

    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Liberation Sans', 'Arial', 'DejaVu Sans']
    plt.rcParams['svg.fonttype'] = 'path'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.it'] = 'Liberation Sans:italic'
    plt.rcParams['mathtext.rm'] = 'Liberation Sans'

    C_TGT  = '#009E73'
    C_EPI  = '#00573F'
    C_AB   = '#0072B2'
    NAVY   = '#1F2D4E'
    CARD_FC, CARD_EC = '#F5F7FA', '#C9D2DC'
    BOX_FC, BOX_EC   = '#FFF4D6', '#D9BC5F'
    LIQ_EC = '#C0574F'
    GREY   = '#5A6577'
    INK    = '#333333'

    # ------------------------------------------------------------ helpers
    def draw_site(ax, x, y, s, color=C_AB, lw=1.8, z=3):
        ax.plot([x, x], [y - 0.60*s, y], color=color, lw=lw,
                solid_capstyle='round', zorder=z)
        w, h = 0.53*s, 0.68*s
        nw, nd = 0.19*s, 0.29*s
        pad = [(x - w, y), (x + w, y), (x + w, y + h), (x + nw, y + h),
               (x, y + h - nd), (x - nw, y + h), (x - w, y + h)]
        ax.add_patch(Polygon(pad, closed=True, fc=color, ec='#013A5C',
                             lw=max(0.6, 0.55*lw), zorder=z))
        ax.add_patch(Polygon([(x - nw, y + h), (x, y + h - nd), (x + nw, y + h)],
                             closed=True, fc='#01496F', ec='none', zorder=z + 1))

    def _blob(x, y, r, seed, rough=0.12, n=9):
        rng = np.random.default_rng(seed)
        ang = np.linspace(0, 2*np.pi, n, endpoint=False)
        rad = r * (1 + rough * (rng.random(n) * 2 - 1))
        pts = np.column_stack([rad*np.cos(ang), rad*np.sin(ang)])
        for _ in range(3):
            out = []
            for i in range(len(pts)):
                p, q = pts[i], pts[(i + 1) % len(pts)]
                out.append(0.75*p + 0.25*q)
                out.append(0.25*p + 0.75*q)
            pts = np.array(out)
        return pts + np.array([x, y])

    def mol(ax, x, y, r, motif_deg=0, seed=0, z=4):
        pts = _blob(x, y, r, seed)
        ax.add_patch(Polygon(pts, closed=True, fc='#00785A', ec='none', zorder=z))
        ov = (pts - [x, y]) * 0.87 + [x - 0.05*r, y + 0.07*r]
        ax.add_patch(Polygon(ov, closed=True, fc=C_TGT, ec='none', zorder=z + 1))
        hi = (pts - [x, y]) * 0.55 + [x - 0.16*r, y + 0.18*r]
        ax.add_patch(Polygon(hi, closed=True, fc='#46C49C', ec='none',
                             alpha=0.9, zorder=z + 2))
        ax.add_patch(Polygon(pts, closed=True, fc='none', ec=C_EPI, lw=1.3,
                             zorder=z + 3))
        a = np.radians(motif_deg)
        base = np.array([[-0.11, r - 0.02], [0.11, r - 0.02], [0.0, r + 0.22]])
        R = np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
        tri = base @ R.T + np.array([x, y])
        ax.add_patch(Polygon(tri, closed=True, fc=C_EPI, ec=C_EPI, zorder=z + 4))

    def badge(ax, x, y, n):
        ax.add_patch(Circle((x, y), 0.52, fc=NAVY, ec='none', zorder=6))
        ax.text(x, y - 0.02, str(n), color='white', fontsize=9.5, weight='bold',
                ha='center', va='center', zorder=7)

    def card(ax, w):
        ax.add_patch(FancyBboxPatch((0.22, 0.18), w - 0.24, 13.0,
                     boxstyle='round,pad=0.02,rounding_size=0.28',
                     fc='#D9DFE7', ec='none', zorder=-1))
        ax.add_patch(FancyBboxPatch((0.12, 0.30), w - 0.24, 13.0,
                     boxstyle='round,pad=0.02,rounding_size=0.28',
                     fc=CARD_FC, ec=CARD_EC, lw=1.2, zorder=0))

    # ------------------------------------------------------------ figure
    W = [7.4, 9.6, 9.0, 9.6]
    fig = plt.figure(figsize=(7.2, 3.0))
    gs = gridspec.GridSpec(1, 4, width_ratios=W,
                           left=0.004, right=0.996, top=0.985, bottom=0.015,
                           wspace=0.28)
    axs = []
    for i in range(4):
        ax = fig.add_subplot(gs[i])
        ax.set_xlim(0, W[i]); ax.set_ylim(0, 13.4)
        ax.set_aspect('equal'); ax.axis('off')
        card(ax, W[i])
        axs.append(ax)
    ax1, ax2, ax3, ax4 = axs

    FS_T, FS_C = 8.3, 6.1

    # ================================================ STEP 1 - sampling
    badge(ax1, 0.72, 12.45, 1)
    ax1.text(1.32, 12.45, 'Molecular\nsampling', fontsize=FS_T, weight='bold',
             color=NAVY, ha='left', va='center', linespacing=1.05)
    cx = 3.9
    Y_CONE, Y_TIP = 6.95, 4.60
    HW = 1.35
    def hw(y):
        if y >= Y_CONE:
            return HW
        return 0.16 + (y - Y_TIP) * (HW - 0.16) / (Y_CONE - Y_TIP)
    _ta = np.linspace(np.pi, 2 * np.pi, 25)
    tube = ([(cx - HW, 10.6), (cx - HW, Y_CONE)]
            + list(zip(cx + 0.16 * np.cos(_ta), Y_TIP + 0.16 * np.sin(_ta)))
            + [(cx + HW, Y_CONE), (cx + HW, 10.6)])
    ax1.add_patch(Polygon(tube, closed=False, fc='none', ec=NAVY, lw=1.5,
                          joinstyle='round', capstyle='round', zorder=2))
    yl = 9.50
    liq = ([(cx - hw(yl) + 0.08, yl), (cx + hw(yl) - 0.08, yl),
            (cx + HW - 0.08, Y_CONE)]
           + list(zip(cx + 0.10 * np.cos(_ta), Y_TIP + 0.10 + 0.10 * np.sin(_ta)))
           + [(cx - HW + 0.08, Y_CONE)])
    liq_clip = Polygon(liq, closed=True, fc='none', ec='none', zorder=1)
    ax1.add_patch(liq_clip)
    bands = ['#F8D6D1', '#F5CAC4', '#F2BDB5', '#EFB0A6', '#ECA499']
    ybs = np.linspace(Y_TIP + 0.02, yl, 6)
    for i in range(5):
        rc = Rectangle((cx - 1.6, ybs[i]), 3.2, ybs[i + 1] - ybs[i],
                       fc=bands[4 - i], ec='none', zorder=1)
        ax1.add_patch(rc)
        rc.set_clip_path(liq_clip)
    ax1.add_patch(Ellipse((cx, 10.6), 2.70, 0.38, fc='#EDF1F5', ec=NAVY, lw=1.5,
                          zorder=2.6))
    ax1.add_patch(Ellipse((cx, yl), 2 * (hw(yl) - 0.10), 0.20, fc='#F8D6D1',
                          ec=LIQ_EC, lw=1.1, zorder=2.4))
    rng = np.random.default_rng(9)
    pts = []
    for _ in range(600):
        if len(pts) >= 8:
            break
        y = 5.9 + 3.1 * rng.random() ** 0.7
        x = cx + (rng.random() * 2 - 1) * (hw(y) - 0.24)
        if all((x - px)**2 + (y - py)**2 > 0.52**2 for px, py in pts):
            pts.append((x, y))
    for (x, y) in pts:
        mol(ax1, x, y, 0.22, motif_deg=float(rng.integers(0, 360)),
            seed=int(rng.integers(1, 10**6)))
    mol(ax1, 1.10, 9.55, 0.22, motif_deg=25, seed=11)
    ax1.text(1.45, 8.75, 'target\nmolecule', fontsize=5.8, color=C_EPI,
             ha='center', va='top', linespacing=1.15)
    ax1.text(cx, 3.55, 'M targets in volume V', fontsize=FS_C, ha='center',
             color=INK)
    ax1.text(cx, 1.55, r'$W \sim \mathrm{Poisson}(M)$', fontsize=10,
             ha='center', va='center', color=NAVY)

    # ========================================== STEP 2 - finite-site binding
    badge(ax2, 0.72, 12.45, 2)
    ax2.text(1.32, 12.45, 'Finite-site binding', fontsize=FS_T, weight='bold',
             color=NAVY, ha='left', va='center')
    ax2.add_patch(FancyBboxPatch((2.10, 5.55), 1.45, 0.28,
                  boxstyle='round,pad=0.02,rounding_size=0.13',
                  fc='#8A97A5', ec=NAVY, lw=0.9, zorder=2))
    ax2.add_patch(Circle((2.83, 5.87), 0.15, fc=NAVY, ec='none', zorder=2))
    draw_site(ax2, 2.83, 6.80, 0.80, lw=2.0)
    mol(ax2, 2.83, 7.85, 0.48, motif_deg=180, seed=4)
    mol(ax2, 2.30, 9.45, 0.33, motif_deg=195, seed=6)
    ax2.add_patch(FancyArrowPatch((2.60, 9.10), (2.60, 8.50),
                  connectionstyle='arc3,rad=-0.30', arrowstyle='-|>',
                  mutation_scale=8, color=GREY, lw=1.1, ls=(0, (3, 2)),
                  alpha=0.85, zorder=3))
    ax2.text(2.10, 10.45, 'incoming target', fontsize=5.8, color=C_EPI,
             ha='center', va='center')
    ax2.add_patch(FancyArrowPatch((1.10, 8.20), (1.10, 7.40), arrowstyle='-|>',
                                  mutation_scale=8, color=NAVY, lw=1.3))
    ax2.add_patch(FancyArrowPatch((1.50, 7.40), (1.50, 8.20), arrowstyle='-|>',
                                  mutation_scale=8, color=NAVY, lw=1.3))
    ax2.text(0.92, 7.80, r'$k_{\mathrm{on}}$', fontsize=6, ha='right',
             va='center', color=NAVY)
    ax2.text(1.68, 7.80, r'$k_{\mathrm{off}}$', fontsize=6, ha='left',
             va='center', color=NAVY)
    tcx, tcy, tr = 7.60, 8.65, 0.92
    tri = [(tcx, tcy + tr),
           (tcx - tr*np.cos(np.pi/6), tcy - tr*np.sin(np.pi/6)),
           (tcx + tr*np.cos(np.pi/6), tcy - tr*np.sin(np.pi/6))]
    ax2.add_patch(Polygon(tri, closed=True, fc='none', ec=NAVY, lw=1.3, zorder=2))
    for (nx, ny), lab in zip(tri, ['M', r'$\Omega$', 'K']):
        ax2.add_patch(Circle((nx, ny), 0.40, fc='white', ec=NAVY, lw=1.3,
                             zorder=4))
        ax2.text(nx, ny, lab, fontsize=7.5, weight='bold', color=NAVY,
                 ha='center', va='center', zorder=5)
    ax2.text(tcx, tcy, r'$\times s$', fontsize=9, weight='bold',
             color='#B03030', ha='center', va='center', zorder=5)
    ax2.text(tcx, 6.95, 'uniform scaling', fontsize=5.8, color=GREY,
             ha='center', va='center')
    ax2.text(4.80, 4.70, r'$B$: bound complexes;  $p \equiv B/\Omega$: occupancy',
             fontsize=5.8, ha='center', color=INK)
    ax2.text(4.80, 4.05, r'$\kappa \equiv KVN_A/\Omega$;  $\xi \equiv M/(KVN_A)$',
             fontsize=6.2, ha='center', color=INK)
    ax2.text(4.80, 3.05, 'exact for any M and Ω (equation 1)',
             fontsize=5.6, color=GREY, ha='center', va='center')
    ax2.add_patch(FancyBboxPatch((0.60, 0.45), 8.40, 2.20,
                  boxstyle='round,pad=0.02,rounding_size=0.22',
                  fc=BOX_FC, ec=BOX_EC, lw=1.5, zorder=3))
    ax2.text(4.80, 2.10, r'$\langle B \rangle \, K V = \langle (M-B)\,(\Omega-B) \rangle$',
             fontsize=9.0, ha='center', va='center', color=NAVY, zorder=4)
    ax2.text(4.80, 1.05, r'$= \frac{M\,\Omega}{KV}\,\frac{Q(M-1,\,\Omega-1)}{Q(M,\,\Omega)}$',
             fontsize=7.5, ha='center', va='center', color=NAVY, zorder=4)

    # ======================================== STEP 3 - partition allocation
    badge(ax3, 0.72, 12.45, 3)
    ax3.text(1.32, 12.45, 'Partition allocation', fontsize=FS_T, weight='bold',
             color=NAVY, ha='left', va='center')
    ax3.add_patch(FancyBboxPatch((0.55, 7.10), 2.40, 2.80,
                  boxstyle='round,pad=0.02,rounding_size=0.18',
                  fc='white', ec=NAVY, lw=1.2, zorder=1))
    for yy in (7.95, 9.05):
        for xx in (1.15, 1.75, 2.35):
            draw_site(ax3, xx, yy, 0.30, lw=1.2)
    mol(ax3, 1.75, 8.52, 0.15, motif_deg=180, seed=7)
    ax3.text(0.55, 6.62, r'analog ($N = 1$):', fontsize=5.3,
             ha='left', va='center', color=GREY)
    ax3.text(0.55, 6.05, r'all $\Omega$ sites', fontsize=5.3,
             ha='left', va='center', color=INK)
    ax3.text(0.55, 5.48, 'in one unit', fontsize=5.3, ha='left',
             va='center', color=INK)
    ax3.add_patch(FancyArrowPatch((3.35, 8.60), (4.85, 8.60), arrowstyle='-|>',
                                  mutation_scale=11, color=NAVY, lw=1.6, zorder=5))
    ax3.text(4.10, 9.05, 'partition', fontsize=5.8, ha='center', va='center',
             color=GREY)
    xw0, yw0, dw = 5.05, 7.40, 0.66
    on = {1, 2, 3, 15}
    for i in range(4):
        for j in range(4):
            xw, yw = xw0 + i*dw + dw/2, yw0 + j*dw + dw/2
            pos = (j*4 + i) in on
            ax3.add_patch(Circle((xw, yw), 0.25,
                                 fc=('#D99400' if pos else 'white'),
                                 ec=NAVY, lw=1.0, zorder=3))
    ax3.text(8.50, 6.62, r'digital ($N > 1$):', fontsize=5.3,
             ha='right', va='center', color=GREY)
    ax3.text(8.50, 6.05, r'$\beta = \Omega/N$ sites', fontsize=5.3,
             ha='right', va='center', color=INK)
    ax3.text(8.50, 5.48, 'per partition', fontsize=5.3, ha='right',
             va='center', color=INK)
    bx, by, br = 7.62, 10.95, 0.84
    ax3.add_patch(Circle((bx, by), br, fc='#D99400', ec=NAVY, lw=1.3, zorder=4))
    for ang, occ in [(0, False), (90, True), (180, False), (270, False)]:
        a = np.radians(ang)
        yx, yy = bx + 0.34*np.cos(a), by + 0.34*np.sin(a)
        draw_site(ax3, yx, yy, 0.28, lw=1.2, z=5)
        if occ:
            mol(ax3, yx, yy + 0.32, 0.13, motif_deg=180, seed=9, z=6)
    ax3.plot([7.06, 7.28], [9.98, 10.22], color='#8A97A5', lw=0.9, zorder=3)
    ax3.plot([7.22, 7.35], [9.82, 10.02], color='#8A97A5', lw=0.9, zorder=3)
    ax3.text(6.62, 10.95, 'one\npartition', fontsize=5.8, ha='right',
             va='center', color=GREY, linespacing=1.2)
    ax3.text(4.70, 4.50, 'C bound complexes are dealt to partitions',
             fontsize=5.6, ha='center', va='center', color=INK)
    ax3.text(4.70, 3.80, 'multivariate hypergeometric', fontsize=7.2,
             ha='center', va='center', color=NAVY)
    ax3.text(4.70, 3.00, r'positive if $\geq$ 1 of its $\beta$ sites is occupied',
             fontsize=5.6, ha='center', va='center', color=GREY)
    # glyph key lives here (the coloured circles are in THIS panel)
    ax3.add_patch(Circle((2.50, 1.30), 0.17, fc='#D99400', ec=NAVY, lw=0.9,
                         zorder=3))
    ax3.text(2.77, 1.30, 'positive', fontsize=FS_C, ha='left', va='center',
             color=INK)
    ax3.add_patch(Circle((5.45, 1.30), 0.17, fc='white', ec=NAVY, lw=0.9,
                         zorder=3))
    ax3.text(5.72, 1.30, 'empty', fontsize=FS_C, ha='left', va='center',
             color=INK)

    # ================================================ STEP 4 - observation
    badge(ax4, 0.72, 12.45, 4)
    ax4.text(1.32, 12.45, 'Observation', fontsize=FS_T, weight='bold',
             color=NAVY, ha='left', va='center')
    ax4.text(4.80, 10.90, 'digital: count positive partitions',
             fontsize=5.8, ha='center', va='center', color=GREY)
    ax4.text(4.80, 10.15, r'$Z \sim \mathrm{Binomial}(N,\, P_{\mathrm{pos}})$',
             fontsize=8.8, ha='center', va='center', color=NAVY)
    ax4.text(4.80, 9.35, 'Type 1: binary count + Poisson inversion',
             fontsize=5.6, ha='center', va='center', color=GREY)
    ax4.text(4.80, 8.80, 'Type 2: direct per-molecule counting',
             fontsize=5.6, ha='center', va='center', color=GREY)
    ax4.plot([1.20, 8.40], [8.15, 8.15], color='#C9D2DC', lw=1.0, zorder=1)
    ax4.text(4.80, 7.50, r'$y = D + (A - D)\,p + \varepsilon$',
             fontsize=8.8, ha='center', va='center', color=NAVY)
    ax4.text(4.80, 6.85, 'analog: single intensity',
             fontsize=5.6, ha='center', va='center', color=GREY)
    ax4.text(4.80, 5.90, r'$\kappa \to 0$:  $M = -N\,\ln(1 - N_{\mathrm{pos}}/N)$',
             fontsize=7.6, ha='center', va='center', color=NAVY)
    ax4.text(4.80, 5.20, 'calibration-free absolute count (dPCR)', fontsize=5.6,
             ha='center', va='center', color=GREY)
    ax4.add_patch(FancyBboxPatch((1.20, 0.55), 7.20, 1.80,
                  boxstyle='round,pad=0.02,rounding_size=0.22',
                  fc=BOX_FC, ec=BOX_EC, lw=1.5, zorder=3))
    ax4.text(4.80, 1.45, r'$M = \kappa\,\Omega\,\xi$',
             fontsize=12, weight='bold', ha='center', va='center', color=NAVY,
             zorder=4)
    ax4.text(4.80, 2.72, 'absolute count (equation 4)',
             fontsize=5.6, color=GREY, ha='center', va='center')

    # ----------------------------------------------- chevrons between cards
    fig.canvas.draw()
    for a, b in ((ax1, ax2), (ax2, ax3), (ax3, ax4)):
        pa, pb = a.get_position(), b.get_position()
        gx = (pa.x1 + pb.x0) / 2
        gy, gh, gw = 0.52, 0.060, 0.010
        fig.patches.append(Polygon(
            [(gx - gw, gy + gh), (gx + gw*0.6, gy), (gx - gw, gy - gh),
             (gx - gw*0.4, gy)],
            closed=True, transform=fig.transFigure, fc=NAVY, ec='none', zorder=8))

    import os as _os
    out = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)),
                        'Fig1a_4steps.svg')
    fig.savefig(out, facecolor='white')
    fig.savefig(out.replace('.svg', '_600dpi.png'), dpi=600, facecolor='white')
    return _fig_to_image(fig)
# =============================================================== panel b
def render_b():
    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-
    """
    Fig 1b redesign: "Unified mapping of bioanalytical technologies in the
    TCS operating space."

    v8.2 — kappa_1 ~= 0.1 / kappa_2 ~= 10 threshold lines added
      (text: asymptotic-analysis thresholds); bars anchor on them.

    Current design — one platform colour and the same uniform #F5F7FA
      background as the panel-a cards; no background grid, binary
      counting/sensing regions, diagonal separator or regime legend.
      The kappa_1 and kappa_2 threshold lines and all platform-specific
      kappa, beta, N and operating-range annotations are retained.

    Fact base (manuscript only — no invented coordinates):
      Caption: x = kappa (binding regime parameter), y = beta (capture sites
               per partition); bubble size = total partitions N.
      Table 1: qPCR   N=1,   beta->inf,   kappa->0 ("strong depletion, yet
                                                    calibration-bound by
                                                    amplification efficiency")
               dPCR   N>1e4, beta->inf,   kappa->0 ("calibration free")
               PICO   N>1e4, beta->inf,   TCS: kappa = 0.1-0.6
               dELISA N>1e5, 1<=beta<inf, finite kappa ("digital readout but
                                                        calibration required")
               DigitISA/Singulex N=1, beta->inf, finite kappa ("direct
                                  p-access, calibration-dependent")
               ELISA 5PL regime  N=1, beta->inf, 1 <~ kappa <~ 10 (BG=1)
               ELISA 4PL regime  N=1, beta->inf, kappa >~ 10 (B=1, Langmuir)
               AAI    N=1, beta~1e5-1e7, kappa >= 100 ("depletion < 1%")

    Range encoding (all ranges verbatim from Table 1 / main text):
      thick bar   = finite numeric range (PICO kappa 0.1-0.6; 5PL 1<~k<~10;
                    4PL kappa>~10 with right arrow; AAI region edges;
                    dELISA kappa 0.579-8.23 = the two analysed Simoa fits:
                    sBetaG 2010 kappa=0.579 [0.347,0.664], IL-6 kappa=8.23
                    [7.98,8.58])
      thin arrow  = open range to a limit (kappa->0: qPCR/dPCR left arrows;
                    beta->inf: up arrows on qPCR, dPCR, PICO, DigitISA, 5PL,
                    4PL; AAI right edge open + arrow = kappa unbounded;
                    DigitISA kappa>0.1 lower-bound bar + right arrow, from
                    the main-text kappa-spectrum sentence "Singulex/DigitISA
                    at kappa > 0.1"; Table 1 gives no upper value)
      AAI shaded region = kappa>=100 x beta~1e5-1e7 (Table 1 row). Note:
                    kappa*beta = KVN_A/N is fixed for a given chemistry, so
                    kappa and beta cannot both diverge; the region, not a
                    diagonal, is the Table-1-consistent shape.
      dELISA 1<=beta<inf: down arrow (beta extends to 1, below the axis
                    floor) + exact range in its sub-label.
    Positions on the log-log map are schematic but regime membership and all
    annotated ranges follow Table 1 / the caption verbatim.

    Bubble-size legend states N only ("N = 1 (one partition)" / "N > 1
    (partitioned)") — it must NOT equate N=1 with analog: Table 1 gives
    DigitISA/Singulex N=1 yet the main text says it "supplies digital"
    (p = Z/Omega, "direct p-access"). DigitISA's digital identity is carried
    by its own sub-label instead.

    Output: Fig1b_重设计.svg (fonts as paths) + Fig1b_重设计.png (600 dpi)
    """
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    from matplotlib.patches import FancyBboxPatch, Rectangle

    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Liberation Sans', 'Arial', 'DejaVu Sans']
    plt.rcParams['svg.fonttype'] = 'path'
    plt.rcParams['mathtext.fontset'] = 'custom'
    plt.rcParams['mathtext.it'] = 'Liberation Sans:italic'
    plt.rcParams['mathtext.rm'] = 'Liberation Sans'

    NAVY   = '#1F2D4E'
    # One platform colour: the map does not define binary counting/sensing classes.
    C_TGT  = '#0072B2'
    C_ON   = '#E69F00'   # unused
    C_AB   = '#0072B2'
    GREY   = '#5A6577'
    CARD_FC = '#F5F7FA'  # same light blue-grey as the panel-a cards
    BOX_FC, BOX_EC = '#FFF4D6', '#D9BC5F'

    def lighten(hexcol, f=0.45):
        r, g, b = mcolors.to_rgb(hexcol)
        return (r + (1 - r) * f, g + (1 - g) * f, b + (1 - b) * f)

    # ------------------------------------------------------------- canvas
    fig, ax = plt.subplots(figsize=(5.8, 5.0), facecolor=CARD_FC)
    ax.set_facecolor(CARD_FC)
    ax.set_xlim(-2.55, 2.75)          # log10 kappa
    ax.set_ylim(3.3, 10.5)            # log10 beta (bottom strip: this work)
    ax.set_xticks(range(-2, 3))
    ax.set_xticklabels(['$10^{-2}$', '$10^{-1}$', '$10^{0}$', '$10^{1}$',
                        '$10^{2}$'], fontsize=8)
    ax.set_yticks(range(5, 11))
    ax.set_yticklabels([f'$10^{{{d}}}$' for d in range(5, 11)], fontsize=8)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_color(NAVY); ax.spines[s].set_linewidth(1.3)
    # minor ticks at log positions (axes are linear over log10 coordinates),
    # matched in style to the other main-text figures (outward, 3.5/2.0 pt)
    _xmin = [np.log10(_m) + _d for _d in range(-2, 2) for _m in range(2, 10)]
    _ymin = [np.log10(_m) + _d for _d in range(5, 10) for _m in range(2, 10)]
    ax.set_xticks(_xmin, minor=True)
    ax.set_yticks(_ymin, minor=True)
    ax.tick_params(colors=NAVY, length=3.5, width=1.0, direction='out')
    ax.tick_params(which='minor', colors=NAVY, length=2.0, width=0.8,
                   direction='out')
    ax.grid(False)

    # kappa thresholds kappa_1 ~= 0.1, kappa_2 ~= 10 (text, asymptotic
    # analysis): vertical dotted lines; Table-1 bars anchor on them
    # (PICO/DigitISA start at kappa_1; 5PL ends at kappa_2; 4PL starts at it)
    for xx in (-1.0, 1.0):
        ax.plot([xx, xx], [3.3, 10.5], ls=(0, (2, 2)), color=NAVY, lw=1.1,
                alpha=0.75, zorder=1)
    ax.text(-1.0, 10.62, r'$\kappa_1 \approx 0.1$', fontsize=7, color=NAVY,
            ha='center', weight='bold')
    ax.text(1.0, 10.62, r'$\kappa_2 \approx 10$', fontsize=7, color=NAVY,
            ha='center', weight='bold')

    # ------------------------------------------------------------ helpers
    def bubble(X, Y, s, color, z=5):
        """s in points^2; pseudo-3D via concentric overlays."""
        ax.scatter([X], [Y], s=s * 1.28, c=[lighten(color, 0.72)],
                   edgecolors='none', zorder=z - 0.1)          # halo
        ax.scatter([X], [Y], s=s, c=[color], edgecolors=NAVY,
                   linewidths=1.2, zorder=z)
        ax.scatter([X - 0.05], [Y + 0.11], s=s * 0.16, c=['white'],
                   edgecolors='none', alpha=0.65, zorder=z + 1)  # highlight

    def kappa_bar(x1, x2, Y, color):
        """finite numeric range (Table 1)."""
        ax.plot([x1, x2], [Y, Y], color=color, lw=3, alpha=0.45,
                solid_capstyle='round', zorder=4)

    def range_arrow(x1, y1, x2, y2, color):
        """open range to a limit (kappa->0, kappa->inf, beta->inf, beta->1)."""
        ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle='-|>', color=color, lw=2.0,
                                    alpha=0.55, mutation_scale=8))

    S_SMALL, S_BIG, S_HUGE = 260, 850, 1350
    _KO = dict(fc=CARD_FC, ec='none', alpha=0.92, pad=0.15)  # line knockout

    # --- counting side (Table 1 rows) ---
    bubble(-2.15, 9.10, S_SMALL, C_TGT)                       # qPCR
    kappa_bar(-2.48, -2.34, 9.10, C_TGT)                      # kappa -> 0
    range_arrow(-2.40, 9.10, -2.55, 9.10, C_TGT)
    range_arrow(-2.15, 9.38, -2.15, 9.67, C_TGT)              # beta -> inf
    ax.text(-1.86, 9.52, 'qPCR', fontsize=8.5, weight='bold', color=NAVY)
    ax.text(-1.84, 9.34, r'$\kappa \rightarrow 0$ · calibration-bound',
            fontsize=6, color=GREY, va='center', bbox=_KO)

    bubble(-1.70, 8.50, S_BIG, C_TGT)                         # dPCR
    kappa_bar(-2.20, -1.90, 8.50, C_TGT)
    range_arrow(-2.15, 8.50, -2.33, 8.50, C_TGT)
    range_arrow(-1.70, 8.90, -1.70, 9.20, C_TGT)              # beta -> inf
    ax.text(-1.38, 8.50, 'dPCR', fontsize=8.5, weight='bold', color=NAVY,
            va='center')
    ax.text(-1.82, 7.92, r'$\kappa \rightarrow 0$ · calibration free',
            fontsize=6, color=GREY, ha='center')

    bubble(-0.72, 8.45, S_BIG, C_TGT)                         # PICO
    kappa_bar(-1.00, -0.22, 8.45, C_TGT)                      # kappa = 0.1-0.6
    range_arrow(-0.72, 8.85, -0.72, 9.15, C_TGT)              # beta -> inf
    ax.text(-0.60, 9.32, 'PICO', fontsize=8.5, weight='bold', color=NAVY,
            ha='center')
    ax.text(-0.95, 7.81, r'TCS: $\kappa = 0.1$–$0.6$', fontsize=6,
            color=GREY, ha='center', bbox=_KO)

    bubble(0.00, 8.02, S_SMALL, C_TGT)                        # DigitISA/Singulex
    kappa_bar(-1.00, 0.40, 8.02, C_TGT)                       # kappa > 0.1
    range_arrow(0.40, 8.02, 0.70, 8.02, C_TGT)
    range_arrow(0.00, 8.295, 0.00, 8.60, C_TGT)               # beta -> inf
    ax.text(0.10, 8.86, 'DigitISA / Singulex', fontsize=8.5, weight='bold',
            color=NAVY, ha='center')
    ax.text(0.00, 7.64, 'digital · direct p-access', fontsize=6, color=GREY,
            ha='center')
    ax.text(0.00, 7.46, 'calibration-dependent', fontsize=6, color=GREY,
            ha='center')

    # --- counting side, cont.: dELISA (digital, finite kappa) ---
    bubble(0.00, 6.30, S_HUGE, C_TGT)                          # dELISA
    kappa_bar(-0.237, 0.915, 6.30, C_TGT)                      # kappa 0.579-8.23
    range_arrow(0.00, 5.76, 0.00, 4.62, C_TGT)                 # beta down to 1
    ax.text(-0.74, 6.48, 'Simoa (dELISA)', fontsize=8.5, weight='bold',
            color=NAVY, ha='right', bbox=_KO)
    ax.text(-0.74, 6.22, r'finite $\kappa$ · $1 \leq \beta < \infty$',
            fontsize=6, color=GREY, ha='right', bbox=_KO)
    ax.text(-0.74, 5.96, 'calibration required', fontsize=6, color=GREY,
            ha='right', bbox=_KO)

    # --- sensing side: ELISA 5PL regime (analog, intermediate kappa) ---
    bubble(0.90, 7.30, S_SMALL, C_AB)                         # ELISA 5PL regime
    kappa_bar(0.00, 1.00, 7.30, C_AB)                         # 1 <~ kappa <~ 10
    range_arrow(0.90, 7.575, 0.90, 7.875, C_AB)               # beta -> inf
    ax.text(1.24, 7.68, 'ELISA (5PL regime)', fontsize=8.5, weight='bold',
            color=NAVY, va='center')
    ax.text(1.24, 7.40, r'$BG = 1$ · $1 \lesssim \kappa \lesssim 10$',
            fontsize=6, color=GREY, va='center')

    # --- sensing side ---
    bubble(1.35, 6.70, S_SMALL, C_AB)                         # ELISA 4PL regime
    kappa_bar(1.00, 1.70, 6.70, C_AB)                         # kappa >~ 10
    range_arrow(1.72, 6.70, 1.94, 6.70, C_AB)
    range_arrow(1.35, 6.975, 1.35, 7.275, C_AB)               # beta -> inf
    ax.text(1.36, 6.13, 'ELISA (4PL regime)', fontsize=8.5, weight='bold',
            color=NAVY, ha='center', bbox=_KO)
    ax.text(1.36, 5.85, r'$B = 1$ · $\kappa \gtrsim 10$ · Langmuir limit',
            fontsize=6, color=GREY, ha='center', bbox=_KO)

    # AAI region: kappa >= 100 right edge crossed by the kappa->inf arrow x beta ~ 1e5-1e7
    ax.add_patch(Rectangle((2.00, 5.50), 0.65, 1.50, fc=C_AB, ec='none',
                 alpha=0.09, zorder=2))
    for xx, yy, dx, dy in ((2.00, 5.50, 0, 1.50), (2.00, 5.50, 0.65, 0),
                           (2.00, 7.00, 0.65, 0),
                           (2.65, 5.50, 0, 1.50)):
        ax.plot([xx, xx + dx], [yy, yy + dy], color=C_AB, lw=1.2, alpha=0.55,
                ls=(0, (4, 2)), zorder=2)
    bubble(2.30, 6.00, S_SMALL, C_AB)                         # AAI
    range_arrow(2.55, 6.62, 2.74, 6.62, C_AB)                 # kappa -> inf
    ax.text(2.375, 6.68, 'AAI', fontsize=8.5, weight='bold', color=NAVY,
            ha='center')
    ax.text(2.375, 6.40, r'$\beta \sim 10^{5}$–$10^{7}$', fontsize=6,
            color=GREY, ha='center')
    ax.text(2.70, 5.30, r'$\kappa \geq 100$ · depletion $< 1\%$', fontsize=6,
            color=GREY, ha='right')

    # -------------------------------------------- this work: droplet dELISA
    # kappa ~ 3.5e3, beta -> 1 (~1 effective site per bead -> tiny Omega ->
    # kappa huge): drawn AAI-style as a dashed box at the bottom strip,
    # immediately left of the AAI region (which shares kappa >= 100 but has
    # beta ~ 1e5-1e7; the contrast in beta is the point).
    bx0, bx1, by0, by1 = 1.40, 2.10, 3.40, 4.34
    ax.add_patch(Rectangle((bx0, by0), bx1 - bx0, by1 - by0, fc=C_TGT,
                           ec='none', alpha=0.09, zorder=2))
    for xx, yy, dx, dy in ((bx0, by0, 0, by1 - by0), (bx0, by0, bx1 - bx0, 0),
                           (bx0, by1, bx1 - bx0, 0), (bx1, by0, 0, by1 - by0)):
        ax.plot([xx, xx + dx], [yy, yy + dy], color=C_TGT, lw=1.2, alpha=0.55,
                ls=(0, (4, 2)), zorder=2)
    bubble(1.75, 3.87, 1500, C_TGT)      # N ~ 2e6 droplets: biggest ball here
    ax.text(1.24, 4.10, 'droplet dELISA (this work)', fontsize=7,
            weight='bold', color=NAVY, ha='right', bbox=_KO)
    ax.text(1.24, 3.88, r'$\kappa \approx 3.5\times10^3$ · $\beta \to 1$',
            fontsize=6, color=GREY, ha='right', bbox=_KO)
    ax.text(1.24, 3.66, r'$\approx$1 effective site per bead', fontsize=6,
            color=GREY, ha='right', bbox=_KO)

    # -------------------------------------------------------------- legend
    lx, ly = 1.28, 10.28
    ax.scatter([lx], [ly], s=40, c=['white'], edgecolors=NAVY,
               linewidths=1, zorder=7)
    ax.text(lx + 0.14, ly, 'N = 1 (one partition)', fontsize=6.5,
            va='center', color=NAVY)
    ax.scatter([lx], [ly - 0.34], s=150, c=['white'], edgecolors=NAVY,
               linewidths=1, zorder=7)
    ax.text(lx + 0.14, ly - 0.34, 'N > 1 (partitioned)', fontsize=6.5,
            va='center', color=NAVY)
    ax.plot([lx - 0.08, lx + 0.02], [ly - 0.68, ly - 0.68], color=GREY,
            lw=2.2, alpha=0.55, solid_capstyle='round', zorder=7)
    ax.annotate('', xy=(lx + 0.10, ly - 0.68), xytext=(lx + 0.02, ly - 0.68),
                arrowprops=dict(arrowstyle='-|>', color=GREY, lw=2.2,
                                alpha=0.55, mutation_scale=8))
    ax.text(lx + 0.14, ly - 0.68, r'$\kappa$ / $\beta$ operating range',
            fontsize=6.5, va='center', color=NAVY, bbox=_KO)

    # ------------------------------------------------------------ axis labels
    ax.set_xlabel(r'binding regime parameter  $\kappa = KVN_A/\Omega$',
                  fontsize=9, color=NAVY, weight='bold')
    ax.set_ylabel(r'capture sites per partition  $\beta = \Omega/N$',
                  fontsize=9, color=NAVY, weight='bold')
    ax.annotate('', xy=(2.75, 3.3), xytext=(2.45, 3.3),
                arrowprops=dict(arrowstyle='-|>', color=NAVY, lw=1.3))
    ax.annotate('', xy=(-2.55, 10.5), xytext=(-2.55, 10.2),
                arrowprops=dict(arrowstyle='-|>', color=NAVY, lw=1.3))

    fig.tight_layout()
    return _fig_to_image(fig)

# =============================================================== panel c
def render_c():
    # Fig 1c (3D xi surface) portrait re-render for the composite rebalance:
    #   figsize (7.0, 8.6), box_aspect (1,1,1.4), base font 14 pt,
    #   legend above the axes, horizontal colorbar below, x ticks thinned
    #   to [0.01, 1, 100] -> /tmp/fig1c.png (600 dpi, aspect ~0.709)
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    plt.rcParams.update({'font.family': 'sans-serif',
                         'font.sans-serif': ['Arial', 'Liberation Sans', 'DejaVu Sans'],
                         'mathtext.fontset': 'custom',
                         'mathtext.rm': 'Arial',
                         'mathtext.it': 'Arial:italic',
                         'mathtext.bf': 'Arial:bold'})
    FS = 14

    p_vals = np.logspace(-4, np.log10(0.9999), 1000)
    kappa_vals = np.logspace(-2, 3, 1000)
    P, K = np.meshgrid(p_vals, kappa_vals)
    Xi = P/(1-P) + P/K
    logXi, logK, z_P = np.log10(Xi), np.log10(K), P
    kappa_curve = np.logspace(-2, 3, 500)
    p_curve = (2*kappa_curve + 1 - np.sqrt(4*kappa_curve**2 + 1)) / 2
    p_elev = p_curve + 0.008
    idx_text = int(np.argmin(np.abs(kappa_curve - 1.0)))

    fig = plt.figure(figsize=(7.0, 8.6))
    ax = fig.add_subplot(111, projection='3d')
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    for k in ('x', 'y', 'z'):
        getattr(ax, f'{k}axis')._axinfo["grid"]['color'] = (0.9, 0.9, 0.9, 0.5)
    surf = ax.plot_surface(logXi, logK, z_P, cmap='coolwarm', edgecolor='none',
                           alpha=0.9, zorder=1)
    ax.plot(np.zeros_like(kappa_curve), np.log10(kappa_curve), p_elev,
            color='red', lw=5, zorder=15, label=r'$\xi = 1$')
    ax.text(0, 0, p_elev[idx_text]+0.1, r'$\xi = 1$', color='red',
            fontsize=FS, weight='bold', zorder=20)
    # regime threshold curves kappa = 0.1 / 10 on the surface (echo panel b)
    for k0, lab in ((0.1, r'$\kappa = 0.1$'), (10.0, r'$\kappa = 10$')):
        p_t = np.logspace(-4, np.log10(0.999), 400)
        xi_t = p_t / (1 - p_t) + p_t / k0
        ax.plot(np.log10(xi_t), np.full_like(p_t, np.log10(k0)), p_t + 0.008,
                color='white', lw=4.2, zorder=14)
        ax.plot(np.log10(xi_t), np.full_like(p_t, np.log10(k0)), p_t + 0.008,
                color='#1F2D4E', lw=1.8, ls=(0, (4, 2)), zorder=15, label=lab)
    # contour projection on the floor for depth
    ax.contour(logXi, logK, z_P, levels=np.arange(0.1, 1.0, 0.1), zdir='z',
               offset=0, cmap='coolwarm', linewidths=0.7, alpha=0.5, zorder=2)
    ax.set_xlabel(r'$\log_{10}\xi$', fontsize=FS)
    xi_ticks = [0.01, 1, 100]
    ax.set_xticks(np.log10(xi_ticks))
    ax.set_xticklabels([f'{v:g}' for v in xi_ticks])
    ax.set_ylabel(r'$\log_{10}\kappa$', fontsize=FS)
    kappa_ticks = [0.001, 0.01, 0.1, 1, 10, 100, 1000]
    ax.set_yticks(np.log10(kappa_ticks))
    ax.set_yticklabels([f'{v:g}' for v in kappa_ticks])
    ax.set_zticks(np.arange(0, 1.01, 0.2))
    ax.set_zticklabels([f'{t:.1f}' for t in np.arange(0, 1.01, 0.2)])
    ax.set_zlim(0, 1)
    ax.tick_params(axis='x', labelsize=FS)
    ax.tick_params(axis='y', labelsize=FS)
    ax.tick_params(axis='z', labelsize=FS)
    ax.set_box_aspect((1, 1, 1.4))
    ax.view_init(elev=30, azim=-110)
    ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.0), fontsize=FS,
              framealpha=1.0)
    cbar = fig.colorbar(surf, ax=ax, shrink=0.5, aspect=12, pad=0.12,
                        orientation='horizontal')
    cbar.set_label(r'$p$', fontsize=FS)
    cbar.ax.tick_params(labelsize=FS)
    return _fig_to_image(fig)

# ================================================================ composite
img_a = render_a()
print(f'panel a: {img_a.size[0]}x{img_a.size[1]} px, aspect {img_a.size[0]/img_a.size[1]:.4f}')
img_b = render_b()
print(f'panel b: {img_b.size[0]}x{img_b.size[1]} px, aspect {img_b.size[0]/img_b.size[1]:.4f}')
img_c = render_c()
print(f'panel c: {img_c.size[0]}x{img_c.size[1]} px, aspect {img_c.size[0]/img_c.size[1]:.4f}')

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Liberation Sans', 'Arial', 'DejaVu Sans']   # restore after render_c
MM = 1 / 25.4

def _place(fig, page_w, page_h, im, x_mm, y_mm, w_mm):
    """Place an image with its native aspect; (x_mm, y_mm) = bottom-left."""
    w_px, h_px = im.size
    h_mm = w_mm * h_px / w_px
    ax = fig.add_axes([x_mm / page_w, y_mm / page_h,
                       w_mm / page_w, h_mm / page_h])
    ax.imshow(np.asarray(im), aspect='auto', interpolation='lanczos')
    ax.axis('off')
    return h_mm

PW = 183.0
wa, wb, wc = 182.6, 119.0, 61.5
ha = wa * img_a.size[1] / img_a.size[0]      # heights from measured aspects
hb = wb * img_b.size[1] / img_b.size[0]
hc = wc * img_c.size[1] / img_c.size[0]
GUT, TOP, BOT = 4.0, 4.0, 1.5   # TOP hosts the panel letter + title band
PH = TOP + ha + GUT + hb + BOT
print(f'abc page: {PW:.1f} x {PH:.1f} mm (limit 183 x 247)')

fig = plt.figure(figsize=(PW * MM, PH * MM), facecolor='white')
ya = PH - TOP - ha
_place(fig, PW, PH, img_a, 0.2, ya, wa)
y2 = ya - GUT                      # top of the b/c row
yb = y2 - hb
_place(fig, PW, PH, img_b, 0.2, yb, wb)
xc = 121.4
yc = y2 - hc                       # top-aligned with b (letters on one line)
_place(fig, PW, PH, img_c, xc, yc, wc)

# panel labels: Nature-style bold lowercase letters (descriptions live in the caption)
panel_labels = [(1.2, PH - 0.4, 'top', 'a'),
                (1.2, y2 + 0.5, 'bottom', 'b'),
                (xc + 1.0, y2 + 0.5, 'bottom', 'c')]
for lx, ly, va, lab in panel_labels:
    fig.text(lx / PW, ly / PH, lab, fontsize=10, weight='bold',
             color='black', ha='left', va=va)

OUT_DIR = os.path.dirname(os.path.abspath(__file__))
fig.savefig(os.path.join(OUT_DIR, 'Fig1_abc.svg'), dpi=300,
            facecolor='white')
fig.savefig(os.path.join(OUT_DIR, 'Fig1_abc.pdf'), facecolor='white')
fig.savefig(os.path.join(OUT_DIR, 'Fig1_abc.png'), dpi=600,
            facecolor='white')
print('saved abc composite')