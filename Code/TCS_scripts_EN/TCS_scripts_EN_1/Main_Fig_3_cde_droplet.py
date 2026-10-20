#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main_Fig_3_cde_droplet.py — Figure 3 panels c-e: home-built droplet-based
digital immunoassay at the deep Langmuir end of the kappa spectrum.
(Rebuilt from Main_Fig_4_droplet.py after the droplet panels were merged
into Figure 3; panel c -> chip, d -> calibration, e -> clinical.)

Layout (v11): two rows. Panel c is a full-width top strip (chip cross-section
schematic with the OCA top view, plus the bright-field and fluorescence
micrographs); panels d and e are large square axes below, matched in
displayed size to the Fig. 3a panel.

Panels:
  c  platform: chip cross-section schematic (glass / 50 um OCA spacer with
     laser-cut channel / glass, two 0.9 mm holes as inlet and outlet) with the
     channel top view drawn from the CAD geometry (Drawing1.dxf: 20 x 70 mm
     tape, 10 x 35 mm channel, 6 mm wells), plus bright-field (bead monolayer)
     and two-channel fluorescence micrographs of the 10,000 fg/mL standard
     (green: all beads, FITC channel; red: positives). Micrograph assets are
     built by build_fig4_panel_a_assets.py (fig4a_brightfield.png /
     fig4a_fluorescence.png).
  d  equilibrium response: mean occupancy lambda vs nominal IL-6 concentration,
     saturation hyperbola fitted by the R2 binomial counting likelihood with
     the background b estimated jointly (SI S2f), counting-statistics
     LoB/LoD/LoQ.
  e  clinical cross-platform validation: droplet vs Roche Cobas 411 (n=12) and
     vs Simoa HD-X (n=10 in range), log-log, identity line.

Raw data: positive/negative droplet-bead counts recorded in the preprint
(ChemRxiv 10.26434/chemrxiv-2024-19rj6, v7) and re-analysed here with the TCS
framework. All fits are recomputed inline; the script is self-contained.

Output: Main_Fig_3_cde_droplet.png (600 dpi) / Main_Fig_3_cde_droplet.svg
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager as _fm
from scipy.optimize import brentq, minimize
from scipy.stats import gaussian_kde, beta as _beta
from matplotlib.transforms import blended_transform_factory
from matplotlib.lines import Line2D
from PIL import Image

from tcs_figstyle import apply_style, OI

_here = os.path.dirname(os.path.abspath(__file__))
_fonts = os.path.join(_here, 'fonts')
if os.path.isdir(_fonts):
    for _f in os.listdir(_fonts):
        if _f.lower().endswith(('.ttf', '.otf')):
            _fm.fontManager.addfont(os.path.join(_fonts, _f))
apply_style(base_size=12)

# ---------------- raw data ----------------
NA = 6.02214076e23
MW = 20900.0          # g/mol, recombinant human IL-6
VS = 160e-6           # L reaction volume
N_NOM = 2e6           # nominal bead number
VD_NOM = VS / N_NOM   # 8e-11 L

# (nominal fg/mL, positive, negative)
STD = [
    (1e6, 26675, 66273),
    (1e5, 5655, 90581),
    (1e4, 776, 114680),
    (1e3, 125, 139732),
    (100, 56, 193471),
    (10, 35, 133020),
    (1, 18, 161295),
]
BLANK = (21, 154409)  # 0 fg/mL

# (total beads observed, positive beads, Cobas pg/mL, Simoa pg/mL or None if exceeding range)
CLIN = [
    (130564, 44046, 1865.00, None),
    (136445, 32237, 1037.00, None),
    (115343, 13018, 558.80, 161.13),
    (165529, 880, 9.93, 3.27),
    (118828, 17596, 532.80, 184.92),
    (172737, 2577, 66.82, 29.43),
    (142637, 3898, 123.60, 52.67),
    (197965, 10288, 264.70, 90.35),
    (170980, 5435, 150.90, 60.28),
    (195295, 1869, 45.37, 17.56),
    (203242, 13412, 327.60, 106.65),
    (171355, 2880, 95.15, 34.34),
]


def lam(k, n):
    return -np.log(1 - k / n)


lam_bg = lam(BLANK[0], BLANK[0] + BLANK[1])

Cs = np.array([s[0] for s in STD])
Ls = np.array([lam(k, k + ng) for _, k, ng in STD])

def hyper(C, lmax, c50):
    return lmax * C / (c50 + C)

# ---- R2 counting-likelihood fit -----------------------------------------
# Joint binomial likelihood of all seven standards plus the blank, with the
# background occupancy b fitted together with (lambda_max, C50) instead of
# being pre-subtracted (SI S2f: joint estimation of the background is both
# safe and recommended). Model: lambda_tot(C) = b + lmax*C/(C50 + C),
# positives binomial with p = 1 - exp(-lambda_tot); fit in log space.
_k_all = np.array([s[1] for s in STD] + [BLANK[0]], dtype=float)
_n_all = np.array([s[1] + s[2] for s in STD] + [BLANK[0] + BLANK[1]],
                  dtype=float)
_C_all = np.array([s[0] for s in STD] + [0.0])

def _nll(th):
    a, c, d = th
    lam_tot = np.exp(d) + np.exp(a) * _C_all / (np.exp(c) + _C_all)
    p = np.clip(1.0 - np.exp(-lam_tot), 1e-15, 1.0 - 1e-15)
    return -np.sum(_k_all * np.log(p) + (_n_all - _k_all) * np.log(1.0 - p))

_res = minimize(_nll, np.log([0.7, 1e6, lam_bg]), method='Nelder-Mead',
                options=dict(xatol=1e-10, fatol=1e-9, maxiter=40000))
_res = minimize(_nll, _res.x, method='Nelder-Mead',
                options=dict(xatol=1e-10, fatol=1e-9, maxiter=40000))
th_hat = _res.x
lmax, c50, b_fit = np.exp(th_hat)

# expected Fisher information in log space (exact for the binomial model);
# its inverse is the joint covariance, marginal over the fitted background.
lam_tot = b_fit + hyper(_C_all, lmax, c50)
p_hat = 1.0 - np.exp(-lam_tot)
_lsp = hyper(_C_all, lmax, c50)
_d = np.vstack([_lsp, -_lsp * c50 / (c50 + _C_all),
                np.full_like(_C_all, b_fit)])
_w = _n_all * np.exp(-2.0 * lam_tot) / (p_hat * (1.0 - p_hat))
_F = (_d * _w) @ _d.T
cov_ln = np.linalg.inv(_F)

# goodness of fit over all seven standards (joint model incl. background)
_mod = b_fit + hyper(Cs, lmax, c50)
r2 = 1 - np.sum((Ls - _mod) ** 2) / np.sum((Ls - Ls.mean()) ** 2)
slope = lmax / c50

p_b = BLANK[0] / (BLANK[0] + BLANK[1])
sig_b = np.sqrt(p_b * (1 - p_b) / (BLANK[0] + BLANK[1])) / (1 - p_b)
LoB = 1.645 * sig_b / slope
LoD = LoB + 1.645 * np.sqrt(2) * sig_b / slope

def cv_total(C):
    lsp = slope * C
    n = 1.5e5
    p = 1 - np.exp(-(lsp + b_fit))
    var = p / (n * (1 - p)) + sig_b ** 2   # sample counting + one blank subtraction
    return np.sqrt(var) / lsp

LoQ = brentq(lambda C: cv_total(C) - 0.20, LoD * 1.01, 1e5)

kappa = VD_NOM / (slope / (1e-12 / MW * NA)) - 1

# ---- parametric bootstrap: sampling distributions of LoB / LoD / LoQ ----
# blank positives are binomial (only 21 of ~1.54e5 beads -> ~15% relative
# uncertainty on the blank occupancy); the hyperbola slope uncertainty is
# propagated from the log-space Fisher covariance of the R2 likelihood fit.
rng_bs = np.random.default_rng(42)
N_BS = 20000
n_blank = BLANK[0] + BLANK[1]
k_bs = rng_bs.binomial(n_blank, p_b, N_BS)
p_bs = np.maximum(k_bs / n_blank, 1e-12)
sig_bs = np.sqrt(p_bs * (1 - p_bs) / n_blank) / (1 - p_bs)
_th = np.exp(rng_bs.multivariate_normal(th_hat[:2], cov_ln[:2, :2], N_BS))
sl_bs = _th[:, 0] / _th[:, 1]

LoB_bs = 1.645 * sig_bs / sl_bs
LoD_bs = LoB_bs + 1.645 * np.sqrt(2) * sig_bs / sl_bs


def _loq_one(sl, sb, lb):
    def f(C):
        lsp = sl * C
        p = 1 - np.exp(-(lsp + lb))
        return np.sqrt(p / (1.5e5 * (1 - p)) + sb ** 2) / lsp - 0.20
    try:
        return brentq(f, 1.0, 1e6, maxiter=100)
    except Exception:
        return np.nan


LoQ_bs = np.array([_loq_one(s, b, b_fit) for s, b in
                   zip(sl_bs[::2], sig_bs[::2])])
LoQ_bs = LoQ_bs[np.isfinite(LoQ_bs) & (LoQ_bs > 0)]

# clinical inversion
dp = np.array([(lambda l: (l - b_fit) * c50 / (lmax - (l - b_fit)) / 1000)(lam(k, n))
               for n, k, _, _ in CLIN])
cb = np.array([c[2] for c in CLIN])
sm = np.array([c[3] if c[3] is not None else np.nan for c in CLIN])
msk = ~np.isnan(sm)

def r2_of(x, y):
    A = np.vstack([x, np.ones_like(x)]).T
    m, b = np.linalg.lstsq(A, y, rcond=None)[0]
    return 1 - np.sum((y - (m * x + b)) ** 2) / np.sum((y - y.mean()) ** 2)

r2_cb = r2_of(cb, dp)
r2_sm = r2_of(sm[msk], dp[msk])

# ---------------- figure ----------------
fig = plt.figure(figsize=(10.4, 7.8))


def plabel(x, y, s):
    fig.text(x, y, s, fontsize=17, fontweight='bold', va='bottom', ha='left')


# --- panel c: platform schematic + micrographs (top strip) ---
axs = fig.add_axes([0.045, 0.715, 0.455, 0.26])
axs.set_xlim(0, 1.66); axs.set_ylim(-0.09, 0.98); axs.axis('off')
GLASS = '#d7e8f5'; GE = '#8aa8bd'; OCAC = '#f2c14e'; LIQ = '#e9f5fc'
chx = (0.14, 0.80)                      # chamber extent# bottom glass
axs.add_patch(plt.Rectangle((0.06, 0.04), 0.88, 0.24, fc=GLASS, ec=GE, lw=1.0))
# OCA spacer flanking the chamber
axs.add_patch(plt.Rectangle((0.06, 0.28), chx[0] - 0.06, 0.12, fc=OCAC,
                            ec='0.55', lw=0.8))
axs.add_patch(plt.Rectangle((chx[1], 0.28), 0.94 - chx[1], 0.12, fc=OCAC,
                            ec='0.55', lw=0.8))
# chamber liquid
axs.add_patch(plt.Rectangle((chx[0], 0.28), chx[1] - chx[0], 0.12, fc=LIQ,
                            ec='none'))
# top glass with two holes
for xa, xb in [(0.06, 0.24), (0.29, 0.65), (0.70, 0.94)]:
    axs.add_patch(plt.Rectangle((xa, 0.40), xb - xa, 0.24, fc=GLASS, ec=GE,
                                lw=1.0))
# beads in the chamber (right of the label zone)
rng = np.random.default_rng(7)
bx = rng.uniform(0.44, chx[1] - 0.04, 14)
by = rng.uniform(0.295, 0.375, 14)
for i, (x, y) in enumerate(zip(bx, by)):
    axs.add_patch(plt.Circle((x, y), 0.0085,
                             fc=OI['vermillion'] if i in (3, 9) else '0.35',
                             ec='none'))
# inlet / outlet arrows through the holes
axs.annotate('', xy=(0.265, 0.42), xytext=(0.265, 0.74),
             arrowprops=dict(arrowstyle='-|>', color='0.2', lw=1.4))
axs.text(0.265, 0.78, 'inlet\n' + r'($\varnothing$0.9 mm)', fontsize=8.5,
         ha='center', va='bottom', color='0.15')
axs.annotate('', xy=(0.675, 0.74), xytext=(0.675, 0.42),
             arrowprops=dict(arrowstyle='-|>', color='0.2', lw=1.4))
axs.text(0.675, 0.78, 'outlet', fontsize=8.5, ha='center', va='bottom',
         color='0.15')
# microscope objective
axs.add_patch(plt.Polygon([(0.415, -0.04), (0.525, -0.04), (0.49, 0.04),
                           (0.45, 0.04)], fc='#c3c9cf', ec='0.45', lw=0.8))
axs.text(0.545, -0.01, '10x objective', fontsize=8.5, ha='left', va='center',
         color='0.3')
# layer labels
axs.text(0.30, 0.16, 'glass slide', fontsize=9, color='#39586e', va='center')
axs.text(0.47, 0.52, 'top glass', fontsize=9, color='#39586e', va='center',
         ha='center')
axs.text(0.155, 0.33, 'sample\nchamber', fontsize=8, color='#4a6f88',
         va='center')

# top view of the OCA cutout (Drawing1.dxf geometry, mm, rotated 90 deg)
ox, oy, smm = 1.06, 0.34, 0.0076
def T(x_mm, y_mm):
    return ox + x_mm * smm, oy + (y_mm - 10) * smm
axs.add_patch(plt.Rectangle(T(0, 0), 70 * smm, 20 * smm, fc='#f5f7f9',
                            ec='0.45', lw=0.9))
axs.add_patch(plt.Rectangle(T(17.5, 5), 35 * smm, 10 * smm, fc=LIQ, ec='none'))
axs.add_patch(plt.Polygon([T(12.5, 7), T(12.5, 13), T(17.5, 15), T(17.5, 5)],
                          fc=LIQ, ec='none'))
axs.add_patch(plt.Polygon([T(57.5, 13), T(57.5, 7), T(52.5, 5), T(52.5, 15)],
                          fc=LIQ, ec='none'))
for cx in (12.5, 57.5):
    axs.add_patch(plt.Circle(T(cx, 10), 3 * smm, fc=LIQ, ec='0.25', lw=0.9))
axs.plot([T(17.5, 5)[0], T(52.5, 5)[0]], [T(17.5, 5)[1]] * 2, color='0.25',
         lw=0.9)
axs.plot([T(17.5, 15)[0], T(52.5, 15)[0]], [T(17.5, 15)[1]] * 2, color='0.25',
         lw=0.9)
for p1, p2 in [((12.5, 13), (17.5, 15)), ((12.5, 7), (17.5, 5)),
               ((57.5, 13), (52.5, 15)), ((57.5, 7), (52.5, 5))]:
    axs.plot([T(*p1)[0], T(*p2)[0]], [T(*p1)[1], T(*p2)[1]], color='0.25',
             lw=0.9)
for cx in (12.5, 57.5):
    axs.add_patch(plt.Circle(T(cx, 10), 0.007, fc='white', ec='0.2', lw=0.8))
axs.text(T(35, 10)[0], oy + 0.115, 'OCA cutout (top view)', fontsize=8,
         ha='center', va='bottom', color='0.35')
axs.text(T(35, 10)[0], oy - 0.125, 'tape 20 x 70 mm\nchannel 35 x 10 mm',
         fontsize=7.5, ha='center', va='top', color='0.35')
# OCA swatch legend
axs.add_patch(plt.Rectangle((1.06, 0.015), 0.028, 0.038, fc=OCAC, ec='0.55',
                            lw=0.8))
axs.text(1.10, 0.034, 'OCA spacer, 50 µm', fontsize=7.5, color='#7a5b12',
         ha='left', va='center')

# micrographs
axb1 = fig.add_axes([0.555, 0.715, 0.195, 0.26])
axb1.imshow(Image.open(os.path.join(_here, 'fig4a_brightfield.png')))
axb1.axis('off')
axb1.set_title('bright field: bead monolayer', fontsize=9.5, pad=3)

axf1 = fig.add_axes([0.770, 0.715, 0.195, 0.26])
axf1.imshow(Image.open(os.path.join(_here, 'fig4a_fluorescence.png')))
axf1.plot([800 - 96, 800 - 28], [800 - 34, 800 - 34], color='white', lw=3.2,
          solid_capstyle='butt')
axf1.text(800 - 62, 800 - 66, '50 µm', fontsize=9.5, color='white',
          ha='center')
axf1.text(14, 800 - 82, 'beads', fontsize=8.5, color='#8cff8c', ha='left',
          weight='bold', bbox=dict(fc='black', alpha=0.45, pad=1.6, ec='none'))
axf1.text(14, 800 - 26, 'positives', fontsize=8.5, color='#ff7a6e', ha='left',
          weight='bold', bbox=dict(fc='black', alpha=0.45, pad=1.6, ec='none'))
axf1.axis('off')
axf1.set_title('fluorescence, 10$^4$ fg/mL IL-6', fontsize=9.5, pad=3)
plabel(0.006, 0.972, 'c')

# --- panel d: standard curve ---
axb = fig.add_axes([0.062, 0.085, 0.43, 0.582])
Cf = np.logspace(-0.3, 6.7, 400)
axb.plot(Cf, b_fit + hyper(Cf, lmax, c50), color=OI['blue'], lw=2,
         label='saturation hyperbola fit')
axb.axhline(b_fit, color='0.35', ls='--', lw=1.2,
            label='fitted background $b$')
axb.axhline(b_fit + lmax, color=OI['vermillion'], ls='-.', lw=1.2,
            label='$\\lambda_{max}$ = %.3f' % lmax)
axb.axvline(c50, color='0.6', ls=':', lw=1.0)
axb.text(4.5e6, 1.28, '$C_{50}\\approx K\\approx$50 pM', fontsize=10,
         color='0.35', ha='right', va='top')
# threshold point estimates (dotted) + bootstrap sampling distributions
# (shaded humps, log-x density, drawn in a bottom band in axes coordinates)
_tr = blended_transform_factory(axb.transData, axb.transAxes)

def _hump(samples, color):
    xs = np.log10(samples)
    kde = gaussian_kde(xs)
    lo, hi = np.percentile(xs, [0.3, 99.7])
    grid = np.linspace(lo, hi, 200)
    dens = kde(grid)
    y0 = 0.015
    ytop = y0 + 0.115 * dens / dens.max()
    axb.fill_between(10 ** grid, y0, ytop, transform=_tr, color=color,
                     alpha=0.35, lw=0, zorder=2)
    axb.plot(10 ** grid, ytop, transform=_tr, color=color, lw=1.0, zorder=3)

# 95% interval bound pairs + boxed interval legend: identical convention to
# the Simoa panels (2)/(4) (peaks unlabelled; entries in a lower-right box)
for _bs, _c, _ls in [(LoB_bs, OI['green'], ':'),
                     (LoD_bs, OI['orange'], '--'),
                     (LoQ_bs, OI['purple'], '-.')]:
    for _b in np.percentile(_bs, [2.5, 97.5]):
        axb.axvline(_b, color=_c, linestyle=_ls, linewidth=1.2, alpha=0.8)
_hump(LoB_bs, OI['green'])
_hump(LoD_bs, OI['orange'])
_hump(LoQ_bs, OI['purple'])
_qb, _qd, _qq = (np.percentile(_a, [2.5, 97.5])
                 for _a in (LoB_bs, LoD_bs, LoQ_bs))
_leg2 = axb.legend(
    [Line2D([0], [0], color=OI['green'], lw=1.5),
     Line2D([0], [0], color=OI['orange'], lw=1.5),
     Line2D([0], [0], color=OI['purple'], lw=1.5)],
    ['LoB [%.0f, %.0f] fg/mL' % tuple(_qb),
     'LoD [%.0f, %.0f] fg/mL' % tuple(_qd),
     'LoQ [%.0f, %.0f] fg/mL' % tuple(_qq)],
    loc='lower right', bbox_to_anchor=(0.985, 0.02), fontsize=9,
    framealpha=0.95, ncol=1, borderaxespad=0.6, handlelength=1.6,
    labelspacing=0.35, handletextpad=0.6)
axb.add_artist(_leg2)
# standards with 95% Clopper-Pearson counting intervals (same convention as
# the Simoa panels): each point's occupancy CI follows from its own counts
_Ks = np.array([s[1] for s in STD], dtype=float)
_Ns = np.array([s[1] + s[2] for s in STD], dtype=float)
_lam_lo = np.array([-np.log(1 - _beta.ppf(0.025, k, n - k + 1))
                    for k, n in zip(_Ks, _Ns)])
_lam_hi = np.array([-np.log(1 - _beta.ppf(0.975, k + 1, n - k))
                    for k, n in zip(_Ks, _Ns)])
_yerr = [Ls - _lam_lo, _lam_hi - Ls]
axb.errorbar(Cs[:4], Ls[:4], yerr=[_yerr[0][:4], _yerr[1][:4]], fmt='o',
             color='0.2', ms=6, ecolor='0.45', elinewidth=1.0, capsize=2.5,
             label='standards', zorder=5)
axb.errorbar(Cs[4:], Ls[4:], yerr=[_yerr[0][4:], _yerr[1][4:]], fmt='o',
             color=OI['orange'], ms=6, ecolor=OI['orange'], elinewidth=1.0,
             capsize=2.5, alpha=0.9, label='below LoD', zorder=5)
axb.set_xscale('log'); axb.set_yscale('log')
axb.set_xlim(0.6, 5e6); axb.set_ylim(8e-5, 1.6)
axb.set_xlabel('nominal IL-6 concentration (fg/mL)')
axb.set_ylabel('mean occupancy $\\lambda$')
_r2lab = ('$R^2>0.9999$' if r2 >= 0.9999 else '$R^2$=%.4f' % r2)
axb.text(0.97, 0.30, '$\\lambda_{sp}=\\lambda_{max}\\,C/(C_{50}+C)$\n' + _r2lab,
         transform=axb.transAxes, fontsize=10.5, va='top', ha='right')
axb.legend(loc='upper left', bbox_to_anchor=(0.02, 0.86), fontsize=9.5,
           handlelength=1.4, framealpha=1.0)
plabel(0.012, 0.678, 'd')

# --- panel d inset: joint uncertainty of (lambda_max, C50) and the
# parameter-degeneracy direction. The log-space Fisher covariance of the R2
# likelihood fit (marginal over the fitted background b) is drawn as a 95%
# confidence ellipse in log-log space; the major axis is the direction along
# which the two parameters can be traded against each other (only the ratio
# lambda_max/C50, the linear-range slope, is tightly pinned).
# Placed in the empty left-middle band of panel d: clear of the density
# humps below, the LoB/LoD/LoQ interval lines on the right, the upper-left
# legend above. Tick labels are suppressed so nothing spills onto the main
# axes; the caption explains the inset.
axin = fig.add_axes([0.077, 0.236, 0.101, 0.128])
_sub = rng_bs.integers(0, len(_th), 1400)
axin.plot(_th[_sub, 1], _th[_sub, 0], '.', color=OI['blue'], ms=1.2,
          alpha=0.20, rasterized=True)
_cl = cov_ln[:2, :2]                           # covariance in ln space
_eval, _evec = np.linalg.eigh(_cl)
_o = np.argsort(_eval)[::-1]
_eval, _evec = _eval[_o], _evec[:, _o]
_tt = np.linspace(0, 2 * np.pi, 200)
_r = np.sqrt(5.991)                          # chi2(df=2) at 95%
_ex = _r * np.sqrt(_eval[0]) * np.cos(_tt)
_ey = _r * np.sqrt(_eval[1]) * np.sin(_tt)
_el = th_hat[:2][:, None] + _evec @ np.vstack([_ex, _ey])
axin.plot(np.exp(_el[1]), np.exp(_el[0]), color=OI['vermillion'], lw=1.4)
_c0 = th_hat[:2]
_maj = _evec[:, 0] * _r * np.sqrt(_eval[0])
_t1 = np.exp(_c0 + _maj); _t2 = np.exp(_c0 - _maj)
axin.annotate('', xy=(_t1[1], _t1[0]), xytext=(_t2[1], _t2[0]),
              arrowprops=dict(arrowstyle='<|-|>', color='0.3', lw=1.0,
                              shrinkA=0, shrinkB=0))
axin.plot([c50], [lmax], '+', color='0.15', ms=6.5, mew=1.5)
axin.set_xscale('log'); axin.set_yscale('log')
axin.tick_params(which='both', direction='in', labelbottom=False,
                 labelleft=False, length=2)
for _s in axin.spines.values():
    _s.set_linewidth(0.8)
# interior mini-labels (tick labels suppressed: they would spill onto the
# main panel's tick labels)
_xl, _xh = axin.get_xlim(); _yl, _yh = axin.get_ylim()
axin.text(10 ** (np.log10(_xh) - 0.06 * np.log10(_xh / _xl)),
          10 ** (np.log10(_yl) + 0.05 * np.log10(_yh / _yl)),
          '$C_{50}$', fontsize=7.5, ha='right', va='bottom', color='0.25')
axin.text(10 ** (np.log10(_xl) + 0.05 * np.log10(_xh / _xl)),
          10 ** (np.log10(_yh) - 0.04 * np.log10(_yh / _yl)),
          '$\\lambda_{max}$', fontsize=7.5, ha='left', va='top',
          rotation=90, color='0.25')
_axr = np.sqrt(_eval[0] / _eval[1])
axin.text(0.95, 0.30, 'degeneracy\ndirection', transform=axin.transAxes,
          fontsize=6.5, color='0.3', ha='right', va='center',
          linespacing=1.15)
print('degeneracy inset: major-axis slope d(ln lmax)/d(ln c50) = %.3f, '
      'axis ratio %.1f:1' % (_evec[0, 0] / _evec[1, 0], _axr))

# --- panel e: clinical validation ---
axd = fig.add_axes([0.552, 0.085, 0.43, 0.582])
lims = [2, 4000]
axd.plot(lims, lims, '--', color='0.55', lw=1.2, label='identity')
axd.plot(cb, dp, 'o', color=OI['blue'], ms=6.5,
         label='vs Roche Cobas 411 ($n$=12)')
axd.plot(sm[msk], dp[msk], 's', color=OI['vermillion'], ms=6.5,
         label='vs Simoa HD-X ($n$=10)')
axd.set_xscale('log'); axd.set_yscale('log')
axd.set_xlim(lims); axd.set_ylim(lims)
axd.set_xlabel('reference concentration (pg/mL)')
axd.set_ylabel('droplet platform (pg/mL)')
axd.text(0.03, 0.97, '$R^2$ = %.3f (Cobas)\n$R^2$ = %.3f (Simoa)' % (r2_cb, r2_sm),
         transform=axd.transAxes, fontsize=10.5, va='top')
axd.legend(loc='lower right', fontsize=9.5, handlelength=1.4)
plabel(0.502, 0.678, 'e')


out_png = os.path.join(_here, 'Main_Fig_3_cde_droplet.png')
out_svg = os.path.join(_here, 'Main_Fig_3_cde_droplet.svg')
fig.savefig(out_png, dpi=600)
fig.savefig(out_svg)
print('saved:', out_png)
print('saved:', out_svg)
print('R2-likelihood fit: lmax = %.4f [%.4f, %.4f], c50 = %.4e [%.3e, %.3e], '
      'b = %.3e (measured blank %.3e), R2 = %.5f' % (
          lmax, *np.percentile(_th[:, 0], [2.5, 97.5]),
          c50, *np.percentile(_th[:, 1], [2.5, 97.5]), b_fit, lam_bg, r2))
_kap_bs = VD_NOM / (sl_bs / (1e-12 / MW * NA)) - 1
print('kappa = %.0f [95%% CI %.0f, %.0f]' % (
    kappa, *np.percentile(_kap_bs, [2.5, 97.5])))
print('kappa =', round(kappa), '| LoB/LoD/LoQ =', round(LoB), round(LoD), round(LoQ))
for _nm, _a in [('LoB', LoB_bs), ('LoD', LoD_bs), ('LoQ', LoQ_bs)]:
    _q = np.percentile(_a, [2.5, 50, 97.5])
    print('%s bootstrap median %.0f [95%% CI %.0f, %.0f]' % (_nm, _q[1], _q[0], _q[2]))
print('R2 Cobas = %.4f | R2 Simoa = %.4f' % (r2_cb, r2_sm))

# ============================================================================
# Workbook export: raw counts, fit results and statistics (main-text Fig. 3c-e)
# ============================================================================
import pandas as _pd
from scipy.stats import beta as _betad

def _cp_lam(_k, _n):
    """Clopper-Pearson 95% interval for the positive fraction, mapped to lambda."""
    _lo = _betad.ppf(0.025, _k, _n - _k + 1) if _k > 0 else 0.0
    _hi = _betad.ppf(0.975, _k + 1, _n - _k) if _k < _n else 1.0
    return -np.log(1 - _lo), -np.log(1 - _hi)


_xlsx = os.path.join(_here, 'SI_Table_6.xlsx')
with _pd.ExcelWriter(_xlsx) as _xw:
    _pd.DataFrame({'SI Table 6: droplet digital immunoassay workbook': [
        'Raw droplet-bead counts, TCS/R2 fit results and counting-statistics',
        'performance bounds of the home-built platform (main-text Fig. 3c-e).',
        '',
        'Sheets:',
        '  Standards_raw     - positive/negative bead counts per standard',
        '  Fit_parameters    - R2 binomial-likelihood fit (lambda_max, C50, b),',
        '                      kappa, slope, R2, Fisher covariance, bootstrap CI',
        '  LoB_LoD_LoQ       - point estimates and parametric-bootstrap 95% CI',
        '  Clinical_12       - twelve clinical plasma samples: raw counts,',
        '                      droplet-platform back-calculation, Roche Cobas 411',
        '                      and Simoa HD-X reference values',
        '  Bootstrap_samples - raw bootstrap draws for LoB / LoD / LoQ (fg/mL)',
        '',
        'Model: lambda_tot(C) = b + lambda_max*C/(C50 + C); positives binomial',
        'with p = 1 - exp(-lambda_tot); fitted jointly in log space (SI S2f).',
        'kappa = V_d / (slope in molecule units) - 1, V_d = 160 uL / 2e6 beads.',
    ]}).to_excel(_xw, sheet_name='README', index=False)

    _rows = []
    for _c, _k, _n in ([(s[0], s[1], s[1] + s[2]) for s in STD]
                       + [(0.0, BLANK[0], BLANK[0] + BLANK[1])]):
        _lo, _hi = _cp_lam(_k, _n)
        _rows.append({'nominal_IL6_fg_per_mL': _c,
                      'positive_beads': int(_k), 'negative_beads': int(_n - _k),
                      'total_beads': int(_n), 'lambda_observed': lam(_k, _n),
                      'lambda_CI95_lower': _lo, 'lambda_CI95_upper': _hi})
    _pd.DataFrame(_rows).to_excel(_xw, sheet_name='Standards_raw', index=False)

    _b_lo, _b_hi = np.exp(th_hat[2] + np.array([-1.96, 1.96])
                          * np.sqrt(cov_ln[2, 2]))
    _rho = np.corrcoef(np.log(_th[:, 0]), np.log(_th[:, 1]))[0, 1]
    _fit = _pd.DataFrame([
        ('lambda_max', lmax, *np.percentile(_th[:, 0], [2.5, 97.5]), '-'),
        ('C50_fg_per_mL', c50, *np.percentile(_th[:, 1], [2.5, 97.5]), '-'),
        ('background_b', b_fit, _b_lo, _b_hi,
         'measured blank lambda = %.3e' % lam_bg),
        ('linear_slope_per_fg_mL', slope, *np.percentile(sl_bs, [2.5, 97.5]),
         '-'),
        ('kappa', kappa, *np.percentile(_kap_bs, [2.5, 97.5]),
         'platform binding-regime parameter'),
        ('R_squared_fit', r2, np.nan, np.nan, 'all seven standards + blank'),
        ('R_squared_vs_Cobas', r2_cb, np.nan, np.nan, 'n = 12'),
        ('R_squared_vs_Simoa', r2_sm, np.nan, np.nan, 'n = 10 in range'),
        ('degeneracy_major_slope', _evec[0, 0] / _evec[1, 0], np.nan, np.nan,
         'd(ln lambda_max)/d(ln C50) along the 95% ellipse major axis'),
        ('degeneracy_axis_ratio', _axr, np.nan, np.nan,
         'major/minor axis ratio of the log-space Fisher 95% ellipse'),
        ('bootstrap_corr_ln_lmax_ln_C50', _rho, np.nan, np.nan,
         'Pearson correlation of the log-space bootstrap samples'),
    ], columns=['quantity', 'estimate', 'CI95_lower', 'CI95_upper', 'note'])
    _fit.to_excel(_xw, sheet_name='Fit_parameters', index=False)
    _cov = _pd.DataFrame(cov_ln, index=['ln_lambda_max', 'ln_C50', 'ln_b'],
                         columns=['ln_lambda_max', 'ln_C50', 'ln_b'])
    _cov.to_excel(_xw, sheet_name='Fit_parameters', index=True,
                  startrow=len(_fit) + 3)
    _ws = _xw.sheets['Fit_parameters']
    _ws.cell(row=len(_fit) + 3, column=1,
             value='log-space Fisher covariance (R2 likelihood)')

    _rows = []
    for _nm, _pt, _a in [('LoB', LoB, LoB_bs), ('LoD', LoD, LoD_bs),
                         ('LoQ', LoQ, LoQ_bs)]:
        _lo, _md, _hi = np.percentile(_a, [2.5, 50, 97.5])
        _rows.append((_nm, _pt, _md, _lo, _hi))
    _pd.DataFrame(
        _rows, columns=['bound', 'point_estimate_fg_per_mL', 'bootstrap_median',
                        'bootstrap_CI95_lower', 'bootstrap_CI95_upper']
    ).to_excel(_xw, sheet_name='LoB_LoD_LoQ', index=False)

    _cl = []
    for (_n, _k, _cb, _sm), _d in zip(CLIN, dp):
        _lo, _hi = _cp_lam(_k, _n)
        _cl.append({'total_beads': int(_n), 'positive_beads': int(_k),
                    'lambda_observed': lam(_k, _n),
                    'lambda_CI95_lower': _lo, 'lambda_CI95_upper': _hi,
                    'droplet_platform_pg_per_mL': _d,
                    'Roche_Cobas411_pg_per_mL': _cb,
                    'Simoa_HDX_pg_per_mL': _sm if _sm is not None else ''})
    _pd.DataFrame(_cl).to_excel(_xw, sheet_name='Clinical_12', index=False)

    _bs = _pd.concat([_pd.Series(_th[:, 0], name='lambda_max'),
                      _pd.Series(_th[:, 1], name='C50_fg_per_mL'),
                      _pd.Series(LoB_bs, name='LoB_fg_per_mL'),
                      _pd.Series(LoD_bs, name='LoD_fg_per_mL'),
                      _pd.Series(LoQ_bs, name='LoQ_fg_per_mL')], axis=1)
    _bs.to_excel(_xw, sheet_name='Bootstrap_samples', index=False)
print('bootstrap correlation rho(ln lambda_max, ln C50) = %.3f' % _rho)
print('saved:', _xlsx)
