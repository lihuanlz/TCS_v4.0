# -*- coding: utf-8 -*-
"""
Supplementary Figure ED6 (droplet morphometry) and accompanying Excel workbook
for the droplet-based digital immunoassay platform (manuscript Fig. 3c-e,
Supplementary Information S12).

Provenance of the embedded data
-------------------------------
Bright-field image of the bead monolayer acquired in the observation chip in
the same experiments as the assay (Fig. 3c). Magnetic beads (nominal diameter
3 um) were marked by hand; each marked bead served as an anchor. For every
anchor the droplet boundary was measured by bead-anchored ray-casting: along
72 rays per bead the droplet rim was located as the 55% half-height point of
the inner flank of the dark rim, and the droplet radius was taken from the
inscribed boundary. Anchors whose measured droplet radius is not at least
1.25x the bead radius are bare beads (no droplet) and are excluded; two
borderline objects were excluded on visual inspection (status column).

Image scale: the drawn 50 um scale bar measures 333 px
(0.15015 um/px), independently verified against the nominal 3 um bead
diameter (140 hand-marked beads, measured median below).

Contents
--------
_SCALE_BEADS : (140, 3) hand-marked scale beads - cx, cy (px), bead radius (px)
_MEAS        : (131, 5) droplet anchors    - cx, cy (px), droplet radius (px),
               bead radius (px), status (1 = bead-bearing droplet, 0 = bare bead)

Outputs
-------
ED_Fig_6.png / .svg : (a) droplet diameter distribution,
                      (b) concentration bias for a volume partition
ED_Fig_6_droplet_morphometry.xlsx : README, Summary_stats, Droplet_measurements,
                                    Scale_beads, Bias_curve
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy import stats
from tcs_figstyle import apply_style, OI, save_fig, panel_label

print("=" * 80)
print("ED Fig. 6 - droplet morphometry of the droplet-based digital immunoassay")
print("=" * 80)

# ----------------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------------
SCALE_BAR_UM = 50.0          # drawn scale bar length (um)
SCALE_BAR_PX = 333.0         # measured length of the scale bar (px)
UMPX = SCALE_BAR_UM / SCALE_BAR_PX   # um per pixel
NOMINAL_BEAD_UM = 3.0        # nominal magnetic bead diameter (um)
VOLUME_ASSAY_UL = 50.0       # emulsion aqueous phase per assay (uL)

# ----------------------------------------------------------------------------
# Embedded measurements (see provenance above)
# ----------------------------------------------------------------------------
_SCALE_BEADS = np.array([
    (507.8, 268.2, 8.20),
    (612.7, 254.2, 7.00),
    (659.0, 270.3, 8.00),
    (684.9, 300.1, 8.50),
    (687.9, 263.1, 6.50),
    (742.3, 360.4, 9.70),
    (692.6, 354.8, 8.60),
    (648.4, 391.4, 6.20),
    (241.6, 383.5, 10.30),
    (161.2, 414.5, 8.20),
    (196.7, 444.4, 8.70),
    (282.1, 424.3, 7.70),
    (298.9, 382.7, 7.30),
    (331.3, 384.4, 7.00),
    (273.2, 331.3, 8.20),
    (306.7, 267.6, 7.80),
    (252.0, 284.9, 7.90),
    (406.3, 401.0, 9.40),
    (475.2, 391.2, 8.70),
    (519.1, 410.8, 9.30),
    (531.3, 374.6, 8.50),
    (511.8, 320.4, 9.40),
    (483.5, 311.6, 8.40),
    (446.8, 367.3, 8.50),
    (408.7, 338.5, 8.20),
    (390.7, 289.7, 9.20),
    (415.1, 242.8, 8.70),
    (375.5, 440.1, 9.40),
    (387.7, 372.7, 9.10),
    (296.4, 459.6, 8.60),
    (249.0, 461.6, 9.00),
    (266.6, 483.6, 11.20),
    (186.5, 482.1, 10.60),
    (106.9, 484.5, 10.20),
    (34.6, 363.9, 10.30),
    (30.2, 323.8, 8.70),
    (70.8, 264.3, 8.50),
    (187.0, 212.5, 10.20),
    (265.6, 197.8, 8.10),
    (382.8, 211.0, 9.50),
    (711.1, 401.0, 9.90),
    (686.5, 418.4, 9.40),
    (642.0, 418.4, 8.70),
    (604.5, 464.3, 9.50),
    (685.9, 456.5, 10.70),
    (724.0, 435.2, 10.10),
    (767.3, 426.2, 9.50),
    (794.4, 463.7, 8.40),
    (833.8, 450.1, 10.00),
    (855.8, 433.9, 9.80),
    (888.0, 413.9, 8.70),
    (912.6, 433.9, 8.80),
    (875.8, 457.2, 9.30),
    (930.7, 471.4, 9.20),
    (943.6, 538.6, 11.00),
    (888.0, 539.9, 11.00),
    (837.0, 546.3, 10.10),
    (811.8, 520.5, 11.30),
    (778.9, 540.5, 10.40),
    (808.0, 562.5, 10.80),
    (755.0, 530.8, 9.40),
    (751.8, 474.0, 10.80),
    (806.7, 425.5, 8.70),
    (788.6, 399.1, 9.80),
    (796.3, 351.3, 12.00),
    (853.2, 355.1, 11.00),
    (876.4, 386.1, 10.40),
    (884.8, 340.3, 9.00),
    (922.9, 370.6, 12.70),
    (855.8, 309.3, 11.00),
    (853.2, 277.6, 11.40),
    (813.1, 298.3, 11.60),
    (780.8, 276.3, 11.00),
    (893.0, 725.6, 10.40),
    (923.8, 725.6, 11.60),
    (905.0, 753.8, 9.40),
    (957.9, 729.0, 9.20),
    (943.4, 778.5, 13.50),
    (884.5, 782.8, 11.10),
    (783.7, 723.9, 10.30),
    (831.5, 739.2, 12.60),
    (711.9, 757.2, 14.20),
    (477.9, 614.5, 12.30),
    (439.4, 607.7, 12.10),
    (420.6, 641.8, 10.30),
    (426.6, 676.0, 12.10),
    (378.8, 647.8, 11.10),
    (395.9, 621.3, 12.30),
    (512.0, 442.0, 10.90),
    (556.5, 459.0, 10.00),
    (477.9, 471.0, 13.30),
    (442.8, 481.3, 13.80),
    (397.6, 527.4, 12.30),
    (310.4, 600.0, 11.60),
    (285.7, 619.6, 11.50),
    (934.0, 332.6, 10.90),
    (188.0, 608.3, 8.90),
    (228.1, 590.5, 11.60),
    (283.1, 574.2, 9.40),
    (325.4, 531.8, 10.60),
    (269.0, 667.8, 10.30),
    (271.2, 701.2, 10.00),
    (237.8, 715.3, 8.70),
    (208.8, 706.4, 10.10),
    (182.0, 660.3, 9.70),
    (136.0, 642.5, 9.50),
    (133.0, 583.1, 10.00),
    (156.8, 544.5, 9.50),
    (212.5, 553.4, 10.00),
    (99.6, 723.5, 10.60),
    (78.1, 723.5, 7.80),
    (77.3, 755.4, 8.90),
    (150.1, 755.4, 9.70),
    (180.6, 733.9, 8.00),
    (168.7, 705.6, 8.90),
    (35.0, 656.6, 9.70),
    (46.9, 577.1, 8.70),
    (124.1, 452.3, 10.70),
    (656.8, 589.3, 9.70),
    (749.1, 579.5, 11.10),
    (720.7, 609.0, 12.60),
    (777.6, 647.3, 10.60),
    (844.4, 641.4, 9.80),
    (815.9, 650.2, 8.80),
    (895.5, 620.7, 12.00),
    (974.1, 583.4, 11.10),
    (969.2, 634.5, 8.50),
    (905.3, 575.6, 9.80),
    (818.9, 778.9, 10.00),
    (785.5, 783.8, 9.80),
    (742.3, 787.7, 11.00),
    (744.2, 758.3, 10.60),
    (571.3, 729.8, 10.40),
    (581.2, 677.7, 10.60),
    (578.2, 614.8, 12.80),
    (617.5, 549.0, 13.20),
    (626.3, 596.2, 10.60),
    (525.2, 586.4, 9.80),
    (541.9, 636.5, 10.60),
    (488.8, 750.4, 11.80)
])                            # (140, 3): cx, cy, r_bead (px)

_MEAS = np.array([
    (506.0, 278.0, 11.599, 8.20, 1),
    (609.0, 260.0, 11.437, 7.00, 1),
    (654.0, 281.0, 11.889, 8.00, 1),
    (695.0, 312.0, 17.099, 8.50, 1),
    (691.0, 251.0, 11.359, 6.50, 1),
    (740.0, 372.0, 19.055, 9.70, 1),
    (684.0, 366.0, 21.541, 8.60, 1),
    (660.0, 403.0, 10.283, 6.20, 1),
    (243.0, 392.0, 24.206, 10.30, 1),
    (165.0, 402.0, 15.544, 8.20, 1),
    (186.0, 443.0, 12.051, 8.70, 1),
    (274.0, 429.0, 10.642, 7.70, 1),
    (291.0, 386.0, 11.658, 7.30, 1),
    (319.0, 372.0, 8.711, 7.00, 0),
    (263.0, 343.0, 11.478, 8.20, 1),
    (295.0, 278.0, 15.963, 7.80, 1),
    (240.0, 283.0, 25.513, 7.90, 1),
    (397.0, 391.0, 5.124, 9.40, 0),
    (472.0, 379.0, 8.737, 8.70, 0),
    (507.0, 409.0, 22.421, 9.30, 1),
    (519.0, 362.0, 7.234, 8.50, 0),
    (442.0, 364.0, 19.730, 8.50, 1),
    (398.0, 350.0, 7.943, 8.20, 0),
    (380.0, 301.0, 21.548, 9.20, 1),
    (427.0, 250.0, 30.159, 8.70, 1),
    (364.0, 443.0, 19.430, 9.40, 1),
    (375.0, 373.0, 18.603, 9.10, 1),
    (292.0, 463.0, 12.795, 8.60, 1),
    (237.0, 453.0, 16.723, 9.00, 1),
    (267.0, 495.0, 15.886, 11.20, 1),
    (192.0, 494.0, 20.294, 10.60, 1),
    (101.0, 496.0, 8.506, 10.20, 0),
    (22.0, 351.0, 13.468, 10.30, 1),
    (37.0, 311.0, 22.789, 8.70, 1),
    (58.0, 271.0, 4.577, 8.50, 0),
    (197.0, 224.0, 5.124, 10.20, 0),
    (261.0, 189.0, 22.830, 8.10, 1),
    (708.0, 413.0, 10.770, 9.90, 0),
    (678.0, 423.0, 9.437, 9.40, 0),
    (640.0, 428.0, 13.503, 8.70, 1),
    (596.0, 452.0, 16.142, 9.50, 1),
    (676.0, 465.0, 13.375, 10.70, 1),
    (715.0, 444.0, 11.068, 10.10, 0),
    (757.0, 432.0, 11.146, 9.50, 0),
    (782.0, 473.0, 6.918, 8.40, 0),
    (828.0, 462.0, 17.039, 10.00, 1),
    (850.0, 424.0, 6.303, 9.80, 0),
    (900.0, 405.0, 11.786, 8.70, 1),
    (913.0, 441.0, 12.077, 8.80, 1),
    (880.0, 469.0, 17.287, 9.30, 1),
    (934.0, 467.0, 11.172, 9.20, 0),
    (951.0, 526.0, 9.685, 11.00, 0),
    (878.0, 545.0, 18.269, 11.00, 1),
    (841.0, 534.0, 16.912, 10.10, 1),
    (776.0, 552.0, 18.703, 10.40, 1),
    (743.0, 535.0, 14.605, 9.40, 1),
    (742.0, 480.0, 12.452, 10.80, 0),
    (802.0, 431.0, 12.547, 8.70, 1),
    (776.0, 387.0, 5.261, 9.80, 0),
    (800.0, 353.0, 32.816, 12.00, 1),
    (845.0, 367.0, 8.839, 11.00, 0),
    (865.0, 398.0, 14.348, 10.40, 1),
    (896.0, 330.0, 16.228, 9.00, 1),
    (910.0, 371.0, 17.099, 12.70, 1),
    (850.0, 321.0, 15.886, 11.00, 1),
    (843.0, 273.0, 12.051, 11.40, 0),
    (805.0, 310.0, 7.259, 11.60, 0),
    (774.0, 284.0, 10.889, 11.00, 0),
    (895.0, 723.0, 5.808, 10.40, 0),
    (935.0, 713.0, 26.937, 11.60, 1),
    (897.0, 751.0, 10.505, 9.40, 0),
    (945.0, 741.0, 15.416, 9.20, 1),
    (932.0, 784.0, 17.876, 13.50, 1),
    (882.0, 788.0, 9.839, 11.10, 0),
    (783.0, 735.0, 21.020, 10.30, 1),
    (832.0, 748.0, 17.202, 12.60, 1),
    (712.0, 746.0, 19.448, 14.20, 1),
    (478.0, 626.0, 8.096, 12.30, 0),
    (438.0, 614.0, 15.246, 12.10, 1),
    (417.0, 653.0, 7.703, 10.30, 0),
    (428.0, 688.0, 12.846, 12.10, 0),
    (385.0, 659.0, 19.731, 11.10, 1),
    (392.0, 610.0, 10.761, 12.30, 0),
    (508.0, 454.0, 9.771, 10.90, 0),
    (551.0, 470.0, 13.067, 10.00, 1),
    (465.0, 472.0, 5.397, 13.30, 0),
    (439.0, 493.0, 25.095, 13.80, 1),
    (301.0, 590.0, 3.620, 11.60, 0),
    (279.0, 630.0, 14.519, 11.50, 1),
    (925.0, 342.0, 7.601, 10.90, 0),
    (180.0, 620.0, 21.096, 8.90, 1),
    (236.0, 602.0, 16.066, 11.60, 1),
    (275.0, 586.0, 13.256, 9.40, 1),
    (313.0, 535.0, 24.633, 10.60, 1),
    (257.0, 669.0, 11.342, 10.30, 0),
    (241.0, 727.0, 18.834, 8.70, 1),
    (199.0, 701.0, 14.451, 10.10, 1),
    (177.0, 669.0, 10.881, 9.70, 0),
    (137.0, 654.0, 15.886, 9.50, 1),
    (132.0, 595.0, 18.227, 10.00, 1),
    (151.0, 554.0, 12.759, 9.50, 1),
    (201.0, 562.0, 18.653, 10.00, 1),
    (87.0, 735.0, 7.053, 10.60, 0),
    (86.0, 735.0, 7.292, 7.80, 0),
    (65.0, 767.0, 8.916, 8.90, 0),
    (138.0, 753.0, 19.252, 9.70, 1),
    (186.0, 745.0, 9.215, 8.00, 0),
    (156.0, 704.0, 16.057, 8.90, 1),
    (24.0, 653.0, 12.298, 9.70, 1),
    (46.0, 565.0, 13.067, 8.70, 1),
    (667.0, 577.0, 19.055, 9.70, 1),
    (741.0, 589.0, 10.633, 11.10, 0),
    (712.0, 620.0, 16.723, 12.60, 1),
    (772.0, 652.0, 10.104, 10.60, 0),
    (832.0, 631.0, 9.215, 9.80, 0),
    (815.0, 638.0, 23.318, 8.80, 1),
    (883.0, 625.0, 8.711, 12.00, 0),
    (974.0, 594.0, 12.555, 11.10, 0),
    (966.0, 622.0, 11.530, 8.50, 1),
    (893.0, 587.0, 16.680, 9.80, 1),
    (822.0, 790.0, 9.924, 10.00, 0),
    (780.0, 787.0, 12.316, 9.80, 1),
    (730.0, 775.0, 8.361, 11.00, 0),
    (732.0, 770.0, 6.046, 10.60, 0),
    (588.0, 689.0, 20.499, 10.60, 1),
    (575.0, 620.0, 18.014, 12.80, 1),
    (607.0, 561.0, 25.965, 13.20, 1),
    (634.0, 608.0, 10.326, 10.60, 0),
    (519.0, 598.0, 20.431, 9.80, 1),
    (530.0, 648.0, 14.605, 10.60, 1),
    (484.0, 762.0, 20.498, 11.80, 1)
])                            # (131, 5): cx, cy, r_droplet, r_bead (px), status


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
def weighted_quantile(x, w, q):
    i = np.argsort(x)
    x, w = x[i], w[i]
    cw = np.cumsum(w)
    cw /= cw[-1]
    return np.interp(q / 100.0, cw, x)


def weighted_mean(x, w):
    return np.sum(x * w) / np.sum(w)


def weighted_cv(x, w):
    m = weighted_mean(x, w)
    return 100.0 * np.sqrt(np.sum(w * (x - m) ** 2) / np.sum(w)) / m


# ----------------------------------------------------------------------------
# 1. Scale validation: 140 hand-marked 3 um beads
# ----------------------------------------------------------------------------
bead_d_um = 2.0 * _SCALE_BEADS[:, 2] * UMPX
print("\n[1] Scale validation (n = %d beads, nominal 3 um)" % len(bead_d_um))
print("    image scale          : %.5f um/px (50 um bar = 333 px)" % UMPX)
print("    bead diameter median : %.2f um (mean %.2f, CV %.0f%%)"
      % (np.median(bead_d_um), bead_d_um.mean(),
         100 * bead_d_um.std(ddof=1) / bead_d_um.mean()))

# ----------------------------------------------------------------------------
# 2. Droplet population: bead-bearing droplets and 1/V bias correction
# ----------------------------------------------------------------------------
status = _MEAS[:, 4].astype(bool)
d_all = 2.0 * _MEAS[:, 2] * UMPX
v_all = np.pi / 6.0 * d_all ** 3          # fL (1 um^3 = 1 fL)
d = d_all[status]
v = v_all[status]
w = 1.0 / v                               # P(bead | droplet) proportional to V

print("\n[2] Segmentation: %d anchors -> %d bead-bearing droplets, "
      "%d bare beads excluded (%.0f%%)"
      % (len(_MEAS), status.sum(), (~status).sum(), 100 * (~status).sum() / len(_MEAS)))
print("    bead-bearing droplets: diameter median %.2f um, CV %.0f%%; "
      "volume median %.0f fL, mean %.0f fL"
      % (np.median(d), 100 * d.std(ddof=1) / d.mean(), np.median(v), v.mean()))

d5c = weighted_quantile(d, w, 5)
d50c = weighted_quantile(d, w, 50)
d95c = weighted_quantile(d, w, 95)
v5c = weighted_quantile(v, w, 5)
v50c = weighted_quantile(v, w, 50)
v95c = weighted_quantile(v, w, 95)
vmeanc = weighted_mean(v, w)
cvd = weighted_cv(d, w)
cvv = weighted_cv(v, w)
n_droplets = VOLUME_ASSAY_UL * 1e-6 / (vmeanc * 1e-15)

print("    full population (1/V-weighted): diameter median %.2f um, CV %.0f%%; "
      "volume median %.0f fL, mean %.0f fL, CV %.0f%%"
      % (d50c, cvd, v50c, vmeanc, cvv))
print("    volume 5-95%% interval: %.0f-%.0f fL (%.1fx); "
      "droplets per 50 uL: ~%.1e" % (v5c, v95c, v95c / v5c, n_droplets))

# ----------------------------------------------------------------------------
# 3. Concentration bias if droplet volume were used as the counting partition
# ----------------------------------------------------------------------------
bias_at = lambda vx: (v50c / vx - 1.0) * 100.0
print("\n[3] Volume-partition bias: %+.0f%% at the 5th percentile, "
      "%+.0f%% at the 95th (relative to the median volume)"
      % (bias_at(v5c), bias_at(v95c)))

# ----------------------------------------------------------------------------
# 4. Figure: Extended Data Fig. 6
# ----------------------------------------------------------------------------
apply_style()
fig, (axa, axb) = plt.subplots(1, 2, figsize=(7.3, 3.6))

# panel a: diameter distribution
bins = np.arange(2.5, 10.75, 0.5)
axa.hist(d, bins=bins, color=OI['sky'], edgecolor='black', linewidth=0.8, zorder=2)
axa.set_xlabel('Droplet diameter (µm)')
axa.set_ylabel('Count')
axa.set_xlim(2.5, 10.5)
axa.set_ylim(0, 21)
axa.yaxis.set_major_locator(MaxNLocator(integer=True))
axa.axvline(NOMINAL_BEAD_UM, color='0.35', ls=':', lw=1.5, zorder=3)
axa.text(3.1, 19.6, 'bead Ø 3 µm', fontsize=9.5, color='0.35', rotation=90, va='top')
axa.axvline(np.median(d), color=OI['vermillion'], ls='--', lw=1.5, zorder=3)
txt = (f'n = {len(d)} bead-bearing\n'
       f'median {np.median(d):.1f} µm ({np.median(v):.0f} fL), '
       f'CV {100*d.std(ddof=1)/d.mean():.0f}%\n'
       f'all droplets (bias-corrected):\n'
       f'median {d50c:.1f} µm ({v50c:.0f} fL), CV {cvd:.0f}%')
axa.text(0.97, 0.95, txt, transform=axa.transAxes, ha='right', va='top',
         fontsize=9, bbox=dict(boxstyle='round,pad=0.35', fc='white',
                               ec='0.6', lw=0.8))
d2v = lambda x: np.pi / 6 * np.maximum(x, 0) ** 3
v2d = lambda x: (6 * np.maximum(x, 0) / np.pi) ** (1 / 3)
sec = axa.secondary_xaxis('top', functions=(d2v, v2d))
sec.set_xticks([10, 25, 50, 100, 250, 500])
sec.set_xticklabels(['10', '25', '50', '100', '250', '500'])
sec.set_xlabel('Droplet volume (fL)')
sec.tick_params(axis='x', direction='out', width=1.0, labelsize=9.5)

# panel b: volume-partition bias
vv = np.logspace(np.log10(10), np.log10(600), 400)
bias = bias_at(vv)
axb.plot(vv, bias, color=OI['vermillion'], lw=2.0, zorder=4)
axb.axhline(0, color='0.35', ls='--', lw=1.2, zorder=2)
kde = stats.gaussian_kde(np.log10(v), weights=w, bw_method=0.35)
dens = kde(np.log10(vv))
dens = dens / dens.max()
ax2 = axb.twinx()
ax2.fill_between(vv, 0, dens, color=OI['sky'], alpha=0.35, lw=0, zorder=1)
ax2.set_ylim(0, 3.2)
ax2.set_yticks([])
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
for vx, lab, dx, dy in ((v5c, '%+.0f%%' % bias_at(v5c), -34, 10),
                        (v95c, '%+.0f%%' % bias_at(v95c), 4, 14)):
    axb.plot([vx], [bias_at(vx)], 'o', ms=6, mfc='white',
             mec=OI['vermillion'], mew=1.6, zorder=5)
    axb.axvline(vx, color='0.55', ls=':', lw=1.2, zorder=2)
    axb.annotate(lab, (vx, bias_at(vx)), textcoords='offset points',
                 xytext=(dx, dy), ha='center', fontsize=10.5,
                 color=OI['vermillion'], zorder=6)
axb.text(0.30, 0.16, 'droplet volume\ndistribution', transform=axb.transAxes,
         fontsize=9.5, color=OI['blue'], ha='center', va='top')
axb.text(0.96, 0.94, '5\u201395%% interval:\n%.0f\u2013%.0f fL' % (v5c, v95c),
         transform=axb.transAxes, fontsize=9.5, ha='right', va='top', color='0.25')
axb.set_xscale('log')
axb.set_xlim(10, 600)
axb.set_xticks([10, 20, 50, 100, 200, 500])
axb.set_xticklabels(['10', '20', '50', '100', '200', '500'])
axb.set_xlabel('Droplet volume (fL)')
axb.set_ylabel('Concentration bias (%)')
axb.set_ylim(-100, 260)
axb.set_yticks([-100, -50, 0, 50, 100, 150, 200, 250])

panel_label(axa, 'a')
panel_label(axb, 'b')
fig.tight_layout()
save_fig(fig, 'ED_Fig_6')
print("\n[4] saved ED_Fig_6.png / .svg")

# ----------------------------------------------------------------------------
# 5. Excel workbook
# ----------------------------------------------------------------------------
readme = pd.DataFrame({'ED Fig. 6 droplet morphometry': [
    'Droplet-based digital immunoassay (manuscript Fig. 3c-e; SI S12).',
    'Bright-field image of the bead monolayer, same experiments as the assay.',
    'Bead-anchored ray-casting: 72 rays per bead, 55% half-height of the inner rim flank.',
    'Anchors with r_droplet < 1.25x r_bead are bare beads and are excluded',
    '(two borderline objects excluded on visual inspection; see status column).',
    'Scale: 50 um bar = 333 px (%.5f um/px), verified by 140 hand-marked 3 um beads.' % UMPX,
    'Bias correction: P(bead | droplet) proportional to volume, weight = 1/V.',
    'Sheets: Summary_stats, Droplet_measurements, Scale_beads, Bias_curve.']})

summary = pd.DataFrame([
    ('image scale', UMPX, 'um/px', '50 um drawn bar = 333 px'),
    ('scale beads n', len(bead_d_um), '', 'hand-marked 3 um beads'),
    ('bead diameter median', np.median(bead_d_um), 'um', 'nominal 3.0 um'),
    ('bead diameter mean', bead_d_um.mean(), 'um', ''),
    ('bead diameter CV', 100 * bead_d_um.std(ddof=1) / bead_d_um.mean(), '%', ''),
    ('anchors measured', len(_MEAS), '', 'hand-marked bead anchors'),
    ('bare beads excluded', int((~status).sum()), '', 'r_droplet < 1.25x r_bead'),
    ('bead-bearing droplets n', int(status.sum()), '', ''),
    ('droplet diameter median (bead-bearing)', np.median(d), 'um', ''),
    ('droplet diameter CV (bead-bearing)', 100 * d.std(ddof=1) / d.mean(), '%', ''),
    ('droplet volume median (bead-bearing)', np.median(v), 'fL', ''),
    ('droplet volume mean (bead-bearing)', v.mean(), 'fL', ''),
    ('droplet diameter median (1/V-weighted)', d50c, 'um', 'full droplet population'),
    ('droplet diameter CV (1/V-weighted)', cvd, '%', ''),
    ('droplet volume median (1/V-weighted)', v50c, 'fL', ''),
    ('droplet volume mean (1/V-weighted)', vmeanc, 'fL', ''),
    ('droplet volume CV (1/V-weighted)', cvv, '%', ''),
    ('droplet volume 5th percentile', v5c, 'fL', ''),
    ('droplet volume 95th percentile', v95c, 'fL', ''),
    ('droplets per 50 uL emulsion', n_droplets, '', 'from 1/V-weighted mean volume'),
    ('volume-partition bias at 5th pct', bias_at(v5c), '%', 'c proportional to 1/V'),
    ('volume-partition bias at 95th pct', bias_at(v95c), '%', ''),
], columns=['quantity', 'value', 'unit', 'note'])

meas = pd.DataFrame({
    'anchor_id': np.arange(1, len(_MEAS) + 1),
    'cx_px': _MEAS[:, 0],
    'cy_px': _MEAS[:, 1],
    'r_bead_px': _MEAS[:, 3],
    'r_droplet_px': _MEAS[:, 2],
    'r_ratio': _MEAS[:, 2] / _MEAS[:, 3],
    'status': np.where(status, 'bead-bearing droplet', 'bare bead (excluded)'),
    'diameter_um': d_all,
    'volume_fL': v_all,
    'weight_1_over_V': np.where(status, 1.0 / v_all, np.nan),
})

beads_df = pd.DataFrame({
    'bead_id': np.arange(1, len(_SCALE_BEADS) + 1),
    'cx_px': _SCALE_BEADS[:, 0],
    'cy_px': _SCALE_BEADS[:, 1],
    'r_px': _SCALE_BEADS[:, 2],
    'diameter_um': bead_d_um,
})

bias_df = pd.DataFrame({'volume_fL': vv, 'concentration_bias_pct': bias})
kde_df = pd.DataFrame({'volume_fL': vv, 'distribution_density_norm': dens})

with pd.ExcelWriter('ED_Fig_6_droplet_morphometry.xlsx',
                    engine='openpyxl') as xw:
    readme.to_excel(xw, sheet_name='README', index=False)
    summary.to_excel(xw, sheet_name='Summary_stats', index=False)
    meas.to_excel(xw, sheet_name='Droplet_measurements', index=False)
    beads_df.to_excel(xw, sheet_name='Scale_beads', index=False)
    bias_df.to_excel(xw, sheet_name='Bias_curve', index=False)
    kde_df.to_excel(xw, sheet_name='Volume_distribution', index=False)
print("[5] saved ED_Fig_6_droplet_morphometry.xlsx")
print("=" * 80)
