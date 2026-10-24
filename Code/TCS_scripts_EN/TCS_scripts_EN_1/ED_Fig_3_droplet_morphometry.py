# -*- coding: utf-8 -*-
"""
Supplementary Figure ED3 (droplet morphometry) and accompanying Excel workbook
for the droplet-based digital immunoassay platform (manuscript Fig. 3c-e,
Supplementary Information S12).

Provenance of the data
----------------------
Bright-field image of the bead monolayer acquired in the observation chip in
the same experiments as the assay (Fig. 3c). On the full 1000 x 1000 px image
all clearly bounded droplets were hand-circled (n = 387) and the magnetic
beads were hand-circled separately (n = 342); circle radii are the
measurements. A droplet is bead-bearing when a bead circle lies fully inside
its circle.

The observation field is enriched for bead-bearing droplets (magnetic
enrichment), a selection whose probability is proportional to droplet volume
V; the full emulsion population is therefore recovered by weighting each
droplet by 1/V. Droplets under the drawn scale bar were not marked.

Image scale: the drawn 50 um scale bar measures 333 px
(0.15015 um/px), independently verified against the nominal 3 um bead
diameter (342 hand-marked beads, median 2.94 um).

Inputs (raw data, same directory)
---------------------------------
brightfield_droplet_annotations.csv : id, cx, cy, r (px) - hand-circled droplets
brightfield_bead_annotations.csv    : id, cx, cy, r (px) - hand-circled beads

Outputs
-------
ED_Fig_3.png / .svg : (a) droplet diameter distribution,
                      (b) concentration bias for a volume partition
ED_Fig_3_droplet_morphometry.xlsx : README, Summary_stats, Droplet_measurements,
                                    Scale_beads, Bias_curve, Volume_distribution
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy import stats
from tcs_figstyle import apply_style, OI, save_fig, panel_label

print("=" * 80)
print("ED Fig. 3 - droplet morphometry of the droplet-based digital immunoassay")
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
# Raw data (hand-circled ROIs, full 1000 x 1000 px bright-field image)
# ----------------------------------------------------------------------------
beads = pd.read_csv('brightfield_bead_annotations.csv')
drop = pd.read_csv('brightfield_droplet_annotations.csv')

# bead-bearing status: a bead circle fully inside the droplet circle
n_beads = np.zeros(len(drop), dtype=int)
for j, r in drop.iterrows():
    dist = np.hypot(beads.cx - r.cx, beads.cy - r.cy)
    n_beads[j] = int(((dist + beads.r) < r.r).sum())
bearing = n_beads > 0

bare = np.ones(len(beads), dtype=bool)
for i, r in beads.iterrows():
    dist = np.hypot(drop.cx - r.cx, drop.cy - r.cy)
    bare[i] = not bool(((dist + r.r) < drop.r).any())

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
# 1. Scale validation: 331 hand-marked 3 um beads
# ----------------------------------------------------------------------------
bead_d_um = 2.0 * beads.r.values * UMPX
print("\n[1] Scale validation (n = %d beads, nominal 3 um)" % len(bead_d_um))
print("    image scale          : %.5f um/px (50 um bar = 333 px)" % UMPX)
print("    bead diameter median : %.2f um (mean %.2f, CV %.0f%%)"
      % (np.median(bead_d_um), bead_d_um.mean(),
         100 * bead_d_um.std(ddof=1) / bead_d_um.mean()))
print("    bare beads (not inside any droplet): %d (%.0f%%)"
      % (bare.sum(), 100 * bare.mean()))

# ----------------------------------------------------------------------------
# 2. Droplet population: enriched field sample and 1/V bias correction
# ----------------------------------------------------------------------------
d_all = 2.0 * drop.r.values * UMPX
v_all = np.pi / 6.0 * d_all ** 3          # fL (1 um^3 = 1 fL)
d, v = d_all, v_all
w = 1.0 / v                               # enrichment probability proportional to V

print("\n[2] Field sample: %d droplets circled; %d bead-bearing (%.0f%%), "
      "%d with 2+ beads"
      % (len(drop), bearing.sum(), 100 * bearing.mean(),
         int((n_beads >= 2).sum())))
print("    field droplets: diameter median %.2f um, CV %.0f%%; "
      "volume median %.0f fL, mean %.0f fL"
      % (np.median(d), 100 * d.std(ddof=1) / d.mean(), np.median(v), v.mean()))
db, vb = d_all[bearing], v_all[bearing]
print("    bead-bearing subset: n = %d, diameter median %.2f um; "
      "volume median %.0f fL" % (bearing.sum(), np.median(db), np.median(vb)))

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

print("    full emulsion (1/V-weighted): diameter median %.2f um, CV %.0f%%; "
      "volume median %.0f fL, mean %.0f fL, CV %.0f%%"
      % (d50c, cvd, v50c, vmeanc, cvv))
print("    volume 5-95%% interval: %.0f-%.0f fL (%.1fx); "
      "droplets per 50 uL: ~%.1e" % (v5c, v95c, v95c / v5c, n_droplets))

# ----------------------------------------------------------------------------
# 3. Concentration bias if droplet volume were used as the counting partition
# ----------------------------------------------------------------------------
bias_at = lambda vx: (v50c / vx - 1.0) * 100.0
print("\n[3] Volume-partition bias: %+.0f%% at the 5th percentile, "
      "%+.0f%% at the 95th (relative to the weighted median volume)"
      % (bias_at(v5c), bias_at(v95c)))

# ----------------------------------------------------------------------------
# 4. Figure: Extended Data Fig. 3
# ----------------------------------------------------------------------------
apply_style()
fig, (axa, axb) = plt.subplots(1, 2, figsize=(7.3, 3.6))

# panel a: diameter distribution
bins = np.arange(2.5, 14.0, 0.5)
counts, _, _ = axa.hist(d, bins=bins, color=OI['sky'], edgecolor='black',
                        linewidth=0.8, zorder=2)
axa.set_xlabel('Droplet diameter (µm)')
axa.set_ylabel('Count')
axa.set_xlim(2.5, 13.5)
axa.set_ylim(0, counts.max() * 1.30)
axa.yaxis.set_major_locator(MaxNLocator(integer=True))
axa.axvline(NOMINAL_BEAD_UM, color='0.35', ls=':', lw=1.5, zorder=3)
axa.text(3.15, counts.max() * 1.22, 'bead Ø 3 µm', fontsize=9.5,
         color='0.35', rotation=90, va='top')
axa.axvline(np.median(d), color=OI['vermillion'], ls='--', lw=1.5, zorder=3)
txt = (f'n = {len(d)} circled droplets\n'
       f'median {np.median(d):.1f} µm ({np.median(v):.0f} fL), '
       f'CV {100*d.std(ddof=1)/d.mean():.0f}%\n'
       f'full emulsion (1/V-corrected):\n'
       f'median {d50c:.1f} µm ({v50c:.0f} fL), CV {cvd:.0f}%')
axa.text(0.97, 0.95, txt, transform=axa.transAxes, ha='right', va='top',
         fontsize=9, bbox=dict(boxstyle='round,pad=0.35', fc='white',
                               ec='0.6', lw=0.8))
d2v = lambda x: np.pi / 6 * np.maximum(x, 0) ** 3
v2d = lambda x: (6 * np.maximum(x, 0) / np.pi) ** (1 / 3)
sec = axa.secondary_xaxis('top', functions=(d2v, v2d))
sec.set_xticks([10, 25, 50, 100, 250, 500, 1000])
sec.set_xticklabels(['10', '25', '50', '100', '250', '500', '1000'])
sec.set_xlabel('Droplet volume (fL)')
sec.tick_params(axis='x', direction='out', width=1.0, labelsize=9.5)

# panel b: volume-partition bias
vv = np.logspace(np.log10(15), np.log10(800), 400)
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
for vx, lab, dx, dy in ((v5c, '%+.0f%%' % bias_at(v5c), 38, 8),
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
axb.set_xlim(15, 800)
axb.set_xticks([20, 50, 100, 200, 500])
axb.set_xticklabels(['20', '50', '100', '200', '500'])
axb.set_xlabel('Droplet volume (fL)')
axb.set_ylabel('Concentration bias (%)')
axb.set_ylim(-100, 300)
axb.set_yticks([-100, -50, 0, 50, 100, 150, 200, 250, 300])

panel_label(axa, 'a')
panel_label(axb, 'b')
fig.tight_layout()
save_fig(fig, 'ED_Fig_3')
print("\n[4] saved ED_Fig_3.png / .svg")

# ----------------------------------------------------------------------------
# 5. Excel workbook
# ----------------------------------------------------------------------------
readme = pd.DataFrame({'ED Fig. 3 droplet morphometry': [
    'Droplet-based digital immunoassay (manuscript Fig. 3c-e; SI S12).',
    'Bright-field image of the bead monolayer, same experiments as the assay.',
    'All clearly bounded droplets (n = %d) and all beads (n = %d) hand-circled'
    ' on the full 1000 x 1000 px image; circle radius is the measurement.'
    % (len(drop), len(beads)),
    'A droplet is bead-bearing when a bead circle lies fully inside it.',
    'The field is enriched for bead-bearing droplets (P proportional to V);',
    'the full emulsion population is recovered by weight = 1/V.',
    'Droplets under the drawn 50 um scale bar were not marked.',
    'Scale: 50 um bar = 333 px (%.5f um/px), verified by %d hand-marked 3 um beads.'
    % (UMPX, len(beads)),
    'Sheets: Summary_stats, Droplet_measurements, Scale_beads, Bias_curve, '
    'Volume_distribution.']})

summary = pd.DataFrame([
    ('image scale', UMPX, 'um/px', '50 um drawn bar = 333 px'),
    ('scale beads n', len(bead_d_um), '', 'hand-marked 3 um beads'),
    ('bead diameter median', np.median(bead_d_um), 'um', 'nominal 3.0 um'),
    ('bead diameter mean', bead_d_um.mean(), 'um', ''),
    ('bead diameter CV', 100 * bead_d_um.std(ddof=1) / bead_d_um.mean(), '%', ''),
    ('bare beads (not in any droplet)', int(bare.sum()), '', ''),
    ('droplets circled n', len(drop), '', 'all clearly bounded droplets in the field'),
    ('bead-bearing droplets n', int(bearing.sum()), '', 'bead circle fully inside droplet'),
    ('droplet diameter median (field)', np.median(d), 'um', 'enriched field sample'),
    ('droplet diameter CV (field)', 100 * d.std(ddof=1) / d.mean(), '%', ''),
    ('droplet volume median (field)', np.median(v), 'fL', ''),
    ('droplet volume mean (field)', v.mean(), 'fL', ''),
    ('bead-bearing diameter median', np.median(db), 'um', ''),
    ('bead-bearing volume median', np.median(vb), 'fL', ''),
    ('droplet diameter median (1/V-weighted)', d50c, 'um', 'full emulsion population'),
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
    'droplet_id': drop.id.values,
    'cx_px': drop.cx.values,
    'cy_px': drop.cy.values,
    'r_px': drop.r.values,
    'diameter_um': d_all,
    'volume_fL': v_all,
    'n_beads_inside': n_beads,
    'status': np.where(bearing, 'bead-bearing', 'empty'),
    'weight_1_over_V': 1.0 / v_all,
})

beads_df = pd.DataFrame({
    'bead_id': beads.id.values,
    'cx_px': beads.cx.values,
    'cy_px': beads.cy.values,
    'r_px': beads.r.values,
    'diameter_um': bead_d_um,
    'status': np.where(bare, 'bare (not in a droplet)', 'inside a droplet'),
})

bias_df = pd.DataFrame({'volume_fL': vv, 'concentration_bias_pct': bias})
kde_df = pd.DataFrame({'volume_fL': vv, 'distribution_density_norm': dens})

with pd.ExcelWriter('ED_Fig_3_droplet_morphometry.xlsx',
                    engine='openpyxl') as xw:
    readme.to_excel(xw, sheet_name='README', index=False)
    summary.to_excel(xw, sheet_name='Summary_stats', index=False)
    meas.to_excel(xw, sheet_name='Droplet_measurements', index=False)
    beads_df.to_excel(xw, sheet_name='Scale_beads', index=False)
    bias_df.to_excel(xw, sheet_name='Bias_curve', index=False)
    kde_df.to_excel(xw, sheet_name='Volume_distribution', index=False)
print("[5] saved ED_Fig_3_droplet_morphometry.xlsx")
print("=" * 80)
