# Target-Capture Scaling (TCS): reproducibility capsule

Companion code for "Scale degeneracy: the physical basis of calibration and
absolute quantification in molecular binding" (preprint: ChemRxiv,
DOI: 10.26434/chemrxiv-2024-19rj6).

The capsule is fully self-contained: every experimental dataset used in the
manuscript (digitized Yalow and Berson 1960 radioimmunoassay values, Simoa
single-molecule array bead counts, digital PCR partition counts, and the
immunofluorescence scale-degeneracy experiment) is embedded directly in the
scripts. No external data download is required.

## Environment

Python 3.11 with the packages listed in `requirements.txt`:

```
pip install -r requirements.txt
```

## Quick start

Run any single script directly, for example:

```
python Main_Fig_1_abc.py
```

Each script writes its outputs (figures as PDF/PNG, tables as XLSX/CSV/TXT)
into the working directory. To reproduce everything, execute the `run` script,
which runs all 26 analysis scripts in sequence and collects the outputs
(`tcs_figstyle.py` is a shared style module imported by the figure scripts,
not an analysis script):

```
bash run
```

## Script to result mapping

| Script | Reproduces |
| --- | --- |
| Main_Fig_1_abc.py | Fig. 1a-c (TCS framework figure) |
| Main_Fig_1_and_ED_Fig_1_ab.py | Fig. 1 and Supplementary Figure ED1a-b (in silico master-equation validation) |
| SD_1_and_Fig_3a_scale_degeneracy.py | Supplementary Data 1 and Fig. 3a (immunofluorescence scale-degeneracy experiment, 72 measurements) |
| SI_Table_4_and_Fig_3b_simoa.py | Fig. 3b (IL-6 Simoa reanalysis), SI Table 4, Supplementary Data |
| Main_Fig_3_cde_droplet.py | Fig. 3c-e (home-built droplet digital immunoassay: chip, equilibrium response, clinical correlation) and Supplementary Data 6 (SI_Table_6.xlsx: raw counts, fit parameters, LoB/LoD/LoQ, clinical values, bootstrap samples) |
| SD_2_ria_yalow1960.py | Supplementary Data 2 (Yalow 1960 RIA reanalysis) and Supplementary Figure ED4 |
| SD_3_and_ED_Fig_5_sbg_simoa.py | Supplementary Data 3 and Supplementary Figure ED5 (2010 Simoa benchmark) |
| SD_5_and_ED_Fig_6_dpcr.py | Supplementary Data 5 and Supplementary Figure ED6 (dPCR analysis) |
| ED_Fig_2a_4pl_5pl_comparison.py | Supplementary Figure ED2a (4PL/5PL versus exact TCS) |
| ED_Fig_3_droplet_morphometry.py | Supplementary Figure ED3 and ED Fig. 3 workbook (droplet morphometry, SI S12) |
| SI_S3_4pl_5pl_verification.py | SI S3: symbolic derivation check of 4PL/5PL from TCS |
| SI_S5_pico_verification.py | SI S5: PICO platform verification |
| SI_Fig_S1c_1_error_scan.py | SI Fig. S1c.1 (error scan) |
| Fig_2_cv_estimator.py | Fig. 2 (exact CV of the count estimator; formerly SI Fig. S2b.1, moved to the main text) |
| SI_Fig_S2b_2_excess_variance.py | SI Fig. S2b.2 (excess variance) |
| SI_S2b_6_7_optimality_lod_tails.py | SI S2b.6-S2b.7 (optimality and LoD tails) |
| SI_S2c_5_2_profile_likelihood.py | SI S2c.5.2 (profile likelihood; Bayesian, long runtime) |
| SI_S2c_5_5_fisher_verification.py | SI S2c.5.5 (Fisher information verification) |
| SI_S2d_5_analog_fisher_scan.py | SI S2d.5 (analog Fisher scan) |
| SI_S2e_10_temperature_identifiability.py | SI S2e.10 (multi-temperature identifiability) |
| SI_S2f_b_uncertainty_verification.py | SI S2f (uncertainty formulas verification) |
| SI_Table_S2b_2_icc_variance.py | SI Table S2b.2 (ICC variance) |
| SI_Fig_S10_1_scatchard.py | SI Fig. S10.1 (Scatchard analysis) |
| SI_Fig_S10_2_and_Table_S10_1.py | SI Fig. S10.2 and Table S10.1 |
| SI_Table_S9_1_and_Fig_S9_1.py | SI Table S9.1 and Fig. S9.1 |
| SI_Table_S9_2_and_Fig_S9_2.py | SI Table S9.2 and Fig. S9.2 |
| SI_S9_deep_purification_biosorption.py | SI S9 deep-purification criterion: κ distribution across 361 biosorption systems from Febrianto et al. (2009) Table 2 (SI_S9_biosorption_kappa_distribution.csv) |

## Runtime notes

Most scripts finish in seconds to a few minutes. The Bayesian inference
scripts (profile likelihood, dPCR MCMC, IL-6 Simoa fitting) use emcee,
dynesty, and arviz and can take considerably longer depending on hardware.
Every script can also be run standalone for inspection of a single result.

## Code availability

Repository: https://github.com/lihuanlz/TCS_v4.0
Archive: https://zenodo.org/records/22970538
