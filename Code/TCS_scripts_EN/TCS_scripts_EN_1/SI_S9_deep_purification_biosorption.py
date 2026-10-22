# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 14:56:00 2026

@author: lihua
"""

"""
S9 Biosorption κ Distribution Calculator (Final)
=================================================
Strictly follows Table 2 of Febrianto et al. (2009) J. Hazard. Mater. 162, 616-645.
Excludes: negative KL, zero KL, missing qmax, negative qmax.

TCS definitions (S9.1):
  K = 1/K_L          (dissociation constant)
  C_site = qmax*m/V  (total binding sites, m/V = 1 g/L)
  κ = K / C_site = 1/(K_L * qmax * m/V)
  R = κ/(1+κ)        (exact residual fraction)
"""

import csv

# ============================================================
# Table 2 data: (biosorbent, solute, KL [L/mg], qmax [mg/g])
# All entries programmatically verified against the original table
# (361/361 rows match the PDF; see verification notes in README).
# Excluded: negative KL (saltbush Pb/Zn, unmodified maize husk Cd/Pb/Zn,
#           Salvinia cucullata Cr), zero KL (marine algae DP95Ca pH 1, Pb and Cu).
# ============================================================
data = [
    # Aspergillus niger Cu(II)
    ("Aspergillus niger", "Cu(II)", 0.0120, 33.11),
    ("Aspergillus niger", "Cu(II)", 0.0150, 33.23),
    ("Aspergillus niger", "Cu(II)", 0.0260, 33.54),
    # Aspergillus niger Pb(II)
    ("Aspergillus niger", "Pb(II)", 0.0210, 34.69),
    ("Aspergillus niger", "Pb(II)", 0.0220, 34.92),
    ("Aspergillus niger", "Pb(II)", 0.0390, 34.92),
    # Cladophora fascicularis Cu(II)
    ("Cladophora fascicularis", "Cu(II)", 0.1406, 47.019),
    ("Cladophora fascicularis", "Cu(II)", 0.1607, 70.53),
    ("Cladophora fascicularis", "Cu(II)", 0.1649, 92.77),
    ("Cladophora fascicularis", "Cu(II)", 0.1882, 110.56),
    # Cassia fistula Ni(II)
    ("Cassia fistula leaves", "Ni(II)", 0.0126, 163.93),
    ("Cassia fistula stem bark", "Ni(II)", 0.0020, 172.41),
    ("Cassia fistula pods bark", "Ni(II)", 0.0345, 196.07),
    # Rhizopus arrhizus
    ("Rhizopus arrhizus", "Cr(VI)", 0.3310, 23.92),
    # Crab shell Cu(II)
    ("Crab shell", "Cu(II)", 0.0022, 163.90),
    ("Crab shell", "Cu(II)", 0.0023, 188.70),
    ("Crab shell", "Cu(II)", 0.0026, 204.10),
    ("Crab shell", "Cu(II)", 0.0039, 208.30),
    ("Crab shell", "Cu(II)", 0.0054, 222.20),
    ("Crab shell", "Cu(II)", 0.0055, 243.90),
    # Crab shell Co(II)
    ("Crab shell", "Co(II)", 0.0015, 212.80),
    ("Crab shell", "Co(II)", 0.0016, 235.50),
    ("Crab shell", "Co(II)", 0.0017, 270.30),
    ("Crab shell", "Co(II)", 0.0018, 285.70),
    ("Crab shell", "Co(II)", 0.0021, 303.00),
    ("Crab shell", "Co(II)", 0.0025, 322.60),
    # Dried activated sludge
    ("Dried activated sludge", "Cd(II)", 0.0210, 84.30),
    ("Dried activated sludge", "Pb(II)", 0.0320, 131.60),
    # Tamarindus indica seeds
    ("Tamarindus indica seeds", "Cr(VI)", 16023.07, 0.098),
    ("Tamarindus indica seeds", "Cr(VI)", 3191.93, 0.055),
    ("Tamarindus indica seeds", "Cr(VI)", 72.43, 0.08),
    ("Tamarindus indica seeds", "Cr(VI)", 1.3600, 0.023),
    # Dried activates sludge
    ("Dried activated sludge", "Cr(II)", 0.0063, 294.10),
    ("Dried activated sludge", "Cr(II)", 0.0055, 95.20),
    ("Dried activated sludge", "Ni(II)", 0.0079, 106.40),
    ("Dried activated sludge", "Ni(II)", 0.0048, 238.10),
    # Palm Tree Leaves
    ("Palm Tree Leaves", "Zn(II)", 0.0560, 14.60),
    # Tea waste
    ("Tea waste", "Cu(II)", 0.0076, 48.00),
    ("Tea waste", "Pb(II)", 0.0494, 65.00),
    # Carica papaya
    ("Carica papaya", "Hg(II)", 0.0040, 155.63),
    # Wool etc Cr(VI) - 7 materials
    ("Wool", "Cr(VI)", 7.15e-3, 41.15),
    ("Olive cake", "Cr(VI)", 4.70e-3, 33.44),
    ("Sawdust", "Cr(VI)", 9.15e-3, 15.28),
    ("Pine needles", "Cr(VI)", 5.44e-3, 21.50),
    ("Almond", "Cr(VI)", 5.46e-3, 10.62),
    ("Coal", "Cr(VI)", 11.50e-3, 6.78),
    ("Cactus", "Cr(VI)", 6.13e-3, 7.08),
    # Coconut copra meal
    ("Coconut copra meal", "Cd(II)", 0.1760, 4.92),
    ("Coconut copra meal", "Cd(II)", 0.1410, 4.68),
    ("Coconut copra meal", "Cd(II)", 0.1190, 2.66),
    ("Coconut copra meal", "Cd(II)", 0.0990, 2.01),
    # Rose waste biomass
    ("Rose waste biomass", "Pb(II)", 0.0280, 151.51),
    ("Rose waste biomass", "Co(II)", 0.0372, 27.62),
    # Cotton cellulose B(III)
    ("Cotton cellulose", "B(III)", 0.3040, 6.88),
    ("Cotton cellulose", "B(III)", 0.7370, 13.06),
    ("Cotton cellulose", "B(III)", 1.3580, 15.41),
    ("Cotton cellulose", "B(III)", 0.8460, 17.79),
    ("Cotton cellulose", "B(III)", 0.8660, 41.49),
    # Spent grain
    ("Spent grain", "Cu(II)", 0.0800, 10.47),
    # Tea factory waste
    ("Tea factory waste", "Ni(II)", 0.0880, 15.26),
    ("Tea factory waste", "Ni(II)", 0.0930, 17.73),
    ("Tea factory waste", "Ni(II)", 0.0980, 18.42),
    # Waste beer yeast
    ("Waste beer yeast", "Pb(II)", 0.0883, 55.71),
    # Green coconut shell powder
    ("Green coconut shell powder", "Cd(II)", 0.0190, 285.70),
    # Crude/HCl/Oxalic tamarind fruit shells
    ("Crude tamarind fruit shells", "Cr(VI)", 0.0500, 74.62),
    ("HCl treated tamarind shells", "Cr(VI)", 0.0700, 140.84),
    ("Oxalic acid treated tamarind shells", "Cr(VI)", 0.0350, 151.51),
    # Leaves of saltbush Cu(II) only (Pb/Zn have negative KL)
    ("Leaves of saltbush", "Cu(II)", 9.0900, 590.92),
    # Rice bran Zn(II)
    ("Rice bran", "Zn(II)", 0.0110, 14.17),
    ("Rice bran", "Zn(II)", 0.0150, 14.84),
    ("Rice bran", "Zn(II)", 0.0190, 18.31),
    # Waste weed Salvinia (only valid row, first has negative KL)
    ("Waste weed Salvinia", "Cr(VI)", 0.0130, 232.60),
    # Saccharomyces cerevisiae
    ("Saccharomyces cerevisiae", "Cr(VI)", 0.1259, 384.61),
    # Spirulina platensis
    ("Spirulina platensis fresh", "Cr(VI)", 0.0270, 188.68),
    ("Spirulina platensis spent", "Cr(VI)", 0.0500, 212.76),
    ("Chlorella vulgaris", "Cr(VI)", 0.0360, 163.93),
    # Caulerpa lentillifera
    ("Caulerpa lentillifera", "Cu(II)", 0.0760, 5.57),
    ("Caulerpa lentillifera", "Cu(II)", 0.0618, 42.37),
    ("Caulerpa lentillifera", "Cu(II)", 0.0290, 2.64),
    ("Caulerpa lentillifera", "Cd(II)", 0.0742, 4.69),
    ("Caulerpa lentillifera", "Cd(II)", 0.0229, 4.34),
    ("Caulerpa lentillifera", "Cd(II)", 0.0106, 2.95),
    ("Caulerpa lentillifera", "Pb(II)", 0.0711, 28.99),
    ("Caulerpa lentillifera", "Pb(II)", 0.0324, 24.19),
    ("Caulerpa lentillifera", "Pb(II)", 0.0144, 15.94),
    ("Caulerpa lentillifera", "Zn(II)", 0.0671, 2.66),
    ("Caulerpa lentillifera", "Zn(II)", 0.0500, 1.65),
    ("Caulerpa lentillifera", "Zn(II)", 0.0201, 1.37),
    # Calcium alginate etc
    ("Calcium alginate", "Ni(II)", 0.0036, 310.40),
    ("Chitosan coated calcium alginate", "Ni(II)", 0.0240, 222.20),
    ("Chitosan coated silica", "Ni(II)", 0.0026, 254.30),
    # Chitosan beads
    ("Chitosan beads", "Cr(III)", 1.320e-4, 30.03),
    ("Chitosan beads", "Cr(VI)", 3.421e-3, 76.92),
    # Cladophora fascicularis (non-living)
    ("Cladophora fascicularis NL", "Cu(II)", 0.0652, 102.55),
    ("Cladophora fascicularis NL", "Cu(II)", 0.1241, 106.81),
    ("Cladophora fascicularis NL", "Cu(II)", 0.1750, 112.97),
    ("Cladophora fascicularis NL", "Pb(II)", 0.0396, 200.42),
    ("Cladophora fascicularis NL", "Pb(II)", 0.0431, 217.99),
    ("Cladophora fascicularis NL", "Pb(II)", 0.0759, 229.92),
    # Distillery sludge
    ("Distillery sludge untreated", "Pb(II)", 0.0188, 71.43),
    ("Distillery sludge autoclaved", "Pb(II)", 0.0170, 81.30),
    ("Distillery sludge HCl", "Pb(II)", 0.0235, 85.47),
    ("Distillery sludge NaOH", "Pb(II)", 0.0219, 91.74),
    ("Distillery sludge HCHO", "Pb(II)", 0.0189, 85.47),
    ("Distillery sludge SDS", "Pb(II)", 0.0228, 74.63),
    # Sour orange
    ("Sour orange residue", "Cu(II)", 0.0660, 52.08),
    ("Sour orange residue NaOH", "Cu(II)", 0.0540, 23.47),
    # Chlorella sorokiniana
    ("Chlorella sorokiniana free", "Cr(III)", 0.1100, 56.56),
    ("Chlorella sorokiniana immobilized", "Cr(III)", 0.3870, 68.51),
    # Rhodococcus opacus
    ("Rhodococcus opacus", "Cr(III)", 0.0015, 714.29),
    # Orange/Grapefruit/Lemon peels
    ("Orange peels", "Cd(II)", 1.779e-3, 123.65),
    ("Orange peels", "Cd(II)", 1.379e-3, 150.63),
    ("Grapefruit peels", "Cd(II)", 5.649e-3, 76.44),
    ("Grapefruit peels", "Cd(II)", 5.382e-3, 110.16),
    ("Lemon peels", "Cd(II)", 2.046e-3, 159.62),
    ("Lemon peels", "Cd(II)", 1.557e-3, 209.08),
    # Fungi
    ("A. niger (fungal)", "Cr(VI)", 0.0026, 17.61),
    ("A. sydoni", "Cr(VI)", 0.0022, 9.07),
    ("P. janthinellum", "Cr(VI)", 0.0042, 9.35),
    # Brown seaweed Cu(II)
    ("Brown seaweed", "Cu(II)", 0.0260, 82.60),
    ("Brown seaweed", "Cu(II)", 0.0240, 88.00),
    ("Brown seaweed", "Cu(II)", 0.0254, 93.90),
    # Bengal gram husk
    ("Bengal gram husk", "Cr(VI)", 0.0090, 91.64),
    # Surfactant modified coconut coir
    ("Surfactant modified coconut coir", "Cr(VI)", 0.1240, 76.30),
    # Ureolytic mixed culture
    ("Living ureolytic mixed culture", "Ni(II)", 8e-5, 12.58),
    ("Non-living ureolytic mixed culture", "Ni(II)", 5.9e-5, 7.41),
    # Vegetal biomass
    ("Vegetal biomass olive pits", "Cd(II)", 0.0440, 9.39),
    # Poplar forests
    ("Poplar forests", "Cu(II)", 1.66e-4, 19.53),
    ("Poplar forests", "Cu(II)", 1.86e-4, 29.76),
    ("Poplar forests", "Cu(II)", 3.574e-4, 29.33),
    # Spirulina
    ("Live spirulina", "Cd(II)", 0.0004, 625.00),
    ("Dead spirulina", "Cd(II)", 0.0013, 355.00),
    # Pristine yeast
    ("Pristine baker's yeast", "Pb(II)", 1.0600, 19.01),
    ("Pristine baker's yeast", "Cd(II)", 1.3200, 3.90),
    # Cystine modified
    ("Cystine modified biomass", "Pb(II)", 1.2500, 45.87),
    ("Cystine modified biomass", "Cd(II)", 1.5200, 11.03),
    # Pretreated A. niger
    ("Pretreated Aspergillus niger", "Cu(II)", 0.1888, 2.61),
    # Lentil shell
    ("Lentil shell", "Cu(II)", 0.0190, 8.98),
    ("Lentil shell", "Cu(II)", 0.0220, 9.51),
    ("Lentil shell", "Cu(II)", 0.0410, 9.59),
    # Wheat shell
    ("Wheat shell", "Cu(II)", 0.0020, 7.39),
    ("Wheat shell", "Cu(II)", 0.0040, 16.08),
    ("Wheat shell", "Cu(II)", 0.0050, 17.42),
    # Rice shell
    ("Rice shell", "Cu(II)", 0.0110, 1.85),
    ("Rice shell", "Cu(II)", 0.0160, 2.31),
    ("Rice shell", "Cu(II)", 0.0180, 2.95),
    # Cladophora Pb(II)
    ("Cladophora fascicularis", "Pb(II)", 0.0359, 198.50),
    ("Cladophora fascicularis", "Pb(II)", 0.0435, 215.80),
    ("Cladophora fascicularis", "Pb(II)", 0.0766, 227.70),
    # Protonated rice bran
    ("Protonated rice bran", "Ni(II)", 9.43e-6, 46.51),
    # Poly(amic acid)
    ("Poly(amic acid) yeast", "Pb(II)", 6.6000, 210.50),
    ("Poly(amic acid) yeast", "Cd(II)", 4.0000, 95.20),
    # Spirogyra
    ("Spirogyra", "Pb(II)", 0.0210, 140.84),
    ("Spirogyra", "Pb(II)", 0.0230, 144.93),
    ("Spirogyra", "Pb(II)", 0.0240, 151.57),
    # Laminaria
    ("Laminaria EC1", "Pb(II)", 0.1568, 352.24),
    ("Laminaria EC2", "Pb(II)", 0.0815, 343.95),
    ("Laminaria KMnO4", "Pb(II)", 0.0109, 391.61),
    ("Laminaria washed", "Pb(II)", 0.0299, 273.50),
    # Peat
    ("Peat", "Pb(II)", 0.4370, 27.80),
    ("Peat", "Pb(II)", 0.4160, 28.30),
    ("Peat", "Pb(II)", 0.4150, 29.50),
    ("Peat", "Pb(II)", 0.3820, 30.80),
    # Bacillus licheniformis
    ("Bacillus licheniformis", "Cr(VI)", 0.0300, 69.35),
    # Lyngbya putealis
    ("Lyngbya putealis", "Cr(VI)", 0.1040, 105.30),
    # Arca shell
    ("Arca shell", "Pb(II)", 0.0690, 30.39),
    ("Arca shell", "Cu(II)", 0.0590, 26.88),
    ("Arca shell", "Ni(II)", 0.2810, 11.75),
    ("Arca shell", "Cs(I)", 0.2280, 4.76),
    ("Arca shell", "Co(II)", 0.2250, 11.53),
    # Maize husk EDTA modified (unmodified has negative KL)
    ("Maize husk EDTA", "Cd(II)", 8.93e-3, 833.33),
    ("Maize husk EDTA", "Pb(II)", 7.36e-3, 714.29),
    ("Maize husk EDTA", "Zn(II)", 0.0400, 769.23),
    # Sargassum filipendula
    ("Sargassum filipendula", "Pb(II)", 4.9600, 1.80),
    # Immobilized yeast
    ("Immobilized S. cerevisiae", "Pb(II)", 0.0140, 30.04),
    # Atlantic Cod
    ("Atlantic Cod fish scale", "As(III)", 5.2000, 0.0248),
    ("Atlantic Cod fish scale", "As(V)", 8.8000, 0.0267),
    # Myriophyllum
    ("Myriophyllum", "Co(II)", 0.3200, 2.30),
    ("Myriophyllum", "Cu(II)", 0.0120, 113.00),
    ("Myriophyllum", "Ni(II)", 0.2500, 3.00),
    ("Myriophyllum", "Zn(II)", 0.8800, 6.80),
    # Bacillus thuringiensis
    ("Bacillus thuringiensis veg", "Ni(II)", 0.0160, 35.46),
    ("Bacillus thuringiensis spore", "Ni(II)", 0.0360, 45.87),
    # Magnetic yeast
    ("Magnetic yeast", "Hg(II)", 0.1470, 48.30),
    ("Magnetic yeast", "Hg(II)", 0.1720, 74.10),
    ("Magnetic yeast", "Hg(II)", 0.1580, 93.40),
    ("Magnetic yeast", "Hg(II)", 0.1190, 133.30),
    # Bacillus jeotgali
    ("Bacillus jeotgali", "Cd(II)", 0.2000, 37.30),
    ("Bacillus jeotgali", "Cd(II)", 0.0700, 47.50),
    ("Bacillus jeotgali", "Cd(II)", 0.0700, 57.90),
    ("Bacillus jeotgali", "Zn(II)", 0.0600, 105.20),
    ("Bacillus jeotgali", "Zn(II)", 0.0800, 222.20),
    # Gelidium
    ("Gelidium", "Cu(II)", 0.015, 33.0),
    ("Gelidium", "Cu(II)", 0.013, 45.0),
    # Algal waste
    ("Algal waste", "Cu(II)", 0.028, 16.7),
    ("Algal waste", "Cu(II)", 0.053, 17.0),
    # Baker's yeast Ni(II)
    ("Baker's yeast", "Ni(II)", 0.2120, 9.01),
    ("Baker's yeast", "Ni(II)", 0.1370, 8.46),
    ("Baker's yeast", "Ni(II)", 0.0940, 7.73),
    ("Baker's yeast", "Ni(II)", 0.0620, 7.37),
    # Cedar sawdust
    ("Cedar sawdust", "Cu(II)", 3.584e-3, 294.12),
    ("Cedar sawdust", "Cu(II)", 6.290e-3, 144.93),
    ("Cedar sawdust", "Cu(II)", 7.750e-3, 106.38),
    # Crushed brick
    ("Crushed brick", "Cu(II)", 5.274e-3, 153.85),
    ("Crushed brick", "Cu(II)", 6.935e-3, 104.17),
    ("Crushed brick", "Cu(II)", 10.243e-3, 68.03),
    # Caladium bicolor Pb(II)
    ("Caladium bicolor", "Pb(II)", 0.0393, 37.17),
    ("Caladium bicolor", "Pb(II)", 0.0380, 52.63),
    ("Caladium bicolor", "Pb(II)", 0.0263, 37.59),
    ("Caladium bicolor", "Pb(II)", 0.0260, 20.45),
    ("Caladium bicolor", "Pb(II)", 0.0250, 8.81),
    ("Caladium bicolor", "Pb(II)", 0.0120, 5.49),
    # Caladium bicolor Cd(II)
    ("Caladium bicolor", "Cd(II)", 0.0289, 42.19),
    ("Caladium bicolor", "Cd(II)", 0.0253, 59.02),
    ("Caladium bicolor", "Cd(II)", 0.0208, 46.95),
    ("Caladium bicolor", "Cd(II)", 0.0141, 22.47),
    ("Caladium bicolor", "Cd(II)", 0.0124, 14.03),
    ("Caladium bicolor", "Cd(II)", 0.0195, 8.58),
    # Pseudomonas Co(II)
    ("Pseudomonas", "Co(II)", 0.1884, 24.75),
    ("Pseudomonas", "Co(II)", 0.2393, 22.39),
    ("Pseudomonas", "Co(II)", 0.3699, 22.39),
    ("Pseudomonas", "Co(II)", 0.4123, 22.98),
    ("Pseudomonas", "Co(II)", 0.5125, 22.39),
    # Pseudomonas Zn(II)
    ("Pseudomonas", "Zn(II)", 0.3593, 18.43),
    ("Pseudomonas", "Zn(II)", 0.4327, 18.43),
    ("Pseudomonas", "Zn(II)", 0.5398, 18.43),
    ("Pseudomonas", "Zn(II)", 0.7141, 18.43),
    ("Pseudomonas", "Zn(II)", 0.8272, 19.06),
    # Pseudomonas Pb(II)
    ("Pseudomonas", "Pb(II)", 0.2264, 45.58),
    ("Pseudomonas", "Pb(II)", 0.2639, 43.51),
    ("Pseudomonas", "Pb(II)", 0.3175, 39.37),
    ("Pseudomonas", "Pb(II)", 0.3778, 37.29),
    ("Pseudomonas", "Pb(II)", 0.4020, 41.44),
    # Pseudomonas Hg(II)
    ("Pseudomonas", "Hg(II)", 0.0862, 80.24),
    ("Pseudomonas", "Hg(II)", 0.1176, 76.22),
    ("Pseudomonas", "Hg(II)", 0.1500, 76.22),
    ("Pseudomonas", "Hg(II)", 0.1949, 78.23),
    ("Pseudomonas", "Hg(II)", 0.2353, 84.25),
    # Sargassum wightii
    ("Sargassum wightii", "Ni(II)", 0.0054, 63.20),
    ("Sargassum wightii", "Ni(II)", 0.0055, 76.10),
    ("Sargassum wightii", "Ni(II)", 0.0065, 81.20),
    ("Sargassum wightii", "Ni(II)", 0.0057, 79.10),
    # Lyngbya exopolysaccharides
    ("Lyngbya exopolysaccharides", "Cr(VI)", 0.1890, 178.60),
    # Turbinaria conoides
    ("Turbinaria conoides", "Pb(II)", 0.0410, 385.10),
    ("Turbinaria conoides", "Pb(II)", 0.0440, 420.10),
    ("Turbinaria conoides", "Pb(II)", 0.0480, 439.40),
    ("Turbinaria conoides", "Pb(II)", 0.0570, 228.00),
    ("Turbinaria conoides", "Pb(II)", 0.0470, 401.40),
    # Chlorella vulgaris Cd/Ni
    ("Chlorella vulgaris", "Cd(II)", 28.3000, 0.02),
    ("Chlorella vulgaris", "Ni(II)", 70.9000, 0.03),
    # Treated sawdust
    ("Treated sawdust", "Cr(VI)", 2.3600, 3.60),
    # PVA-Sargassum
    ("PVA-Sargassum", "Cu(II)", 28.0700, 0.21),
    ("Freely suspended Sargassum", "Cu(II)", 2.8600, 0.96),
    # Rhyzopus oryzae
    ("Rhyzopus oryzae viable", "Cu(II)", 0.1290, 19.40),
    ("Rhyzopus oryzae NaOH", "Cu(II)", 0.0280, 43.70),
    # Gelidium Zn/Cr
    ("Gelidium", "Zn(II)", 2.6, 13.0),
    ("Gelidium", "Cr(III)", 2.1, 18.0),
    # Algal waste Zn/Cr
    ("Algal waste", "Zn(II)", 4.2, 7.1),
    ("Algal waste", "Cr(III)", 3.3, 11.8),
    # Nopal
    ("Nopal", "Pb(II)", 2.868e-4, 5.43),
    ("Nopal", "Pb(II)", 6.118e-4, 11.29),
    ("Nopal", "Pb(II)", 20.41e-4, 12.76),
    ("Nopal", "Pb(II)", 22.80e-4, 13.39),
    ("Nopal", "Pb(II)", 9.894e-4, 12.34),
    # Lignin
    ("Lignin", "Pb(II)", 0.2081, 89.51),
    ("Lignin", "Cu(II)", 0.4309, 22.87),
    ("Lignin", "Cd(II)", 0.1421, 25.40),
    ("Lignin", "Zn(II)", 0.2046, 11.25),
    ("Lignin", "Ni(II)", 0.2511, 5.99),
    # Azadirachta
    ("Azadirachta indica bark", "Zn(II)", 0.0572, 33.49),
    # Water plants Cr
    ("Water lily", "Cr(III)", 0.0310, 6.11),
    ("Water lily", "Cr(VI)", 0.0270, 5.11),
    ("Water hyacinth", "Cr(III)", 0.0500, 6.61),
    ("Water hyacinth", "Cr(VI)", 0.9170, 0.34),
    ("Green taro", "Cr(III)", 0.0340, 6.07),
    ("Green taro", "Cr(VI)", 0.0420, 1.42),
    ("Mangrove leaves", "Cr(III)", 0.0320, 6.54),
    ("Mangrove leaves", "Cr(VI)", 0.0440, 5.72),
    ("Reed mat", "Cr(III)", 0.0420, 7.18),
    ("Reed mat", "Cr(VI)", 0.0390, 1.66),
    # Agaricus
    ("Agaricus bisporus", "Cr(VI)", 0.0700, 8.00),
    ("Agaricus bisporus", "Cr(VI)", 0.0500, 11.28),
    ("Agaricus bisporus", "Cr(VI)", 0.0900, 13.79),
    # Sugar beet
    ("Sugar beet pulp", "Pb(II)", 0.4759, 0.37),
    ("Sugar beet pulp", "Cd(II)", 0.1237, 0.13),
    # Cladonia
    ("Cladonia rangiformis", "Cu(II)", 0.1051, 7.69),
    # Wheat bran
    ("Wheat bran", "Cu(II)", 2.9600, 12.58),
    ("Wheat bran", "Cu(II)", 7.5700, 10.61),
    ("Wheat bran", "Zn(II)", 4.9500, 16.02),
    ("Wheat bran", "Zn(II)", 5.4500, 6.67),
    # Phosphorylated orange waste
    ("Phosphorylated orange waste", "In(III)", 0.1220, 81.22),
    ("Phosphorylated orange waste", "Ga(III)", 0.2712, 49.01),
    ("Phosphorylated orange waste", "Cu(II)", 0.0452, 67.35),
    ("Phosphorylated orange waste", "Pb(II)", 0.0367, 251.05),
    ("Phosphorylated orange waste", "Fe(III)", 0.1201, 173.14),
    ("Phosphorylated orange waste", "Zr(IV)", 0.0361, 113.11),
    ("Phosphorylated orange waste", "As(III)", 0.0174, 91.40),
    ("Phosphorylated orange waste", "As(V)", 0.0358, 73.42),
    # Ficus
    ("Ficus religiosa", "Cr(VI)", 0.0120, 26.25),
    ("Ficus religiosa", "Pb(II)", 0.0220, 37.45),
    # Birch wood etc
    ("Birch wood", "Cu(II)", 0.2000, 4.90),
    ("Fucus vesiculosus", "Cu(II)", 1.1000, 23.40),
    ("Pleurozium schreberi", "Cu(II)", 1.1500, 11.10),
    # Cupriavidus
    ("Cupriavidus", "Pb(II)", 0.0211, 50.10),
    ("Cupriavidus", "Cu(II)", 0.0363, 19.00),
    ("Cupriavidus", "Cd(II)", 0.0167, 19.60),
    # Mimosa inoculated
    ("Mimosa pudica inoculated", "Pb(II)", 0.0383, 485.00),
    ("Mimosa pudica inoculated", "Cu(II)", 0.0990, 25.40),
    ("Mimosa pudica inoculated", "Cd(II)", 0.0307, 42.90),
    # Mimosa without
    ("Mimosa pudica w/o", "Pb(II)", 0.0413, 26.10),
    ("Mimosa pudica w/o", "Cu(II)", 0.0614, 22.70),
    ("Mimosa pudica w/o", "Cd(II)", 0.0392, 25.30),
    # Reed
    ("Reed", "Pb(II)", 1.49e-3, 17.1552),
    # Marine algae DP95Ca (pH1 has KL=0, excluded)
    ("Marine algae DP95Ca", "Pb(II)", 0.1252, 158.99),
    ("Marine algae DP95Ca", "Pb(II)", 4.3258, 269.88),
    ("Marine algae DP95Ca", "Pb(II)", 11.8159, 307.54),
    ("Marine algae DP95Ca", "Pb(II)", 23.7369, 324.28),
    ("Marine algae DP95Ca", "Cu(II)", 0.0031, 11.44),
    ("Marine algae DP95Ca", "Cu(II)", 1.1457, 62.90),
    ("Marine algae DP95Ca", "Cu(II)", 1.7485, 76.88),
    ("Marine algae DP95Ca", "Cu(II)", 2.6251, 82.60),
    # Marine algae ER95Ca
    ("Marine algae ER95Ca", "Pb(II)", 0.0425, 10.46),
    ("Marine algae ER95Ca", "Pb(II)", 0.4421, 87.87),
    ("Marine algae ER95Ca", "Pb(II)", 1.2413, 207.12),
    ("Marine algae ER95Ca", "Pb(II)", 1.7155, 244.78),
    ("Marine algae ER95Ca", "Pb(II)", 1.8584, 263.60),
    ("Marine algae ER95Ca", "Cu(II)", 0.0976, 4.45),
    ("Marine algae ER95Ca", "Cu(II)", 0.0630, 28.59),
    ("Marine algae ER95Ca", "Cu(II)", 0.6295, 60.36),
    ("Marine algae ER95Ca", "Cu(II)", 1.2890, 67.35),
    ("Marine algae ER95Ca", "Cu(II)", 2.0696, 70.53),
    # Ca-alginate
    ("Ca-alginate resin", "Pb(II)", 0.8730, 670.20),
    # Rice husk
    ("Rice husk", "Cu(II)", 0.1000, 29.00),
    ("Rice husk", "Cu(II)", 0.0500, 22.00),
    ("Rice husk", "Cu(II)", 0.0300, 18.00),
    ("Rice husk", "Pb(II)", 0.0900, 108.00),
    ("Rice husk", "Pb(II)", 0.0200, 105.00),
    ("Rice husk", "Pb(II)", 0.0100, 96.00),
    # Fucus
    ("Fucus vesiculosos", "Cr(III)", 0.0362, 62.91),
    ("Fucus vesiculosos", "Cr(VI)", 0.0339, 42.63),
    ("Fucus spiralis", "Cr(III)", 0.0340, 60.83),
    ("Fucus spiralis", "Cr(VI)", 0.0283, 35.35),
    # Ulva
    ("Ulva lactuca", "Cr(III)", 0.0381, 36.91),
    ("Ulva lactuca", "Cr(VI)", 0.0379, 27.55),
    ("Ulva spp.", "Cr(III)", 0.0265, 53.03),
    ("Ulva spp.", "Cr(VI)", 0.0231, 30.15),
    # Palmaria
    ("Palmaria palmata", "Cr(III)", 0.0950, 29.63),
    ("Palmaria palmata", "Cr(VI)", 0.1662, 33.79),
    # Polysiphonia
    ("Polysiphonia lanosa", "Cr(III)", 0.0258, 33.79),
    ("Polysiphonia lanosa", "Cr(VI)", 0.0469, 45.75),
    # Termitomyces: SKIP (qmax missing in table, 0.998 is R2)
    # Activated sludge
    ("Activated sludge", "Cd(II)", 0.3200, 28.10),
    ("Activated sludge", "Cu(II)", 0.1201, 19.06),
    ("Activated sludge", "Ni(II)", 0.1697, 7.78),
    ("Activated sludge", "Pb(II)", 0.1000, 142.96),
    ("Activated sludge", "Zn(II)", 0.6098, 15.69),
    # Grape bagasse
    ("Grape bagasse", "Cd(II)", 0.0527, 53.84),
    ("Grape bagasse", "Pb(II)", 0.0074, 42.27),
    # Pinus
    ("Pinus sylvestris", "Cr(VI)", 0.0450, 238.10),
    ("Pinus sylvestris", "Cr(VI)", 0.1165, 263.16),
    ("Pinus sylvestris", "Cr(VI)", 0.2166, 256.41),
    # Pseudomonas veronii
    ("Pseudomonas veronii", "Cd(II)", 0.0526, 54.00),
]

# ============================================================
# Calculate
# ============================================================
valid = [(n, s, kl, qm) for n, s, kl, qm in data if kl > 0 and qm > 0]
n = len(valid)

results = []
for name, solute, KL, qmax in valid:
    kappa = 1.0 / (KL * qmax)  # m/V = 1 g/L
    R = kappa / (1 + kappa)
    results.append((name, solute, KL, qmax, kappa, R))

results.sort(key=lambda x: x[4])
kappas = [r[4] for r in results]

print(f"Total valid systems: {n}")
print(f"Source: Table 2 of Febrianto et al. (2009) J. Hazard. Mater. 162, 616-645")
print(f"Excluded: negative KL (5), zero KL (2), missing qmax (1)")
print(f"\n{'='*65}")
print(f"κ distribution ({n} systems, m/V = 1 g/L)")
print(f"{'='*65}")

bins = [0, 0.001, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, float('inf')]
labels = ['<0.001', '0.001-0.01', '0.01-0.05', '0.05-0.1', '0.1-0.5', '0.5-1', '1-2', '2-5', '≥5']
for i in range(len(bins)-1):
    count = sum(1 for k in kappas if bins[i] <= k < bins[i+1])
    pct = count/n*100
    bar = '#' * int(pct/2)
    print(f"  {labels[i]:>12}: {count:3d} ({pct:5.1f}%) {bar}")

n_lt_001 = sum(1 for k in kappas if k < 0.01)
n_lt_01 = sum(1 for k in kappas if k < 0.1)
n_lt_1 = sum(1 for k in kappas if k < 1)
n_ge_1 = sum(1 for k in kappas if k >= 1)

print(f"\nκ < 0.01:  {n_lt_001:3d} ({n_lt_001/n*100:.1f}%)  → R < ~1%  (deep purification)")
print(f"κ < 0.1:   {n_lt_01:3d} ({n_lt_01/n*100:.1f}%)  → R < ~10%")
print(f"κ < 1:     {n_lt_1:3d} ({n_lt_1/n*100:.1f}%)  → R < 50%")
print(f"κ ≥ 1:     {n_ge_1:3d} ({n_ge_1/n*100:.1f}%)  → R ≥ 50%  (unsuitable for deep purification)")

print(f"\n{'='*65}")
print(f"Systems with κ < 0.01 (capable of <1% residual)")
print(f"{'='*65}")
for name, solute, KL, qmax, k, R in results:
    if k < 0.01:
        print(f"  {name + ' ' + solute:50s}: KL={KL:10.4f}, qmax={qmax:8.2f}, κ={k:.6f}, R={R:.4f} ({R:.2%})")

# Save CSV
with open('SI_S9_biosorption_kappa_distribution.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['Biosorbent', 'Solute', 'K_L', 'qmax', 'kappa', 'R_exact', 'R_percent'])
    for name, solute, KL, qmax, k, R in results:
        w.writerow([name, solute, KL, qmax, f'{k:.6f}', f'{R:.6f}', f'{R*100:.2f}%'])

print(f"\nResults saved to SI_S9_biosorption_kappa_distribution.csv")
