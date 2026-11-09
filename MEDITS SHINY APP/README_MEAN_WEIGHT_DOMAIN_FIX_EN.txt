MEDITS Ocean v9 - Mean Weight Bathymetric Domain Fix
===================================================

Issue fixed
The previous Mweight output from the Mannini et al. pipeline used ALL species TB
records irrespective of the selected graph domain, so changing Total / Shelf /
Slope changed biomass and density but not mean weight.

Fix
- The original Mweight output remains available and unchanged.
- A new Mweight_domain result is calculated in run_medits_from_raw.
- Its TB records are matched ONCE to validated TA hauls via the existing
  country + GSA + year + haul key, then assigned to depth stratum A-E.
- Total comprises strata A-E; Shelf A-C; Slope D-E, matching the same
  bathymetric classification already used by MEDITS for this survey.
- Mean weight (g) = sum of haul total catch weight in grams / sum of
  corresponding haul total number of individuals. This is a descriptive
  pooled-individual mean, NOT biomass-index / density-index.
- Missing or zero-total-number domains are NA (not artificial zeros).
- A diagnostic column compares the valid-haul Total with the historical
  full-TB value, helping detect unmatched hauls without changing data.

Views and exports
- Single-GSA Mean weight plot: now responds to Total/Shelf/Slope.
- Mean weight PNG and CSV: show/export the selected domain.
- Three-panel GSA comparison, PNG and CSV: use domain-specific Mean weight.
- Complete ZIP with all GSA indices: exports Mweight_domain alongside Mweight.
- PDF report chart, tabular summary and scientific description: updated.
- Biomass and density calculations, LFD, ALK, QC and other maps are unchanged.

Scientific interpretation
Changing the depth view will not necessarily remove an anomalous peak.
In the supplied 2016-2018 Malta GSA15 ARIS/FOL data all valid positive-catch
hauls are in the slope, so the approximately 41.69 g peak in 2017 persists.
This is a real value in the supplied TB data, not proof that it is biologically
valid. Check Malta records before any manual correction.

Before running
1. Back up the previously working MEDITS Ocean v9 folder.
2. Extract this package to a NEW folder, copy the three reference CSVs.
3. In RStudio run source('TEST_MEAN_WEIGHT_DOMAIN_FIX.R') for source-code
   parsing checks (requires a local R installation).
4. Run the app and re-run MEDITS analysis once to build Mweight_domain.
5. Compare Total/Shelf/Slope in both the individual plot and GSA comparison.