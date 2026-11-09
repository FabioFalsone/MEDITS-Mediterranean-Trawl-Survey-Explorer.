MEDITS OCEAN v9 – Original v4 LFD graphics restored
=====================================================

This version starts from the MEDITS Ocean v8 ZIP that was confirmed working.
Only the display of the four LFD outputs was changed.

Precisely restored from the original v4 application:
- Total, Female and Male: ggplotFL::plotLengths() graphs in two year blocks,
  stacked with ggpubr::ggarrange(), blue vertical line for mean length,
  red vertical line for L95, original faceted formatting.
- Undetermined: original ggplot2 bar plots faceted by year.
- Screen heights and PNG export sizes are brought closer to the original.

Retained from v8:
- Corrected or observed Female/Male LFD toggle (default LOESS corrected).
- The global LOESS curve from the existing assessment and the balancing audit.
- ALK from age-readings CSV, marine-themed interface, colour-enhanced maps,
  survey indices and their calculation functions, PDFs and CSV downloads.
- No automatic biological data modification.

An exceptional failure inside ggplotFL falls back to v8's safe plot renderer
for that selection only, and R prints a descriptive warning. This prevents
one sparse year from crashing the entire app, while ordinary cases use the
real original plotLengths design.

INSTALLATION
1. Unzip to a NEW folder. Do not merge into the working v8 folder.
2. Copy Sp_Medits.csv, ASFIS_2020.csv and MEDITS_Strata.csv from your setup.
3. From that folder in RStudio, run source('TEST_MEDITS_V9.R').
4. Run shiny::runApp('Medits_app_GSA.R',launch.browser=TRUE).

VALIDATION
R is not installed in the file-creation environment, so full rendering
could not be executed here. TEST_MEDITS_V9.R exercises the true legacy
plotLengths code with synthetic frequencies for Total, Shelf, Slope, sex
correction, undetermined, and PNG output. Compare actual outputs with v4.