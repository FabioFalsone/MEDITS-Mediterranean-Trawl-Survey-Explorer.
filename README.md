# MEDITS - Mediterranean Trawl Survey Explorer

R Shiny application for exploring **MEDITS Mediterranean trawl survey** information, including biomass and density indices, GSA comparisons, spatial maps, length-frequency distributions (LFD), age–length keys (ALK), and biological quality-control checks.

## Application

The application is located in **[`MEDITS SHINY APP/`](MEDITS%20SHINY%20APP/)**.

### Run locally

1. Install R and the packages imported by the R scripts (including FLR dependencies where required).
2. Supply the lookup files specified by the application locally; do not commit third-party survey data.
3. In RStudio, open `MEDITS SHINY APP/Medits_app_GSA.R` or run:

```r
shiny::runApp("MEDITS SHINY APP", launch.browser = TRUE)
```

4. Upload your own TA, TB, and TC input files inside the app. An optional age-reading CSV may be used for ALK.

## Features

- Aggregated and individual GSA trends, including total, shelf and slope domains.
- Biomass, density and mean-weight series, spatial maps, LFDs and optional ALKs.
- Manual biological checks of lengths, weights and ages, with correction logs.
- Exports to plots, tables and reports.
- No bootstrap computations in this version.

## Methods and attribution

This application builds upon MEDITS analytical methods and R scripts developed by Alessandro Mannini at the Joint Research Centre (JRC) of the European Commission. The original scientific and software contributions are acknowledged. Modifications and additional Shiny interfaces, visualisations, and quality-control tools belong to this application.

**Scientific reference**

Mannini, A. (2020). *The JRC MEDITS R script*. Publications Office of the European Union. JRC119776. https://doi.org/10.2760/5799

**Source-code rights:** Citation does not replace license compliance. Original third-party code may carry GNU GPL license conditions; users should preserve original copyright/license notices and review code provenance before redistribution.

## Data and licensing

**No survey TA/TB/TC datasets are included.** Keep raw observations, identifiable site information, and licensed third-party reference files out of public repositories. No license is asserted for any third-party scientific code; confirm rights before publishing publicly.