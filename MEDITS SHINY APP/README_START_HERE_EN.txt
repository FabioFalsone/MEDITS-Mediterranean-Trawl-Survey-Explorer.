MEDITS Ocean Dashboard v8 – LFD two-block restoration, LOESS sex allocation,
contrasting maps, and optional ALK.

Based strictly on the confirmed working MEDITS Ocean V7 package.
The calculation engine Medits_function_GSAMedits_gsa.R is unchanged.

WHAT IS NEW
1. Length-frequency plots: original two-block annual layout is restored, with
   blue vertical lines for mean length and red lines for the 95th percentile
   (L95). The plotting code uses the standardized LFD data tables directly to
   avoid earlier failures in FLQuant conversion.
2. Female/male LFD: choose between observed frequencies and LOESS-corrected
   frequencies. Corrected female = observed female + observed undetermined *
   pooled LOESS SR(length). Corrected male = observed male + observed
   undetermined * (1 - pooled LOESS SR(length)). This same pooled curve is
   reused for all annual series and depth domains, and for independent GSAs.
   The original total, F, M, U model objects remain unchanged.
3. A warning appears if observed female+male exceeds total for any length-year
   cell. Allocation audit CSVs preserve this information. When the pooled LOESS
   curve is unavailable, the app displays observed F/M without allocating U.
4. Map points: bigger circles, dark outlines and a higher-contrast blue/teal/
   yellow/orange/red scale; log1p biomass and density are unchanged.
5. Optional ALK tab: upload age-readings CSV in sidebar and use the ALK tab
   under Length & sex. Displays observed P(Age | Length) for all sexes, females
   or males; can pool years or select one year. ALK is based on observed aged
   individuals; it is NOT raised to area, and it does NOT modify MEDITS LFDs.

AGE READINGS CSV
Use AGE_READINGS_TEMPLATE.csv. Required: AGE, LENGTH_CLASS (or LENGTH), COUNTRY, AREA, YEAR, SPECIES
(3-alpha e.g. HKE). Recommended: SEX (F/M/U) for sex-specific ALKs.
Each row represents one aged individual, unless COUNT is provided explicitly.
Length must be in the intended unit, normally MEDITS mm. Age must be integer.
The ALK is filtered by selected GSA, country and species when those columns exist.
Do not equate AGE and MATURITY; they are different biological quantities.

INSTALLATION
1. Extract this ZIP into a NEW folder (do not overwrite your working v7).
2. Copy Sp_Medits.csv, ASFIS_2020.csv and MEDITS_Strata.csv from your local
   working MEDITS installation into the extracted folder.
3. Open Medits_app_GSA.R in RStudio, press Run App.
4. Run source('TEST_MEDITS_V8.R') before processing your own survey data.
5. Verify annual biomass/density indices against the working v7; these must
   match exactly. Verify F_corrected + M_corrected against Total in QC CSV.

LIMITS
The ZIP has been inspected statically. R and the user survey CSVs are not
available in the execution environment, so runtime and output agreement must
be verified in RStudio. Only the plotting/helper layers were changed.