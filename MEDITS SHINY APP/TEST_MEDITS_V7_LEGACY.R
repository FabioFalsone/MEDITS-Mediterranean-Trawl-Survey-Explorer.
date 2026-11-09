# MEDITS Ocean v7 - independent smoke checks for map and LFD plotting.
# Run from the extracted application folder: source('TEST_MEDITS_V7.R')
# Requires R packages already used by the MEDITS app: shiny and ggplot2.
for (f in c('Medits_app_GSA.R','Medits_function_GSAMedits_gsa.R',
            'medits_display_helpers.R','medits_qc_helpers.R',
            'medits_report_helpers.R')) {
  stopifnot(file.exists(f))
  parse(file=f)
}
stopifnot(file.exists('www/medits_theme.css'))
source('medits_display_helpers.R')

lfd <- data.frame(Year=c(2023,2023,2024,2024,2025),
                  Length=c(120,125,120,130,135),
                  Frequency=c(5,8,0,12,9))
map <- data.frame(YEAR=c(2023,2024,2025),
                  LAT=c(36.1,36.2,36.3),LONG=c(13.1,13.2,13.3),
                  DEPTH=c(130,300,520),BIOMASS=c(0,5,1),DENSITY=c(0,120,35))
example <- list(LFD=lfd,LFDF=lfd,LFDM=lfd,LFDU=lfd,
                LFD_shelf=lfd,LFDF_shelf=lfd,LFDM_shelf=lfd,LFDU_shelf=lfd,
                LFD_slope=lfd,LFDF_slope=lfd,LFDM_slope=lfd,LFDU_slope=lfd,
                BIbubble_map=map,DIbubble_map=map,report_gsas=c(15,16))
for (view in c('total','shelf','slope')) {
  for (sex in c('C','F','M','U')) {
    p <- medits_lfd_display_plot(example,view,sex)
    stopifnot(inherits(p,'ggplot'))
    invisible(ggplot2::ggplot_build(p))
  }
  for (metric in c('BIOMASS','DENSITY')) {
    d <- medits_spatial_prepare(example,metric,view)
    stopifnot(is.data.frame(d),nrow(d)>0,all(is.finite(d$log_value)))
  }
}
stopifnot(nrow(medits_spatial_prepare(example,'BIOMASS','total'))==3L)
stopifnot(nrow(medits_spatial_prepare(example,'BIOMASS','shelf'))==1L)
stopifnot(nrow(medits_spatial_prepare(example,'DENSITY','slope'))==2L)
cat('PASS: MEDITS v7 map/LFD smoke tests completed.\n')
cat('NOTE: this verifies display helpers with synthetic data, not the full assessment with TA/TB/TC.\n')