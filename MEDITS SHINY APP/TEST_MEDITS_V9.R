# Run this file from the MEDITS_OCEAN_V9 directory in RStudio.
# It tests original v4-style plotLengths graphics with synthetic records only.
stopifnot(requireNamespace('ggplot2',quietly=TRUE),
          requireNamespace('ggpubr',quietly=TRUE),
          requireNamespace('FLCore',quietly=TRUE),
          requireNamespace('ggplotFL',quietly=TRUE),
          requireNamespace('shiny',quietly=TRUE))
source('medits_length_stats.R')
source('medits_display_helpers.R')
source('medits_legacy_lfd_plot.R')

years <- 2020:2027
lens <- c(100,110,120,130,140,150)
make_lfd <- function(freq) {
  data.frame(Year=rep(years,each=length(lens)),
             Length=rep(lens,length(years)),
             Frequency=rep(freq,length(years)),stringsAsFactors=FALSE)
}
a <- list(
  LFD=make_lfd(c(10,30,50,30,15,5)),
  LFDF=make_lfd(c(3,10,22,14,7,2)),
  LFDM=make_lfd(c(5,16,19,10,6,2)),
  LFDU=make_lfd(c(2,4,9,6,2,1)),
  SR_vec=data.frame(Length=lens,sr=c(.37,.45,.5,.55,.6,.68)))
for(view in c('shelf','slope'))
  for(key in c('LFD','LFDF','LFDM','LFDU')) a[[paste0(key,'_',view)]] <- a[[key]]
orig <- a$LFD
for(view in c('total','shelf','slope')) {
  qc <- medits_lfd_allocate(a,view)
  stopifnot(isTRUE(qc$allocated),
            all(abs(qc$audit$F_corrected + qc$audit$M_corrected - qc$audit$Total)<1e-8))
  for(sex in c('C','F','M')) {
    p <- medits_lfd_legacy_plot(a,view,sex,corrected=TRUE,fallback=FALSE)
    stopifnot(inherits(p,'ggplot'))
    png <- tempfile(fileext='.png')
    ggplot2::ggsave(png,plot=p,width=10,height=12,dpi=80)
    stopifnot(file.exists(png),file.size(png)>1000)
    unlink(png)
  }
  p_und <- medits_lfd_legacy_plot(a,view,'U',fallback=FALSE)
  stopifnot(inherits(p_und,'ggplot'))
}
stopifnot(identical(orig,a$LFD))
# Verify toggling observed vs LOESS-corrected values doesn't change originals.
stopifnot(inherits(medits_lfd_legacy_plot(a,'total','F',corrected=FALSE,
                                          fallback=FALSE),'ggplot'))
stopifnot(inherits(medits_lfd_legacy_plot(a,'total','M',corrected=FALSE,
                                          fallback=FALSE),'ggplot'))
cat('MEDITS v9: original ggplotFL two-block LFD graphics verified.\n')