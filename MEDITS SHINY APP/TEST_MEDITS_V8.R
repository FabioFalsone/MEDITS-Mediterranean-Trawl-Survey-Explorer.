# Synthetic-only checks: no real MEDITS data or original inputs are modified.
source('medits_display_helpers.R',local=FALSE)
stopifnot(requireNamespace('ggplot2',quietly=TRUE),
          requireNamespace('ggpubr',quietly=TRUE),
          requireNamespace('shiny',quietly=TRUE))
get_data <- function(f) data.frame(Year=rep(2024:2025,each=3),
  Length=rep(c(20,25,30),2),Frequency=rep(f,2))
a <- list(
  LFD=get_data(c(20,30,40)),
  LFDF=get_data(c(8,11,20)),
  LFDM=get_data(c(7,12,16)),
  LFDU=get_data(c(5,7,4)),
  SR_vec=data.frame(Length=c(20,25,30),sr=c(.3,.5,.8)),
  report_gsas=c('15','16'))
for(view in c('total','shelf','slope')) {
  if(view!='total') {
    for(k in c('LFD','LFDF','LFDM','LFDU')) a[[paste0(k,'_',view)]] <- a[[k]]
  }
  z <- medits_lfd_allocate(a,view)
  stopifnot(isTRUE(z$allocated),nrow(z$audit)==6L,
    all(abs(z$audit$F_corrected+z$audit$M_corrected-z$audit$Total)<1e-10),
    all(abs(z$audit$F_allocated+z$audit$M_allocated-z$audit$U_obs)<1e-10))
  for(sex in c('C','F','M','U')) {
    p <- medits_lfd_display_plot(a,view,sex)
    stopifnot(inherits(p,'ggplot'))
    ggplot2::ggplot_build(p)
  }
}
ages <- data.frame(COUNTRY='ITA',AREA=rep(c(15,16),each=5),YEAR=2025,
    SPECIES='HKE',SEX=c('F','M','M','F','U','F','F','M','F','M'),
    LENGTH_CLASS=c(20,20,25,25,30,20,20,25,25,30),
    AGE=c(1,2,2,3,4,1,2,2,3,4))
aged <- medits_alk_prepare(ages,species='HKE',gsa=c('15','16'),country='ITA')
for(sex in c('ALL','F','M')) {
  x <- medits_alk_probabilities(aged,sex=sex,year='2025')
  s <- aggregate(Probability~Length,data=x,FUN=sum)
  stopifnot(all(abs(s$Probability-1)<1e-10))
  ggplot2::ggplot_build(medits_alk_plot(x,sex))
}
map <- data.frame(YEAR=2025,LONG=c(12,12.2),LAT=c(35.8,36),DEPTH=c(100,350),
   SQKM=c(.01,.01),BIOMASS=c(0,35),DENSITY=c(0,150))
a$BIbubble_map <- map
stopifnot(nrow(medits_spatial_prepare(a,'BIOMASS','shelf'))==1L,
          nrow(medits_spatial_prepare(a,'BIOMASS','slope'))==1L,
          nrow(medits_spatial_prepare(a,'DENSITY','total'))==2L)
cat('MEDITS V8 synthetic display / ALK / allocation checks: OK\n')