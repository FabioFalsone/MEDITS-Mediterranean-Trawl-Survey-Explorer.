# MEDITS display-only helpers. No biological estimates are changed.
# These functions are used by the Shiny app and by TEST_MEDITS_V7.R.

# Reliable LFD drawing directly from original standardised result tables.
# FLQuant conversion is intentionally bypassed for the UI only: no model values
# are re-raised, reallocated or otherwise altered. Zero and one-bin years safe.
# Exact observed LFDs are preserved. Allocation is display/export-only.
# Year/Length values are already standardised by run_medits_from_raw().
medits_lfd_standard_table <- function(a,view='total',sex='C') {
  k <- switch(sex,C='LFD',F='LFDF',M='LFDM',U='LFDU')
  if(view %in% c('shelf','slope')) k <- paste0(k,'_',view)
  d <- a[[k]]
  empty <- data.frame(Year=integer(),Length=numeric(),Frequency=numeric())
  if(!is.data.frame(d) || !nrow(d) ||
     !all(c('Year','Length','Frequency') %in% names(d))) return(empty)
  x <- data.frame(Year=suppressWarnings(as.integer(as.character(d$Year))),
                  Length=suppressWarnings(as.numeric(as.character(d$Length))),
                  Frequency=suppressWarnings(as.numeric(as.character(d$Frequency))))
  x <- x[is.finite(x$Year) & is.finite(x$Length) & x$Length>0 &
           is.finite(x$Frequency) & x$Frequency>=0,,drop=FALSE]
  if(!nrow(x)) return(empty)
  z <- stats::aggregate(Frequency~Year+Length,data=x,FUN=sum)
  z[order(z$Year,z$Length),,drop=FALSE]
}

# Pooled LOESS curve, already fitted in the assessment calculation and
# shared by the independently run GSAs. Do not refit it per year/domain.
medits_lfd_allocate <- function(a,view='total') {
  whole <- medits_lfd_standard_table(a,view,'C')
  F <- medits_lfd_standard_table(a,view,'F')
  M <- medits_lfd_standard_table(a,view,'M')
  U <- medits_lfd_standard_table(a,view,'U')
  if(!nrow(whole)) return(list(total=whole,F_observed=F,M_observed=M,
                                U_observed=U,F_corrected=F,M_corrected=M,
                                audit=data.frame(),allocated=FALSE,
                                message='No total LFD for this domain'))
  tab <- whole; names(tab)[3] <- 'Total'
  for(z in list(list(F,'F_obs'),list(M,'M_obs'),list(U,'U_obs'))) {
    if(nrow(z[[1]])) {
      sub <- z[[1]]; names(sub)[3] <- z[[2]]
      tab <- merge(tab,sub,by=c('Year','Length'),all.x=TRUE,sort=FALSE)
    } else tab[[z[[2]]]] <- 0
  }
  for(nm in c('F_obs','M_obs','U_obs')) tab[[nm]][is.na(tab[[nm]])] <- 0
  sr <- a$SR_vec
  usable <- is.data.frame(sr) && all(c('Length','sr') %in% names(sr)) &&
    length(unique(suppressWarnings(as.numeric(as.character(sr$Length)))[
      is.finite(suppressWarnings(as.numeric(as.character(sr$Length)))) &
      is.finite(suppressWarnings(as.numeric(as.character(sr$sr))))])) >= 2L
  if(usable) {
    x <- suppressWarnings(as.numeric(as.character(sr$Length)))
    y <- suppressWarnings(as.numeric(as.character(sr$sr)))
    ok <- is.finite(x) & is.finite(y)
    curve <- stats::aggregate(y[ok],by=list(Length=x[ok]),FUN=mean)
    names(curve)[2] <- 'SR'
    curve <- curve[order(curve$Length),,drop=FALSE]
    tab$SR <- stats::approx(curve$Length,curve$SR,xout=tab$Length,rule=2,ties='ordered')$y
    tab$SR <- pmin(1,pmax(0,tab$SR))
  } else tab$SR <- NA_real_
  tab$Negative_residual <- pmax(0,tab$F_obs+tab$M_obs-tab$Total)
  tab$U_available <- if(usable) tab$U_obs else 0
  tab$F_allocated <- tab$U_available * ifelse(is.finite(tab$SR),tab$SR,0)
  tab$M_allocated <- tab$U_available - tab$F_allocated
  tab$F_corrected <- tab$F_obs + tab$F_allocated
  tab$M_corrected <- tab$M_obs + tab$M_allocated
  tab$Balance_error <- tab$F_corrected+tab$M_corrected-tab$Total
  to_lfd <- function(field) {
    out <- tab[,c('Year','Length',field),drop=FALSE]
    names(out)[3] <- 'Frequency'
    out[order(out$Year,out$Length),,drop=FALSE]
  }
  list(total=whole,F_observed=F,M_observed=M,U_observed=U,
       F_corrected=to_lfd('F_corrected'),M_corrected=to_lfd('M_corrected'),
       audit=tab,allocated=usable,
       message=if(usable) 'Undetermined allocated using pooled LOESS F/(F+M) curve' else
         'LOESS curve unavailable: display shows observed F/M only (no allocation)')
}

# Reproduces the former two-block plotLengths layout (year panels, blue mean,
# red L95), without the fragile FLQuant reshaping in sparse one-bin years.
medits_lfd_display_plot <- function(a,view='total',sex=c('C','F','M','U'),
                                    corrected=TRUE) {
  sex <- match.arg(sex)
  x <- medits_lfd_allocate(a,view)
  d <- switch(sex,C=x$total,F=if(corrected)x$F_corrected else x$F_observed,
              M=if(corrected)x$M_corrected else x$M_observed,U=x$U_observed)
  label <- switch(sex,C='Total',F='Female',M='Male',U='Undetermined')
  shiny::validate(shiny::need(is.data.frame(d) && nrow(d)>0 && any(d$Frequency>0),
                        paste('No positive',label,'LFD frequencies in selected view.')))
  years <- sort(unique(d$Year[d$Frequency>0]))
  d <- d[d$Year %in% years,,drop=FALSE]
  yrs_a <- years[seq_len(ceiling(length(years)/2))]
  yrs_b <- setdiff(years,yrs_a)
  classes <- sort(unique(d$Length))
  gap <- diff(classes); gap <- gap[is.finite(gap) & gap>0]
  barw <- if(length(gap)) min(gap)*0.88 else 1
  make_block <- function(yr) {
    sub <- d[d$Year %in% yr,,drop=FALSE]
    if(!nrow(sub)) return(NULL)
    sub$Year_panel <- factor(as.character(sub$Year),levels=as.character(yr))
    pieces <- lapply(split(sub,sub$Year),function(z){
      z <- z[z$Frequency>0,,drop=FALSE]
      if(!nrow(z)) return(NULL)
      z <- z[order(z$Length),,drop=FALSE]
      total <- sum(z$Frequency)
      avg <- sum(z$Length*z$Frequency)/total
      l95 <- if(nrow(z)==1L) z$Length[1] else
        stats::approx(cumsum(z$Frequency)/total,z$Length,xout=.95,
                      rule=2,ties='ordered')$y
      data.frame(Year_panel=as.character(z$Year[1]),mean=avg,l95=l95)
    })
    sums <- do.call(rbind,Filter(Negate(is.null),pieces))
    sums$Year_panel <- factor(sums$Year_panel,levels=as.character(yr))
    ggplot2::ggplot(sub,ggplot2::aes(x=Length,y=Frequency)) +
      ggplot2::geom_col(width=barw,fill='#3a7fa3',colour='#245775',linewidth=.12) +
      ggplot2::geom_vline(data=sums,ggplot2::aes(xintercept=mean),
                          inherit.aes=FALSE,colour='#2368d7',linewidth=.7) +
      ggplot2::geom_vline(data=sums,ggplot2::aes(xintercept=l95),
                          inherit.aes=FALSE,colour='#d83a44',linewidth=.7) +
      ggplot2::facet_wrap(~Year_panel,ncol=4,scales='free_y',drop=FALSE) +
      ggplot2::labs(x='Length (mm)',y='Standardised frequency (n/km²)') +
      ggplot2::theme_bw(base_size=11) +
      ggplot2::theme(strip.text=ggplot2::element_text(face='bold'),
                     panel.grid.minor=ggplot2::element_blank(),
                     panel.spacing=grid::unit(1.1,'lines'))
  }
  p1 <- make_block(yrs_a)
  if(length(yrs_b)) {
    p2 <- make_block(yrs_b)
    p <- ggpubr::ggarrange(p1,p2,ncol=1,nrow=2,heights=c(1,1))
  } else p <- p1
  p
}


# Prepare an existing per-haul map dataset, independently of the plotting layer.
# Only observed hauls with finite coordinates/values and positive swept area
# are returned. These are not the annual area-stratified MEDITS indices.
medits_spatial_prepare <- function(a,metric=c('BIOMASS','DENSITY'),view='total') {
  metric <- match.arg(metric)
  field <- if(metric=='BIOMASS') 'BIbubble_map' else 'DIbubble_map'
  d <- a[[field]]
  shiny::validate(shiny::need(is.data.frame(d) && nrow(d)>0,
                              paste('No spatial records available for',metric,
                                    '- rerun analysis with the updated functions file.')))
  required <- c('LONG','LAT','DEPTH','YEAR',metric)
  shiny::validate(shiny::need(all(required %in% names(d)),
                              paste('Missing spatial columns:',
                                    paste(setdiff(required,names(d)),collapse=', '))))
  d <- as.data.frame(d)
  d$LONG <- suppressWarnings(as.numeric(as.character(d$LONG)))
  d$LAT <- suppressWarnings(as.numeric(as.character(d$LAT)))
  d$DEPTH <- suppressWarnings(as.numeric(as.character(d$DEPTH)))
  d$YEAR <- suppressWarnings(as.integer(as.character(d$YEAR)))
  d[[metric]] <- suppressWarnings(as.numeric(as.character(d[[metric]])))
  ok <- is.finite(d$YEAR) & is.finite(d$LAT) & is.finite(d$LONG) &
        is.finite(d[[metric]]) & d[[metric]]>=0 &
        abs(d$LAT)<=90 & abs(d$LONG)<=180 & d$LAT!=0 & d$LONG!=0
  if(identical(view,'shelf')) ok <- ok & is.finite(d$DEPTH) & d$DEPTH>=1 & d$DEPTH<=200
  if(identical(view,'slope')) ok <- ok & is.finite(d$DEPTH) & d$DEPTH>200 & d$DEPTH<=800
  d <- d[ok,,drop=FALSE]
  shiny::validate(shiny::need(nrow(d)>0,
                              paste('No valid',metric,'map observations in the selected bathymetric domain.')))
  d$log_value <- log1p(d[[metric]])
  d
}

# Age–length keys from optional, individual age-reading CSV.
# These are sample-based probabilities P(Age | Length), NOT survey-raised LFDs.
medits_alk_prepare <- function(d,species=NULL,gsa=NULL,country=NULL) {
  if(!is.data.frame(d) || !nrow(d)) stop('Age readings CSV is empty.')
  key <- function(x) toupper(gsub('[^A-Z0-9]','',x))
  nm <- names(d)
  find <- function(candidates,required=FALSE) {
    i <- match(key(candidates),key(nm)); i <- i[!is.na(i)]
    if(length(i)) return(nm[i[1]])
    if(required) stop(paste0('Age CSV missing required column (one of): ',
                             paste(candidates,collapse=', ')))
    NULL
  }
  length_col <- find(c('LENGTH_CLASS','LENGTH_MM','LENGTH','LEN','LENGTH_CM'),TRUE)
  age_col <- find(c('AGE','AGE_YEARS','AGE_GROUP','AGE_CLASS'),TRUE)
  sex_col <- find(c('SEX','GENDER'))
  area_col <- find(c('AREA','GSA'),TRUE)
  country_col <- find('COUNTRY',TRUE)
  species_col <- find(c('SPECIES_CODE','SPECIES','3A_CODE'),TRUE)
  year_col <- find('YEAR',TRUE)
  count_col <- find(c('COUNT','N_AGED','N_READINGS'))
  filter <- rep(TRUE,nrow(d))
  if(!is.null(species_col) && length(species) && nzchar(species))
    filter <- filter & trimws(toupper(as.character(d[[species_col]]))) %in% species
  if(!is.null(area_col) && length(gsa))
    filter <- filter & trimws(as.character(d[[area_col]])) %in% gsa
  if(!is.null(country_col) && length(country))
    filter <- filter & trimws(toupper(as.character(d[[country_col]]))) %in% country
  d <- d[filter,,drop=FALSE]
  if(!nrow(d)) stop('No age readings match the selected species / GSA / country.')
  asnum <- function(z) suppressWarnings(as.numeric(gsub(',','.',as.character(z),fixed=TRUE)))
  yy <- if(is.null(year_col)) rep(NA_integer_,nrow(d)) else
    suppressWarnings(as.integer(as.character(d[[year_col]])))
  ll <- asnum(d[[length_col]])
  aa <- asnum(d[[age_col]])
  ss <- if(is.null(sex_col)) rep('U',nrow(d)) else toupper(trimws(as.character(d[[sex_col]])))
  n <- if(is.null(count_col)) rep(1,nrow(d)) else asnum(d[[count_col]])
  ok <- is.finite(ll) & ll>0 & is.finite(aa) & aa>=0 &
    abs(aa-round(aa))<1e-8 & is.finite(n) & n>0 & abs(n-round(n))<1e-8
  ans <- data.frame(Year=yy,Length=ll,Age=as.integer(round(aa)),
     Sex=ss,Count=n,stringsAsFactors=FALSE)
  ans <- ans[ok,,drop=FALSE]
  if(!nrow(ans)) stop('No valid age–length pairs: inspect AGE, LENGTH_CLASS and COUNT.')
  ans
}

medits_alk_probabilities <- function(d,sex='ALL',year='ALL') {
  if(!is.data.frame(d) || !nrow(d)) stop('No age readings available for ALK.')
  if(!identical(sex,'ALL')) d <- d[as.character(d$Sex)==as.character(sex),,drop=FALSE]
  if(!is.null(year) && !identical(year,'ALL'))
    d <- d[!is.na(d$Year) & as.character(d$Year)==as.character(year),,drop=FALSE]
  if(!nrow(d)) stop('No age readings for this sex and year selection.')
  tab <- stats::aggregate(Count~Length+Age,data=d,FUN=sum)
  names(tab)[3] <- 'N_aged'
  length_total <- stats::aggregate(N_aged~Length,data=tab,FUN=sum)
  names(length_total)[2] <- 'N_at_length'
  tab <- merge(tab,length_total,by='Length',sort=TRUE)
  tab$Probability <- tab$N_aged/tab$N_at_length
  tab$Group <- if(identical(sex,'ALL')) 'All sexes' else
    if(identical(sex,'F')) 'Female' else 'Male'
  tab[order(tab$Length,tab$Age),c('Group','Length','Age','N_aged','N_at_length','Probability')]
}

medits_alk_plot <- function(tab,sex='ALL') {
  shiny::validate(shiny::need(is.data.frame(tab) && nrow(tab)>0,
                              'No age–length probabilities to plot.'))
  ggplot2::ggplot(tab,ggplot2::aes(x=Length,y=Age,fill=Probability))+
    ggplot2::geom_tile(colour='white',linewidth=0.25)+
    ggplot2::scale_fill_gradientn(colours=c('#f3f9f8','#9ad8cb','#278b91','#093d5f'),
                                 limits=c(0,1),name='P(age | length)')+
    ggplot2::scale_y_continuous(breaks=sort(unique(tab$Age)))+
    ggplot2::labs(title=paste('Age–length key ·',unique(tab$Group)[1]),
       subtitle='Observed aged individuals; selected year, GSA and country filters applied',
       x='Length (original age-CSV units)',y='Age (years)')+
    ggplot2::theme_bw(base_size=12)+
    ggplot2::theme(plot.title=ggplot2::element_text(face='bold'),
                   panel.grid.minor=ggplot2::element_blank())
}