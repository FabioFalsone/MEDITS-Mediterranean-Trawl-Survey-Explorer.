# MEDITS report generator - self-contained, uses the packages already in the app.
# Figures are generated from the SAME reactive plot functions used on screen.
# The interpretations below read NUMERIC outputs, not pixels or a remote AI model.

medits_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
medits_fmt <- function(x, digits=2) {
  if (!length(x) || !is.finite(x[1])) return("NA")
  format(round(x[1], digits), big.mark=",", trim=TRUE, scientific=FALSE)
}
medits_years <- function(x) sort(unique(medits_num(x)[is.finite(medits_num(x))]))
medits_view <- function(d, view) {
  if (is.null(d) || !is.data.frame(d)) return(data.frame())
  if ("macro" %in% names(d)) d <- d[!is.na(d$macro) & d$macro == view, , drop=FALSE]
  d
}
medits_xy <- function(d, xc, yc) {
  if (!is.data.frame(d) || !all(c(xc,yc) %in% names(d))) return(data.frame(year=numeric(), value=numeric()))
  z <- data.frame(year=medits_num(d[[xc]]), value=medits_num(d[[yc]]))
  z <- z[is.finite(z$year) & is.finite(z$value), , drop=FALSE]
  z[order(z$year), , drop=FALSE]
}
medits_trend_text <- function(df, label, unit="", explain="") {
  if (!is.data.frame(df) || !nrow(df)) return(paste("No valid annual values available for", label, "."))
  df <- aggregate(value ~ year, data=df, FUN=mean)
  df <- df[order(df$year), , drop=FALSE]
  n <- nrow(df)
  if (n == 1L) return(paste0("One observed year (", df$year[1], "): ", medits_fmt(df$value[1]), " ", unit, ". ", explain))
  k <- min(5L, floor(n / 2))
  early <- mean(head(df$value,k))
  late <- mean(tail(df$value,k))
  delta <- late - early
  ratio <- if (abs(early) > 1e-10) 100*delta/abs(early) else NA_real_
  direction <- if (delta > 0) "higher" else if (delta < 0) "lower" else "equal"
  comparison <- if (is.finite(ratio)) paste0(" (", medits_fmt(abs(ratio), 1), "% ", if (delta>=0) "increase" else "decrease", ")") else " (percentage change undefined: early mean close to zero)"
  slope <- tryCatch(unname(coef(stats::lm(value ~ year, data=df))[2]), error=function(e) NA_real_)
  paste0("The first ",k," observed years have an average ",label," of ", medits_fmt(early)," ", unit,
         "; the last ",k," observed years average ", medits_fmt(late)," ",unit," (",direction,comparison,"). ",
         if (is.finite(slope)) paste0("Descriptive linear slope: ",medits_fmt(slope,3)," ",unit,"/year. ") else "",
         "Observed years: ",n," (",df$year[1],"-",df$year[n],"). ",explain)
}
medits_lfd_summary <- function(d, name) {
  if (!is.data.frame(d) || !all(c("Year","Length","Frequency") %in% names(d)) || !nrow(d)) return(data.frame())
  d <- data.frame(Year=medits_num(d$Year), Length=medits_num(d$Length), Frequency=medits_num(d$Frequency))
  d <- d[is.finite(d$Year) & is.finite(d$Length) & is.finite(d$Frequency) & d$Frequency >= 0, , drop=FALSE]
  if (!nrow(d)) return(data.frame())
  result <- lapply(split(d, d$Year), function(z) {
    count <- sum(z$Frequency,na.rm=TRUE)
    data.frame(Year=z$Year[1], Total_frequency=count,
      Weighted_mean_length_mm=if (count > 0) sum(z$Length*z$Frequency,na.rm=TRUE)/count else NA_real_,
      Modal_length_mm=if (count>0) z$Length[which.max(z$Frequency)] else NA_real_)
  })
  z <- do.call(rbind,result)
  z[order(z$Year), ,drop=FALSE]
}
medits_sr_year <- function(sr) {
  needed <- c("Year", "Frequency.x", "Frequency.y")
  if (!is.data.frame(sr) || !all(needed %in% names(sr)) || !nrow(sr)) return(data.frame())
  z <- data.frame(Year=medits_num(sr$Year), Female=medits_num(sr$Frequency.x), Male=medits_num(sr$Frequency.y))
  z <- z[is.finite(z$Year), , drop=FALSE]
  if (!nrow(z)) return(data.frame())
  z$Female[!is.finite(z$Female)] <- 0; z$Male[!is.finite(z$Male)] <- 0
  z <- stats::aggregate(cbind(Female,Male) ~ Year,data=z,FUN=sum)
  z$Female_proportion <- ifelse((z$Female+z$Male)>0,z$Female/(z$Female+z$Male), NA_real_)
  z[order(z$Year),,drop=FALSE]
}
medits_report_data <- function(a, view) {
  years <- medits_years(a$TAn$year)
  ta <- a$TAn
  # Stessi conteggi univoci di cale TA del grafico Survey period.
  ta_hauls <- medits_survey_unique_hauls(ta)
  stratum <- if (nrow(ta_hauls) && "strata" %in% names(ta_hauls)) {
    z <- data.frame(Year=medits_num(ta_hauls$year),
                    Stratum=as.character(ta_hauls$strata),Hauls=1)
    z <- z[is.finite(z$Year) & !is.na(z$Stratum),,drop=FALSE]
    if (nrow(z)) {
      z <- stats::aggregate(Hauls ~ Year+Stratum,data=z,FUN=sum)
      z <- reshape2::dcast(z,Year~Stratum,value.var="Hauls",fill=0)
      cols <- setdiff(names(z),"Year")
      z$Total_hauls <- rowSums(z[,cols,drop=FALSE],na.rm=TRUE)
      z
    } else data.frame()
  } else data.frame()
  survey <- medits_survey_period_summary(ta)
  sample <- data.frame()
  tc <- a$TCn
  if (is.data.frame(tc) && all(c("year","nblon") %in% names(tc))) {
    sample <- data.frame(Year=medits_num(tc$year),N=medits_num(tc$nblon),Sex=if("sex" %in% names(tc)) as.character(tc$sex) else "U")
    sample <- sample[is.finite(sample$Year) & is.finite(sample$N),,drop=FALSE]
    if(nrow(sample)) {
      sample$Sex <- toupper(trimws(sample$Sex))
      sample$Sex[!sample$Sex %in% c("F","M")] <- "U"
      sample$N_female <- ifelse(sample$Sex=="F",sample$N,0)
      sample$N_male <- ifelse(sample$Sex=="M",sample$N,0)
      sample$N_undetermined <- ifelse(sample$Sex=="U",sample$N,0)
      sample <- aggregate(cbind(N,N_female,N_male,N_undetermined)~Year, data=sample, FUN=sum)
      names(sample)[names(sample)=="N"] <- "N_total_raw"
      sample <- sample[order(sample$Year),,drop=FALSE]
    }
  }
  spatial_year <- function(metric) {
    d <- if (metric=="BIOMASS") a$BIbubble_map else a$DIbubble_map
    if (!is.data.frame(d) || !all(c("YEAR","DEPTH",metric) %in% names(d))) return(data.frame())
    z <- data.frame(Year=medits_num(d$YEAR),Depth=medits_num(d$DEPTH),Value=medits_num(d[[metric]]))
    z <- z[is.finite(z$Year) & is.finite(z$Value) & z$Value >= 0,,drop=FALSE]
    if (view=="shelf") z <- z[is.finite(z$Depth) & z$Depth>=1 & z$Depth<=200,,drop=FALSE]
    if (view=="slope") z <- z[is.finite(z$Depth) & z$Depth>200,,drop=FALSE]
    if (!nrow(z)) return(data.frame())
    x <- do.call(rbind,lapply(split(z,z$Year), function(t) data.frame(
      Year=t$Year[1],N_valid_hauls=nrow(t),N_positive=sum(t$Value>0),
      Pct_zero=100*mean(t$Value==0),Mean_per_haul=mean(t$Value),Max_per_haul=max(t$Value))))
    x[order(x$Year),,drop=FALSE]
  }
  biomass <- if (is.data.frame(a$BIOMASS)) a$BIOMASS else data.frame()
  density <- if (is.data.frame(a$DENSITY)) a$DENSITY else data.frame()
  # Guard against old cached results containing artificially inserted zero years.
  if (nrow(biomass) && "year" %in% names(biomass))
    biomass <- biomass[medits_num(biomass$year) %in% years,,drop=FALSE]
  if (nrow(density) && "year" %in% names(density))
    density <- density[medits_num(density$year) %in% years,,drop=FALSE]
  mw <- if (exists("medits_mean_weight_view",mode="function"))
    medits_mean_weight_view(a,view) else if(is.data.frame(a$Mweight)) a$Mweight else data.frame()
  ml <- medits_view(a$ML,view)
  sr <- medits_view(a$LFDSR,view)
  lfdt <- medits_view(a$LFD,view)
  lfdf <- medits_view(a$LFDF,view)
  lfdm <- medits_view(a$LFDM,view)
  lfd_t <- medits_lfd_summary(lfdt,"total")
  lfd_f <- medits_lfd_summary(lfdf,"female")
  lfd_m <- medits_lfd_summary(lfdm,"male")
  sr_y <- medits_sr_year(sr)
  ab_cols <- c("year","shelf_biomass","stdev_shelf","slope_biomass","stdev_slope","total_biomass","stdev_tot")
  de_cols <- c("year","shelf_density","stdev_shelf","slope_density","stdev_slope","total_density","stdev_tot")
  choose <- function(x,nm) {if (!is.data.frame(x)|| !nrow(x)) return(data.frame()); x[,intersect(nm,names(x)),drop=FALSE]}
  list(years=years,haul_counts=stratum,survey=survey,sample=sample,
    biomass=choose(biomass,ab_cols),density=choose(density,de_cols),
    mean_weight=mw,mean_length=ml,
    spatial_biomass=spatial_year("BIOMASS"), spatial_density=spatial_year("DENSITY"),
    lfd_total=lfdt, lfd_female=lfdf, lfd_male=lfdm,
    lfd_summary_total=lfd_t,lfd_summary_female=lfd_f,lfd_summary_male=lfd_m,
    sex_ratio_year=sr_y,sex_ratio_detail=sr)
}

medits_report_figure_catalog <- function(a, view, t, map_year, include_interpretation=TRUE) {
  yi <- paste0("GSA ",paste(a$report_gsas,collapse=", "))
  cap <- c(
    "Haul positions"="Georeferenced, valid unique TA haul positions in every observed survey year, faceted with fixed geographic limits. Each panel reports the number of mapped hauls in the selected depth domain; coastlines are outlines only.",
    "Survey period"="Annual MEDITS survey calendar. Each row represents one surveyed year; red ticks indicate distinct sampling dates and grey horizontal segments show the first-to-last sampling window. Numbers at right are valid unique TA hauls across all selected GSA(s) and depth strata. Coloured vertical bands indicate calendar quarters.",
    "Biomass spatial map"="Spatial distribution of haul-level biomass (kg/km2). Colours represent log(1 + biomass); zeros are retained.",
    "Density spatial map"="Spatial distribution of haul-level density (individuals/km2). Colours represent log(1 + density); zeros are retained.",
    "Hauls by stratum"="Number of survey hauls by year and MEDITS depth stratum (A-E) and macro-stratum, where available.",
    "Biomass trend"="Annual stratified biomass index (kg/km2) for the selected macro-stratum, with dispersion estimates.",
    "Density trend"="Annual stratified density index (individuals/km2) for the selected macro-stratum, with dispersion estimates.",
    "Mean weight"="Annual mean individual weight (g), calculated as the sum of TB catch weights divided by the sum of TB individuals for valid TA hauls in the selected depth domain (Total/Shelf/Slope). This is not the ratio of stratified biomass and density indices.",
    "Mean length"="Annual weighted mean length (mm) from the standardized length-frequency data; not separated by macro-stratum unless the source data contain a macro field.",
    "Total sample size"="Numbers of individuals measured per year using unraised TC records; across all sampled depth strata.",
    "Female sample size"="Numbers of females recorded in TC by year; raw measured counts, not standardized abundance.",
    "Male sample size"="Numbers of males recorded in TC by year; raw measured counts, not standardized abundance.",
    "LFD total"="Annual standardized length-frequency distributions; blue vertical lines mark mean lengths and red lines the 95th-percentile lengths. Interpret the frequency scale as derived by the app, not raw individuals.",
    "LFD females"="Annual standardized female length-frequency distributions, with mean and 95th-percentile markers.",
    "LFD males"="Annual standardized male length-frequency distributions, with mean and 95th-percentile markers.",
    "Sex ratio"="Length-specific female proportion F/(F + M) from standardized sex-specific frequency distributions; dotted horizontal line indicates 0.5."
  )
  desc <- setNames(rep("Numerical interpretation not available for this figure.",length(cap)), names(cap))
  desc["Haul positions"] <- if(is.data.frame(t$haul_counts) && "Year" %in% names(t$haul_counts)) {
    z <- t$haul_counts
    paste0("Survey effort spans ",nrow(z)," observed years (",min(z$Year),"-",max(z$Year),"). ",
           "Annual valid haul counts across ALL depth strata range from ",
           min(z$Total_hauls)," to ",max(z$Total_hauls),
           ". Figure labels refer to mapped hauls in the selected depth stratum, ",
           "so counts can be lower when a position is missing or the stratum is restricted.")
  } else "No valid annual haul counts available."
  desc["Survey period"] <- if(nrow(t$survey)) {
    z <- t$survey
    paste0("Sampled years: ",nrow(z)," (",min(z$Year),"-",max(z$Year),"). ",
           "Annual valid hauls range from ",min(z$Total_hauls)," to ",max(z$Total_hauls),
           "; overall total: ",sum(z$Total_hauls),
           ". Each red tick represents a distinct survey date, not a single haul. ",
           "Changes in survey timing or sampling effort may affect between-year comparisons.")
  } else "No valid survey haul records available."
  for (item in c("Biomass spatial map","Density spatial map")) {
    z <- if(item=="Biomass spatial map") t$spatial_biomass else t$spatial_density
    desc[item] <- if(nrow(z)) {last<-z[nrow(z),];paste0("In ",last$Year,", ",last$N_valid_hauls," valid haul-level measurements are plotted; ",medits_fmt(last$Pct_zero,1),"% are zero observations. The mean at sampled haul positions is ",medits_fmt(last$Mean_per_haul),"; this is not the stratified abundance index.")} else "No valid spatial observations for the selected view."
  }
  desc["Hauls by stratum"] <- if(nrow(t$haul_counts)) medits_trend_text(medits_xy(t$haul_counts,"Year","Total_hauls"),"survey hauls", "hauls", "Design and effort changes can affect year-to-year comparability.") else "No haul counts available."
  ycol <- if(view=="shelf") "shelf" else if(view=="slope") "slope" else "total"
  desc["Biomass trend"] <- medits_trend_text(medits_xy(t$biomass,"year",paste0(ycol,"_biomass")),"biomass", "kg/km2", "This is descriptive; neither statistical significance nor causal mechanisms are inferred.")
  desc["Density trend"] <- medits_trend_text(medits_xy(t$density,"year",paste0(ycol,"_density")),"density", "individuals/km2", "This is descriptive; neither statistical significance nor causal mechanisms are inferred.")
  if ("Mean_weight(g)" %in% names(t$mean_weight)) desc["Mean weight"] <- medits_trend_text(medits_xy(t$mean_weight,"Year","Mean_weight(g)"),"mean weight","g", "Calculated from TB records matched to valid TA hauls in the displayed depth domain.")
  if("Mean_Length_mm" %in% names(t$mean_length)) desc["Mean length"] <- medits_trend_text(medits_xy(t$mean_length,"Year","Mean_Length_mm"),"mean length","mm", "The current underlying calculation is not macro-stratum-specific.")
  if(nrow(t$sample)) {
    desc["Total sample size"] <- medits_trend_text(medits_xy(t$sample,"Year","N_total_raw"),"measured specimens","individuals", "Raw TC counts are not standardized density estimates.")
    desc["Female sample size"] <- medits_trend_text(medits_xy(t$sample,"Year","N_female"),"measured females","individuals", "Raw TC counts.")
    desc["Male sample size"] <- medits_trend_text(medits_xy(t$sample,"Year","N_male"),"measured males","individuals", "Raw TC counts.")
  }
  for (nm in c("LFD total","LFD females","LFD males")) {
    z <- switch(nm,"LFD total"=t$lfd_summary_total,"LFD females"=t$lfd_summary_female,"LFD males"=t$lfd_summary_male)
    desc[nm] <- if(nrow(z)) medits_trend_text(medits_xy(z,"Year","Weighted_mean_length_mm"),"weighted mean length","mm", "Based on standardized LFD frequencies; descriptive only.") else "No valid standardized frequencies in this group."
  }
  desc["Sex ratio"] <- if(nrow(t$sex_ratio_year)) {
    z<-t$sex_ratio_year;good<-z[is.finite(z$Female_proportion),,drop=FALSE]
    if(nrow(good)) paste0("Across ",nrow(good)," years with nonzero standardized sex-specific frequencies, yearly female proportions range from ",medits_fmt(min(good$Female_proportion))," to ",medits_fmt(max(good$Female_proportion)),". Ratios are derived from standardized LFD; length-specific curves can be unstable in sparsely sampled classes.") else "No defined yearly sex ratios."
  } else "No defined yearly sex ratios."
  lapply(names(cap),function(nm) list(title=nm,caption=unname(cap[[nm]]),interpretation=if(include_interpretation) unname(desc[[nm]]) else NULL))
}

# Pivot detailed long-form year/length data to compact matrices for the PDF.
# The original CSV / objects remain untouched. Missing bins remain NA, not zero.
medits_report_lfd_matrix <- function(d, value="Frequency") {
  if (!is.data.frame(d) || !all(c("Year","Length",value) %in% names(d)) || !nrow(d)) return(data.frame())
  z <- data.frame(Year=medits_num(d$Year),Length=medits_num(d$Length),Value=medits_num(d[[value]]))
  z <- z[is.finite(z$Year) & is.finite(z$Length) & is.finite(z$Value),,drop=FALSE]
  if (!nrow(z)) return(data.frame())
  z <- stats::aggregate(Value ~ Year + Length,data=z,FUN=sum)
  wide <- reshape2::dcast(z,Year ~ Length,value.var="Value",fill=NA_real_)
  if (ncol(wide)>2) {
    j <- order(suppressWarnings(as.numeric(names(wide)[-1])))+1L
    wide <- wide[,c(1L,j),drop=FALSE]
  }
  names(wide)[-1] <- paste0("L ",names(wide)[-1])
  wide[order(wide$Year),,drop=FALSE]
}

# Sex-ratio matrix reconstructed from female and male frequencies by length.
# 0/0 is undefined and stays NA rather than being represented as a true zero.
medits_report_sr_matrix <- function(d) {
  if (!is.data.frame(d) || !all(c("Year","Length","Frequency.x","Frequency.y") %in% names(d))) return(data.frame())
  z <- data.frame(Year=medits_num(d$Year),Length=medits_num(d$Length),
                  Fem=medits_num(d$Frequency.x),Male=medits_num(d$Frequency.y))
  z <- z[is.finite(z$Year) & is.finite(z$Length),,drop=FALSE]
  if (!nrow(z)) return(data.frame())
  z$Fem[!is.finite(z$Fem)] <- 0
  z$Male[!is.finite(z$Male)] <- 0
  z <- stats::aggregate(cbind(Fem,Male) ~ Year + Length,data=z,FUN=sum)
  z$SR <- ifelse(z$Fem+z$Male>0,z$Fem/(z$Fem+z$Male),NA_real_)
  medits_report_lfd_matrix(z,"SR")
}

medits_report_table_catalog <- function(t, detail=TRUE) {
  entries <- list(
    list("Survey timing and total hauls by year", t$survey, "Annual number of valid, unique TA hauls (all selected depth strata); distinct sampling days, first/last survey date. One haul is not one date."),
    list("Hauls by stratum and year", t$haul_counts, "Unique valid survey hauls in depth strata A-E, all-strata total. These are effort counts, not species catches."),
    list("Annual biomass indices", t$biomass, "Area-stratified indices in kg/km2 with the original variability estimates for shelf, slope, and total."),
    list("Annual density indices", t$density, "Area-stratified indices in individuals/km2 with the original variability estimates."),
    list("Mean weight by year", t$mean_weight, "Domain-specific individual mean weight (g): sum of TB catch weights divided by sum of TB individuals for valid TA hauls in the selected bathymetric domain."),
    list("Mean length by year", t$mean_length, "Mean length (mm), calculated from the source standardized LFD; not necessarily filtered by macro-stratum."),
    list("Measured specimens by year", t$sample, "Raw measured individuals from TC; female, male and undetermined totals."),
    list("Biomass map yearly summary", t$spatial_biomass, "Observed haul positions in the selected display stratum: zero fractions and per-haul values are NOT stratified indices."),
    list("Density map yearly summary", t$spatial_density, "Observed haul positions in the selected display stratum: zero fractions and per-haul values are NOT stratified indices."),
    list("LFD summary - all", t$lfd_summary_total, "Total standardized frequency and weighted mean/modal lengths by year."),
    list("LFD summary - females", t$lfd_summary_female, "Female standardized frequency and weighted mean/modal lengths by year."),
    list("LFD summary - males", t$lfd_summary_male, "Male standardized frequency and weighted mean/modal lengths by year."),
    list("Annual female proportion", t$sex_ratio_year, "Aggregated standardized sex-specific frequencies: F/(F+M); it is not a proportion from raw individuals."))
  if(detail) entries <- c(entries,list(
    list("Detailed LFD - all", medits_report_lfd_matrix(t$lfd_total), "Matrix of standardized frequencies by year and length class (mm). Blank/NA: no record for that cell."),
    list("Detailed LFD - females", medits_report_lfd_matrix(t$lfd_female), "Matrix of female standardized frequencies by year and length class (mm)."),
    list("Detailed LFD - males", medits_report_lfd_matrix(t$lfd_male), "Matrix of male standardized frequencies by year and length class (mm)."),
    list("Detailed sex ratio by length", medits_report_sr_matrix(t$sex_ratio_detail), "Year-by-length F/(F+M) from standardized sex-specific frequencies; NA where denominator is zero.")))
  out <- list()
  for(e in entries) {
    d <- e[[2]]
    if(!is.data.frame(d) || !nrow(d)) next
    # drop helper columns and non-atomic fields only; do not drop rows.
    d <- d[,vapply(d,function(x) is.atomic(x) && !is.list(x),logical(1)),drop=FALSE]
    if(!ncol(d)) next
    out[[length(out)+1L]] <- list(title=e[[1]],data=d,note=e[[3]])
  }
  out
}

medits_report_printable <- function(d) {
  d <- as.data.frame(d)
  for (nm in names(d)) {
    x <- d[[nm]]
    if (is.numeric(x)) {
      d[[nm]] <- ifelse(is.finite(x),
                       ifelse(x != 0 & abs(x) < 0.001,
                              formatC(x,format="e",digits=2),
                              format(round(x,3),trim=TRUE,scientific=FALSE)),
                       "NA")
    } else d[[nm]] <- ifelse(is.na(x),"NA",as.character(x))
    d[[nm]] <- substr(d[[nm]],1,20)
  }
  header <- gsub("_"," ",names(d))
  # Short readable headers, with metric units explained in each table note.
  alias <- c(Year="Year",year="Year",First_survey_day="First date",
             Last_survey_day="Last date",Days_sampled="Days",Total_hauls="Total hauls",
             shelf_biomass="Shelf biomass",slope_biomass="Slope biomass",
             total_biomass="Total biomass",shelf_density="Shelf density",
             slope_density="Slope density",total_density="Total density",
             stdev_shelf="SD shelf",stdev_slope="SD slope",stdev_tot="SD total",
             N_valid_hauls="N hauls",N_positive="N positive",Pct_zero="Zero (%)",
             Mean_per_haul="Mean/haul",Max_per_haul="Max/haul")
  jj <- which(names(d) %in% names(alias))
  header[jj] <- unname(alias[names(d)[jj]])
  names(d) <- vapply(header,function(x) substr(x,1,20),character(1))
  d
}

medits_report_table_pages <- function(d,max_rows=31L,max_cols=11L) {
  if(!is.data.frame(d)||!nrow(d)) return(list())
  # Keep 'Year' as the first column in every horizontal block of a wide matrix.
  # For annual (<=35 observations), one page usually suffices.
  n <- ncol(d)
  columns <- if(n<=max_cols) list(seq_len(n)) else {
    rest <- seq.int(2L,n)
    lapply(seq.int(1L,length(rest),by=max_cols-1L), function(i)
      c(1L, rest[i:min(length(rest),i+max_cols-2L)]))
  }
  pages <- list()
  for (cs in columns) for(i in seq(1L,nrow(d),by=max_rows)) {
    to <- min(i+max_rows-1L,nrow(d))
    pages[[length(pages)+1L]] <- list(data=d[i:to,cs,drop=FALSE],rows=c(i,to),columns=cs)
  }
  pages
}

medits_report_wrap <- function(s,width=110) {
  if(!length(s)||is.null(s)||is.na(s)||!nzchar(s)) return(character())
  strwrap(gsub("[[:space:]]+"," ",s),width=width)
}
medits_report_footer <- function(page) {
  grid::grid.text(paste0("MEDITS | Report page ",page), x=.955,y=.022,just="right",
                  gp=grid::gpar(fontsize=8,col="grey45"))
}
medits_report_paragraph <- function(text_string,y, size=10.5, line=.027, width=120) {
  lines <- medits_report_wrap(text_string,width)
  if(length(lines)>0) for(i in seq_along(lines))
    grid::grid.text(lines[i], x=.055,y=y-(i-1)*line,just=c("left","top"),
                    gp=grid::gpar(fontsize=size,col="grey25"))
  invisible(length(lines))
}
medits_report_index_pages <- function(items,lines_per_page=30L) {
  split(items,ceiling(seq_along(items)/lines_per_page))
}

medits_report_write_pdf <- function(file, a, view, plot_reactives, map_year=NULL,
                                    include_detail=TRUE,include_interpretation=TRUE) {
  if(!length(a$report_years)) stop("No years available for the report")
  if(is.null(map_year)||!is.finite(medits_num(map_year))) map_year <- max(a$report_years,na.rm=TRUE)
  tables_data <- medits_report_data(a,view)
  figs <- medits_report_figure_catalog(a,view,tables_data,map_year,include_interpretation)
  tbls <- medits_report_table_catalog(tables_data,detail=include_detail)
  available_figs <- names(plot_reactives)
  fig_keys <- vapply(figs,function(x) x$title,"")
  figs <- figs[fig_keys %in% available_figs]
  fig_keys <- vapply(figs,function(x) x$title,"")
  tbl_pgs <- lapply(tbls,function(x) medits_report_table_pages(x$data))
  entries <- c(paste0("Figure ",seq_along(figs),". ",fig_keys),
               if(length(tbls)) paste0("Table ",seq_along(tbls),". ",vapply(tbls,`[[`,"", "title")) else character())
  toc_pages <- max(1L,ceiling(length(entries)/27L))
  first_fig_page <- 2L + toc_pages
  fig_page <- if(length(figs)) first_fig_page+seq_along(figs)-1L else integer()
  first_table_page <- first_fig_page+length(figs)
  tbl_first <- integer(length(tbls))
  page_cursor <- first_table_page
  for(j in seq_along(tbls)) {
    tbl_first[j] <- page_cursor
    page_cursor <- page_cursor+length(tbl_pgs[[j]])
  }
  pages <- c(fig_page,tbl_first)
  # The index lists one entry per figure/table and first page of each table.
  grDevices::pdf(file, width=12,height=9,onefile=TRUE,paper="special",useDingbats=FALSE)
  on.exit(grDevices::dev.off(),add=TRUE)
  grid::grid.newpage()
  grid::grid.text("MEDITS - STANDARD ANALYSIS REPORT",y=.83,gp=grid::gpar(fontsize=25,fontface="bold"))
  grid::grid.text(paste0("Species: ",a$report_species),y=.69,gp=grid::gpar(fontsize=19))
  grid::grid.text(paste0("GSA: ",paste(a$report_gsas,collapse=", ")),y=.63,gp=grid::gpar(fontsize=19))
  grid::grid.text(paste0("Years: ",paste(range(a$report_years,na.rm=TRUE),collapse="-")),y=.57,gp=grid::gpar(fontsize=19))
  grid::grid.text(paste0("Survey: ",a$report_survey," | Display view: ",view),y=.50,gp=grid::gpar(fontsize=15))
  grid::grid.text(paste0("Generated: ",format(Sys.Date(),"%Y-%m-%d")),y=.42,gp=grid::gpar(fontsize=12))
  valid_years <- medits_years(a$report_years)
  missing_years <- setdiff(seq(min(valid_years),max(valid_years)), valid_years)
  if(length(missing_years)) grid::grid.text(
    paste0("No survey records for: ",paste(missing_years,collapse=", ")),
    y=.355,gp=grid::gpar(fontsize=11,col="grey35"))
  medits_report_paragraph("Notes: area-stratified biomass/density indices reflect the selected macro-stratum in figures; sample sizes and mean weight are calculated from all depth strata. Mean length and LFD are not macro-stratum-filtered unless their source objects include a 'macro' field. Spatial map averages are not stratified indices. Descriptive change is not a significance test or a causal interpretation.",.27,size=10,width=116,line=.030)
  medits_report_footer(1)
  for(k in seq_len(toc_pages)) {
    grid::grid.newpage()
    grid::grid.text(if(k==1L) "LIST OF FIGURES AND TABLES" else "LIST OF FIGURES AND TABLES (continued)",
                    x=.06,y=.94,just="left",gp=grid::gpar(fontsize=20,fontface="bold"))
    lo <- (k-1L)*27L+1L; hi <- min(k*27L,length(entries))
    if(lo <= hi) for (j in lo:hi) {
      iy <- j-lo
      line <- entries[j]
      # Limit long titles so the page number remains visible.
      grid::grid.text(substr(line,1,100),x=.07,y=.875-iy*.029,just="left",gp=grid::gpar(fontsize=11))
      grid::grid.text(as.character(pages[j]),x=.95,y=.875-iy*.029,just="right",gp=grid::gpar(fontsize=11))
    }
    medits_report_footer(1L+k)
  }
  issues <- character()
  page <- first_fig_page
  for(i in seq_along(figs)) {
    ff <- figs[[i]]; nm <- ff$title
    grid::grid.newpage()
    grid::grid.text(paste0("Figure ",i,". ",nm),x=.05,y=.966,just="left",
                    gp=grid::gpar(fontsize=17,fontface="bold"))
    p <- tryCatch(plot_reactives[[nm]](),error=function(e) {
      issues <<- c(issues,paste0(nm,": ",conditionMessage(e))); NULL
    })
    if(inherits(p,"ggplot") || inherits(p,"gg")) {
      tryCatch({
        # Reuse app plot objects, with adequate margins for caption + interpretation.
        pg <- ggplot2::ggplotGrob(p)
        grid::pushViewport(grid::viewport(x=.5,y=.585,width=.965,height=.725))
        grid::grid.draw(pg)
        grid::popViewport()
      }, error=function(e){
        issues <<- c(issues,paste0(nm,": ",conditionMessage(e)))
        grid::grid.text(paste0("Unable to render this figure: ",substr(conditionMessage(e),1,110)),y=.6)
      })
    } else {
      issues <- c(issues,paste0(nm,": plot did not return a ggplot"))
      grid::grid.text("This figure could not be generated.",y=.58,gp=grid::gpar(fontsize=13))
    }
    medits_report_paragraph(paste0("Caption. ",ff$caption),.225,size=10.5,width=125,line=.027)
    if(!is.null(ff$interpretation)) medits_report_paragraph(paste0("Automatic description. ",ff$interpretation),.135,size=9.2,width=143,line=.021)
    medits_report_footer(page); page <- page + 1L
  }
  if(length(tbls)) for(j in seq_along(tbls)) {
    tb <- tbls[[j]]
    tbpages <- tbl_pgs[[j]]
    for(ii in seq_along(tbpages)) {
      one <- tbpages[[ii]]
      grid::grid.newpage()
      grid::grid.text(paste0("Table ",j,". ",tb$title,if(length(tbpages)>1L) paste0(" (part ",ii,"/",length(tbpages),")") else ""),
                      x=.055,y=.947,just="left",gp=grid::gpar(fontsize=17,fontface="bold"))
      medits_report_paragraph(tb$note,.897,size=9.5,width=127,line=.026)
      grid::grid.text(paste0("Rows ",one$rows[1],"-",one$rows[2]," of ",nrow(tb$data),
                              "; columns ",paste(one$columns,collapse=", ")," of ",ncol(tb$data)),
                      x=.055,y=.805,just="left",gp=grid::gpar(fontsize=8,col="grey40"))
      dat <- medits_report_printable(one$data)
      nc <- ncol(dat); nr <- nrow(dat)
      fs <- if(nc >= 10) 7.4 else if(nc >= 8) 8.0 else 9.0
      grob <- gridExtra::tableGrob(dat,rows=NULL,theme=gridExtra::ttheme_minimal(
        base_size=fs,padding=grid::unit(c(.9,1.8),"mm"),
        core=list(bg_params=list(fill=rep(c("white","grey96"),length.out=nr),col=NA)),
        colhead=list(fg_params=list(fontface="bold",cex=.95),
                     bg_params=list(fill="grey85",col=NA))))
      # tableGrob uses absolute natural widths by default: explicit npc widths
      # are necessary to really occupy the printable page width.
      if(length(grob$widths)==nc) {
        weights <- pmax(5, pmin(22,nchar(names(dat))+2))
        if(nc>1) weights[1] <- max(6,min(weights[1],9))
        grob$widths <- grid::unit(weights/sum(weights),"npc")
      }
      # Uniform row heights, set relative to available table viewport.
      if(length(grob$heights)==nr+1L) {
        grob$heights <- grid::unit(rep(1/(nr+1),nr+1),"npc")
      }
      # Short tables use proportional height; 31-row tables fill most of the page.
      table_h <- min(.74,max(.18,(nr+1)*.0225))
      table_top <- .785
      grid::pushViewport(grid::viewport(x=.5,y=table_top-table_h/2,
                                        width=.94,height=table_h))
      grid::grid.draw(grob)
      grid::popViewport()
      medits_report_footer(page); page <- page + 1L
    }
  }
  if(length(issues)) {
    warning("Report pages with unavailable graphics: ",paste(issues,collapse=" | "))
  }
  invisible(list(figures=length(figs),tables=length(tbls),pages=page-1L,issues=issues))
}


# ============================================================================
# Survey period - VISTA SOLO REPORT PDF
# Un anno per riga, asse X calendario (gen-dic), date campionate e n. cale.
# Le date sono allineate a un anno di riferimento bisestile (2000): la stessa
# data giorno/mese occupa la stessa ascissa in anni bisestili e non bisestili.
# NON viene usata nella scheda Shiny (che mantiene i pannelli e l'altezza dinamica).
# ============================================================================
medits_make_survey_period_report_plot <- function(TAn, gsa="", state="") {
  summary <- medits_survey_period_summary(TAn)
  if (!nrow(summary)) stop("Survey period report: no valid TA hauls")
  summary <- summary[order(summary$Year), , drop=FALSE]
  n_years <- nrow(summary)
  summary$y <- seq.int(n_years, 1L)

  # Deduplicazione delle cale: esattamente la stessa logica del report e dell'app.
  # Conteggi di cale, non di giorni di campionamento.
  z <- medits_survey_unique_hauls(TAn)
  dates <- data.frame(Year=integer(), x=numeric(), y=numeric())
  if (nrow(z) && all(c("year", "month", "day") %in% names(z))) {
    yy <- suppressWarnings(as.integer(as.character(z$year)))
    mm <- suppressWarnings(as.integer(as.character(z$month)))
    dd <- suppressWarnings(as.integer(as.character(z$day)))
    dt <- suppressWarnings(as.Date(sprintf("%04d-%02d-%02d", yy, mm, dd),
                                   format="%Y-%m-%d"))
    cal <- suppressWarnings(as.Date(sprintf("2000-%02d-%02d", mm, dd),
                                    format="%Y-%m-%d"))
    day_x <- as.numeric(cal - as.Date("2000-01-01")) + 1
    valid <- !is.na(dt) & is.finite(day_x) & !is.na(yy)
    # Verifica sulle singole righe prima di deduplicare.
    dates <- unique(data.frame(Year=yy[valid], x=day_x[valid]))
    dates$y <- summary$y[match(dates$Year, summary$Year)]
    dates <- dates[is.finite(dates$x) & is.finite(dates$y), , drop=FALSE]
  }

  quarters <- data.frame(
    Quarter=factor(c("Q1 - winter", "Q2 - spring", "Q3 - summer", "Q4 - autumn"),
                   levels=c("Q1 - winter", "Q2 - spring", "Q3 - summer", "Q4 - autumn")),
    xmin=c(0.5,91.5,182.5,274.5),
    xmax=c(91.5,182.5,274.5,366.5),
    ymin=rep(0.5,4), ymax=rep(n_years+0.5,4)
  )
  month_first <- as.Date(sprintf("2000-%02d-01", seq_len(12)))
  month_pos <- as.numeric(month_first - as.Date("2000-01-01")) + 1

  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data=quarters,
      ggplot2::aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=Quarter),
      inherit.aes=FALSE,alpha=0.23,colour=NA) +
    ggplot2::geom_vline(xintercept=c(91.5,182.5,274.5),
                        colour="grey75",linewidth=0.27) +
    ggplot2::geom_vline(xintercept=370,colour="grey75",linewidth=0.3) +
    ggplot2::scale_fill_manual(
      name="Calendar quarter",
      values=c("Q1 - winter"="#8EC8E0", "Q2 - spring"="#9DDA99",
               "Q3 - summer"="#EED890", "Q4 - autumn"="#EBB29F"),drop=FALSE)

  if (nrow(dates)) {
    first <- stats::aggregate(x ~ Year, data=dates, FUN=min)
    names(first)[names(first)=="x"] <- "first_x"
    last <- stats::aggregate(x ~ Year, data=dates, FUN=max)
    names(last)[names(last)=="x"] <- "last_x"
    window <- merge(first,last,by="Year",all=TRUE)
    window$y <- summary$y[match(window$Year,summary$Year)]
    p <- p +
      ggplot2::geom_segment(
        data=window,
        ggplot2::aes(x=first_x,xend=last_x,y=y,yend=y),
        inherit.aes=FALSE,linewidth=1.10,colour="#7A7A7A",lineend="round") +
      ggplot2::geom_point(
        data=dates,ggplot2::aes(x=x,y=y),
        inherit.aes=FALSE,shape="|",size=2.6,colour="#AD2929")
  }

  p +
    ggplot2::geom_text(
      data=summary,
      ggplot2::aes(x=379,y=y,label=paste0("n = ",Total_hauls)),
      inherit.aes=FALSE,hjust=0,size=3.05,colour="#303030") +
    ggplot2::scale_x_continuous(
      breaks=month_pos,labels=month.abb,limits=c(0,423),
      expand=ggplot2::expansion(mult=0)) +
    ggplot2::scale_y_continuous(
      breaks=summary$y,labels=summary$Year,
      limits=c(0.5,n_years+0.5),
      expand=ggplot2::expansion(mult=0)) +
    ggplot2::labs(
      title=paste0("Survey period - GSA",gsa," ",state),
      subtitle="Annual survey calendar | red ticks = distinct survey dates; grey bars = sampling window",
      x="Calendar month",y=NULL) +
    ggplot2::theme_minimal(base_size=11) +
    ggplot2::theme(
      plot.title=ggplot2::element_text(face="bold",size=15,hjust=0.5),
      plot.subtitle=ggplot2::element_text(size=9,hjust=0.5,colour="grey25"),
      axis.text.y=ggplot2::element_text(size=9,colour="grey15"),
      axis.text.x=ggplot2::element_text(size=9,colour="grey25"),
      axis.title.x=ggplot2::element_text(size=10),
      panel.grid.minor=ggplot2::element_blank(),
      panel.grid.major.x=ggplot2::element_line(colour="white",linewidth=0.18),
      panel.grid.major.y=ggplot2::element_line(colour="grey86",linewidth=0.22),
      legend.position="bottom",
      legend.title=ggplot2::element_text(size=9),
      legend.text=ggplot2::element_text(size=8),
      legend.key.width=grid::unit(1.35,"lines"),
      plot.margin=ggplot2::margin(t=8,r=6,b=3,l=6,unit="pt"))
}