# medits_functions.R
library(data.table)
library(dplyr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(FLCore)
library(ggplotFL)
library(FLa4a)
library(ggpubr)
#library(ss3diags)
library(lattice)
library(mgcv) # per wireframe

mround <- function(x, base){
  base*round(x/base)
}

# Sum already-standardised LFD contributions by original MEDITS strata.
# A/B/C are shelf (10-200 m); D/E are slope (200-800 m).
# Do not re-raise TC or fabricate frequency from a pooled distribution.
medits_lfd_for_strata <- function(x, strata=c("A","B","C","D","E")) {
  empty <- data.frame(Year=integer(),Length=numeric(),Frequency=numeric())
  if (is.null(x) || !is.data.frame(x) || !nrow(x) ||
      !all(c("stratum","year","LC","lfd") %in% names(x))) return(empty)
  z <- x[as.character(x$stratum) %in% strata,,drop=FALSE]
  if(!nrow(z)) return(empty)
  yy <- suppressWarnings(as.integer(as.character(z$year)))
  ll <- suppressWarnings(as.numeric(as.character(z$LC)))
  ff <- suppressWarnings(as.numeric(as.character(z$lfd)))
  ok <- is.finite(yy) & is.finite(ll) & ll>0 & is.finite(ff)
  if(!any(ok)) return(empty)
  a <- stats::aggregate(ff[ok],by=list(Year=yy[ok],Length=ll[ok]),FUN=sum)
  names(a)[3] <- "Frequency"
  a[order(a$Year,a$Length),,drop=FALSE]
}

# Partitioning is additive across sex categories only when the source LFDs
# were raised on the same per-stratum scale. Flag substantial mismatches.
medits_lfd_undetermined <- function(total,female,male) {
  empty <- data.frame(Year=integer(),Length=numeric(),Frequency=numeric())
  if(!is.data.frame(total) || !nrow(total)) return(empty)
  a <- total[,c("Year","Length","Frequency"),drop=FALSE]
  names(a)[3] <- "All"
  if(is.data.frame(female) && nrow(female)) {
    b <- female[,c("Year","Length","Frequency"),drop=FALSE]
    names(b)[3] <- "Female"
    a <- merge(a,b,by=c("Year","Length"),all.x=TRUE)
  } else a$Female <- 0
  if(is.data.frame(male) && nrow(male)) {
    b <- male[,c("Year","Length","Frequency"),drop=FALSE]
    names(b)[3] <- "Male"
    a <- merge(a,b,by=c("Year","Length"),all.x=TRUE)
  } else a$Male <- 0
  a$Female[is.na(a$Female)] <- 0
  a$Male[is.na(a$Male)] <- 0
  residual <- a$All-a$Female-a$Male
  large_neg <- is.finite(residual) & residual < -1e-5 * pmax(1,abs(a$All))
  if(any(large_neg)) warning(sum(large_neg),
    " LFD cells: female + male exceeds total beyond numeric tolerance; inspect raw TC/raising")
  out <- data.frame(Year=a$Year,Length=a$Length,Frequency=pmax(0,residual))
  out[order(out$Year,out$Length),,drop=FALSE]
}

medits_lfd_sex_ratio <- function(female,male) {
  if(!is.data.frame(female) || !nrow(female) ||
     !is.data.frame(male) || !nrow(male)) return(data.frame())
  z <- merge(female,male,by=c("Year","Length"),all=TRUE,
             suffixes=c(".x",".y"))
  z$Frequency.x[is.na(z$Frequency.x)] <- 0
  z$Frequency.y[is.na(z$Frequency.y)] <- 0
  denom <- z$Frequency.x+z$Frequency.y
  z$SR <- ifelse(denom>0,z$Frequency.x/denom,NA_real_)
  z
}

safe_ggsave <- function(filename, plot, ...) {
  tryCatch(
    {
      ggsave(filename = filename, plot = plot, ...)
    },
    error = function(e) {
      message("safe_ggsave - errore nel salvataggio di ", filename, ": ", e$message)
      NULL
    }
  )
}

# Restored survey/map helpers from the last working MEDITS implementation.
# They use TAn (filtered valid TA hauls) and do not change biological estimates.
# MEDITS: coastline boundaries rendered as independent segments.
# No ggplot2 `group` aesthetic, avoiding polygon grouping/faceting mismatches.
medits_coast_segments <- function(xlims, ylims, database = "world") {
  empty <- data.frame(LONG = numeric(0), LAT = numeric(0),
                      LONG_END = numeric(0), LAT_END = numeric(0))
  coast <- tryCatch(
    maps::map(database = database, xlim = xlims, ylim = ylims,
              plot = FALSE, fill = FALSE),
    error = function(e) {
      warning("Coastline unavailable: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(coast) || length(coast$x) < 2L) return(empty)
  n <- length(coast$x)
  i <- seq_len(n - 1L)
  # NA separators mark the ends of independent islands / coastline pieces.
  ok <- is.finite(coast$x[i]) & is.finite(coast$y[i]) &
        is.finite(coast$x[i + 1L]) & is.finite(coast$y[i + 1L])
  i <- i[ok]
  data.frame(LONG = coast$x[i], LAT = coast$y[i],
             LONG_END = coast$x[i + 1L], LAT_END = coast$y[i + 1L])
}

# Survey period: conteggio delle CALE (TA valide), non dei giorni di campionamento.
# Una cala e' identificata da anno, Paese, GSA, nave e numero di cala.
# Calcolato sull'intera area batimetrica delle GSA/Paesi selezionati.
medits_survey_unique_hauls <- function(TAn) {
  if (!is.data.frame(TAn) || !nrow(TAn)) return(data.frame())
  if (!"year" %in% names(TAn)) stop("Survey period: manca la colonna 'year' in TAn")
  z <- as.data.frame(TAn)
  if ("validity" %in% names(z)) {
    z <- z[!is.na(z$validity) & as.character(z$validity) == "V", , drop = FALSE]
  }
  z$year <- suppressWarnings(as.integer(as.character(z$year)))
  z <- z[!is.na(z$year), , drop = FALSE]
  # Le righe di TA sono cale. Eliminiamo eventuali duplicati con la stessa chiave.
  key <- intersect(c("year", "country", "area", "vessel", "haul_number"), names(z))
  if (length(key) > 1L) {
    z <- z[!duplicated(z[, key, drop = FALSE]), , drop = FALSE]
  }
  z
}

medits_survey_period_summary <- function(TAn) {
  z <- medits_survey_unique_hauls(TAn)
  if (!nrow(z)) return(data.frame(Year=integer(), Total_hauls=integer(),
                                Days_sampled=integer(), First_survey_day=character(),
                                Last_survey_day=character()))
  yr <- sort(unique(z$year))
  out <- data.frame(Year=as.integer(yr),
                    Total_hauls=as.integer(tabulate(match(z$year, yr), nbins=length(yr))),
                    Days_sampled=0L, First_survey_day=NA_character_,
                    Last_survey_day=NA_character_)
  if (all(c("month", "day") %in% names(z))) {
    mm <- suppressWarnings(as.integer(as.character(z$month)))
    dd <- suppressWarnings(as.integer(as.character(z$day)))
    dates <- suppressWarnings(as.Date(sprintf("%04d-%02d-%02d", z$year, mm, dd),
                                      format="%Y-%m-%d"))
    for (k in seq_along(yr)) {
      present <- sort(unique(dates[z$year == yr[k] & !is.na(dates)]))
      out$Days_sampled[k] <- length(present)
      if (length(present)) {
        out$First_survey_day[k] <- format(min(present), "%Y-%m-%d")
        out$Last_survey_day[k] <- format(max(present), "%Y-%m-%d")
      }
    }
  }
  out
}

medits_make_survey_period_plot <- function(TAn, gsa="", state="", ncol=5L) {
  summary <- medits_survey_period_summary(TAn)
  if (!nrow(summary)) stop("Survey period: nessuna cala valida disponibile")
  z <- medits_survey_unique_hauls(TAn)

  # L'etichetta riporta il numero di CALE per anno, non il numero di date distinte.
  summary$Year_panel <- sprintf("%d\nn = %d hauls", summary$Year, summary$Total_hauls)
  summary$Year_panel <- factor(summary$Year_panel, levels=summary$Year_panel)
  year_label <- setNames(as.character(summary$Year_panel), as.character(summary$Year))
  plot_dates <- data.frame(
    Year=suppressWarnings(as.integer(as.character(z$year))),
    month=suppressWarnings(as.integer(as.character(z$month))),
    day=suppressWarnings(as.integer(as.character(z$day))))
  plot_dates$date <- suppressWarnings(as.Date(
    sprintf("%04d-%02d-%02d", plot_dates$Year, plot_dates$month, plot_dates$day),
    format="%Y-%m-%d"))
  plot_dates <- plot_dates[!is.na(plot_dates$date), , drop=FALSE]
  plot_dates <- unique(plot_dates[,c("Year","month","date"),drop=FALSE])
  plot_dates$dayofyear <- as.integer(format(plot_dates$date,"%j"))
  plot_dates$Year_panel <- factor(year_label[as.character(plot_dates$Year)],
                                  levels=levels(summary$Year_panel))
  quarters <- data.frame(Quarter=c("q1_winter","q2_spring","q3_summer","q4_fall"),
                         ymin=c(0,91,182,274), ymax=c(91,182,274,366))

  ggplot2::ggplot() +
    ggplot2::geom_rect(data=quarters,
      ggplot2::aes(xmin=0.5,xmax=12.5,ymin=ymin,ymax=ymax,fill=Quarter),
      inherit.aes=FALSE,alpha=0.34) +
    ggplot2::geom_point(data=plot_dates,
      ggplot2::aes(x=month,y=dayofyear), colour="firebrick",size=2.1) +
    # Forza la presenza del pannello anche se i giorni sono assenti o invalidi.
    ggplot2::geom_blank(data=summary,
      ggplot2::aes(x=6,y=180),inherit.aes=FALSE) +
    ggplot2::facet_wrap(~Year_panel,ncol=as.integer(ncol),drop=FALSE) +
    ggplot2::scale_x_continuous(breaks=c(1,3,6,9,12)) +
    ggplot2::scale_y_continuous(breaks=c(1,91,182,274,366)) +
    ggplot2::coord_cartesian(xlim=c(1,12),ylim=c(0,366),expand=FALSE) +
    ggplot2::scale_fill_manual("Season quarter",
      values=c(q1_winter="lightblue",q2_spring="lightgreen",
               q3_summer="lightyellow",q4_fall="lightsalmon")) +
    ggplot2::labs(title=paste0("Survey period - GSA",gsa," ",state),
                  subtitle="n = total valid survey hauls per year (all depth strata); red dots = distinct survey dates",
                  x="Month",y="Day of year") +
    ggplot2::theme_bw(base_size=12) +
    ggplot2::theme(legend.position="bottom",
                   strip.text=ggplot2::element_text(size=10,face="bold"),
                   plot.title=ggplot2::element_text(size=15,face="bold",hjust=.5),
                   plot.subtitle=ggplot2::element_text(size=9,hjust=.5),
                   axis.text=ggplot2::element_text(size=8))
}

run_medits_from_raw <- function(
    TAraw, TBraw, TCraw,
    areacode,       # vettore di GSA, es. c("16","17")
    countrycode,    # vettore di country, es. c("ITA","HRV")
    sspp,           # codice 3-alpha, es. "HKE"
    namesurvey  = "MEDITS",
    data_format = "DCF",      # solo etichetta
    global_sex_ratio = NULL,   # optional shared sex-ratio curve for GSA estimates
    output_base_dir = file.path(getwd(), paste0(namesurvey, "_output"))
){
  ## -----------------------------------------------------------------
  ## 1) Conversione TA/TB/TC in formato "medits_*"
  ## -----------------------------------------------------------------
  
  # Caso 1: input in formato DCF (colonne MAIUSCOLE, es. COUNTRY, AREA, ...)
  is_dcf_ta    <- all(c("COUNTRY","AREA","VESSEL","YEAR","HAUL_NUMBER") %in% names(TAraw))
  is_medits_ta <- all(c("country","area","vessel","year","haul_number") %in% names(TAraw))
  
  if (is_dcf_ta) {
    message("run_medits_from_raw: riconosciuto formato DCF (TAraw/TBraw/TCraw con colonne maiuscole). Eseguo la conversione in medits_*.")
    
    # ------ TA -> medits_ta ----------
    medits_ta <- data.frame(
      id      = NA_integer_,
      country = TAraw$COUNTRY,
      area    = TAraw$AREA,
      vessel  = TAraw$VESSEL,
      gear    = TAraw$GEAR,
      rigging = TAraw$RIGGING,
      door    = TAraw$DOORS,
      year    = TAraw$YEAR,
      month   = TAraw$MONTH,
      day     = TAraw$DAY,
      haul_number      = TAraw$HAUL_NUMBER,
      codend_closing   = TAraw$CODEND_CLOSING,
      shooting_time    = TAraw$SHOOTING_TIME,
      shooting_quadrant  = TAraw$SHOOTING_QUADRANT,
      shooting_latitude  = TAraw$SHOOTING_LATITUDE,
      shooting_longitude = TAraw$SHOOTING_LONGITUDE,
      shooting_depth     = TAraw$SHOOTING_DEPTH,
      hauling_time       = TAraw$HAULING_TIME,
      hauling_quadrant   = TAraw$HAULING_QUADRANT,
      hauling_latitude   = TAraw$HAULING_LATITUDE,
      hauling_longitude  = TAraw$HAULING_LONGITUDE,
      hauling_depth      = TAraw$HAULING_DEPTH,
      hauling_duration   = TAraw$HAUL_DURATION,
      validity           = TAraw$VALIDITY,
      course             = TAraw$COURSE,
      recorded_species   = TAraw$RECORDED_SPECIES,
      distance           = TAraw$DISTANCE,
      vertical_opening   = TAraw$VERTICAL_OPENING,
      wing_opening       = TAraw$WING_OPENING,
      geometrical_precision = TAraw$GEOMETRICAL_PRECISION,
      bridles_length     = TAraw$BRIDLES_LENGTH,
      warp_length        = TAraw$WARP_LENGTH,
      warp_diameter      = TAraw$WARP_DIAMETER,
      hydrological_station = TAraw$HYDROLOGICAL_STATION,
      observations       = TAraw$OBSERVATIONS,
      type_of_file               = TAraw$TYPE_OF_FILE,
      bottom_temperature_beginning = TAraw$BOTTOM_TEMPERATURE_BEGINNING,
      bottom_temperature_end       = TAraw$BOTTOM_TEMPERATURE_END,
      measuring_system_temp        = TAraw$MEASURING_SYSTEM,
      number_of_the_stratum        = TAraw$NUMBER_OF_THE_STRATUM,
      part_of_the_codend           = TAraw$PART_OF_THE_CODEND,
      upload_id          = NA_integer_,
      bottom_salinity_beginning = TAraw$BOTTOM_SALINITY_BEGINNING,
      bottom_salinity_end       = TAraw$BOTTOM_SALINITY_END,
      measuring_system_sal      = TAraw$MEASURING_SYSTEM_SALINITY
    )
    
    # ------ TB -> medits_tb ----------
    medits_tb <- data.frame(
      id      = NA_integer_,
      country = TBraw$COUNTRY,
      area    = TBraw$AREA,
      vessel  = TBraw$VESSEL,
      year    = TBraw$YEAR,
      haul_number    = TBraw$HAUL_NUMBER,
      codend_closing = TBraw$CODEND_CLOSING,
      partit         = TBraw$PART_OF_THE_CODEND,
      catfau         = TBraw$FAUNISTIC_CATEGORY,
      genus          = TBraw$GENUS,
      species        = TBraw$SPECIES,
      lref           = TBraw$NAME_OF_THE_REFERENCE_LIST,
      ptot           = TBraw$TOTAL_WEIGHT_IN_THE_HAUL,
      nbtot          = TBraw$TOTAL_NUMBER_IN_THE_HAUL,
      nbfem          = TBraw$NB_OF_FEMALES,
      nbmal          = TBraw$NB_OF_MALES,
      nbind          = TBraw$NB_OF_UNDETERMINED,
      tf             = NA_character_,
      month          = TBraw$MONTH,
      day            = TBraw$DAY,
      upload_id      = NA_integer_
    )
    
    # ------ TC -> medits_tc ----------
    medits_tc <- data.frame(
      id      = NA_integer_,
      country = TCraw$COUNTRY,
      area    = TCraw$AREA,
      vessel  = TCraw$VESSEL,
      year    = TCraw$YEAR,
      haul_number    = TCraw$HAUL_NUMBER,
      codend_closing = TCraw$CODEND_CLOSING,
      partit         = TCraw$PART_OF_THE_CODEND,
      genus          = TCraw$GENUS,
      species        = TCraw$SPECIES,
      codlon         = TCraw$LENGTH_CLASSES_CODE,
      pfrac          = TCraw$WEIGHT_OF_THE_FRACTION,
      pechan         = TCraw$WEIGHT_OF_THE_SAMPLE_MEASURED,
      sex            = TCraw$SEX,
      nbsex          = TCraw$NO_OF_INDIVIDUAL_OF_THE_ABOVE_SEX_MEASURED,
      length_class   = TCraw$LENGTH_CLASS,
      maturity       = TCraw$MATURITY,
      nblon          = TCraw$NUMBER_OF_INDIVIDUALS_IN_THE_LENGTH_CLASS_AND_MATURITY_STAGE,
      matsub         = TCraw$MATSUB,
      tf             = NA_character_,
      month          = TCraw$MONTH,
      day            = TCraw$DAY,
      catfau         = TCraw$FAUNISTIC_CATEGORY,
      upload_id      = NA_integer_
    )
    
    # Aggiungo comunque una colonna name_of_survey
    medits_ta$name_of_survey <- namesurvey
    medits_tb$name_of_survey <- namesurvey
    medits_tc$name_of_survey <- namesurvey
    
    # Caso 2: input Ã¨ giÃ  in formato medits_* (colonne minuscole)
  } else if (is_medits_ta) {
    message("run_medits_from_raw: input giÃ  in formato medits_ta/medits_tb/medits_tc. Salto la conversione DCF.")
    
    medits_ta <- as.data.frame(TAraw)
    medits_tb <- as.data.frame(TBraw)
    medits_tc <- as.data.frame(TCraw)
    
    # mi assicuro di avere le colonne 'id' e 'upload_id' come nella definizione originale
    if (!"id"        %in% names(medits_ta))  medits_ta$id        <- NA_integer_
    if (!"id"        %in% names(medits_tb))  medits_tb$id        <- NA_integer_
    if (!"id"        %in% names(medits_tc))  medits_tc$id        <- NA_integer_
    if (!"upload_id" %in% names(medits_ta))  medits_ta$upload_id <- NA_integer_
    if (!"upload_id" %in% names(medits_tb))  medits_tb$upload_id <- NA_integer_
    if (!"upload_id" %in% names(medits_tc))  medits_tc$upload_id <- NA_integer_
    
    # name_of_survey, se manca lo metto
    if (!"name_of_survey" %in% names(medits_ta))
      medits_ta$name_of_survey <- namesurvey
    if (!"name_of_survey" %in% names(medits_tb))
      medits_tb$name_of_survey <- namesurvey
    if (!"name_of_survey" %in% names(medits_tc))
      medits_tc$name_of_survey <- namesurvey
    
  } else {
    stop("run_medits_from_raw: formato di TAraw non riconosciuto. Colonne attese in maiuscolo (DCF) o minuscolo (medits).")
  }
  
  # Definisco gsa e state (usati piÃ¹ avanti e nel return)
  gsa   <- paste(areacode,    collapse = "_")
  state <- paste(countrycode, collapse = "_")
  
  ## -----------------------------------------------------------------
  ## 2) Parte specie / gen / spec / lfstep
  ## -----------------------------------------------------------------
  
  alpha_code <- sspp
  
  mdts  <- data.table::fread("Sp_Medits.csv")
  asfis <- data.table::fread("ASFIS_2020.csv")
  mdts  <- mdts[, c(2,3,7)]
  asfis <- asfis[, c(3,4)]
  colnames(mdts)[2] <- "Scientific_name"
  sp_list <- merge(mdts, asfis, by = "Scientific_name", all = TRUE)
  sp_list <- subset(sp_list[!is.na(sp_list$code),])
  sp_list <- subset(sp_list[!is.na(sp_list$`3A_CODE`),])
  sp_list$gen  <- stringr::str_sub(sp_list$code,  start = 1, end = 4)
  sp_list$spec <- stringr::str_sub(sp_list$code,  start = 5, end = 7)
  
  gen    <- unique(sp_list$gen[  sp_list$`3A_CODE` == sspp])
  spec   <- unique(sp_list$spec[ sp_list$`3A_CODE` == sspp])
  catfau <- unique(sp_list$catfau[sp_list$`3A_CODE` == sspp])
  
  lfstep <- unique(ifelse(catfau %in% c("B","Bst"), 1, 5))
  
  # max_SR come nel tuo codice originale
  if (sspp %in% c("ARA","ARS","DPS","HKE","MUT","MUR")){
    max_SR <- 1
  } else if (sspp %in% c("NEP")) {
    max_SR <- 0
  } else {
    max_SR <- 0.5
  }
  ## -----------------------------------------------------------------
  ## 3) Estrazione dati TA/TB/TC (SENZA NAME_OF_SURVEY)
  ## -----------------------------------------------------------------
  
  # TA
  TAn <- subset(
    medits_ta,
    area    %in% areacode &
      country %in% countrycode &
      validity == "V"
  )
  droplevels(TAn)
  
  # TB: filtro anche per specie selezionata
  TBn <- subset(
    medits_tb,
    area    %in% areacode &
      country %in% countrycode &
      genus   %in% gen &
      species %in% spec
  )
  droplevels(TBn)
  
  # TC: idem
  TCn <- subset(
    medits_tc,
    area    %in% areacode &
      country %in% countrycode &
      genus   %in% gen &
      species %in% spec
  )
  droplevels(TCn)
  
  ## -----------------------------------------------------------------
  ## 4) Elaborazioni successive
  ## -----------------------------------------------------------------
  
  # Correzioni specifiche NEP Italia
  TCn$length_class = ifelse(TCn$country=="ITA" & alpha_code=="NEP" & TCn$length_class>100,
                            TCn$length_class/10, TCn$length_class)
  TBn$ptot         = ifelse(TBn$country=="ITA" & alpha_code=="NEP" & TBn$ptot>10000,
                            TBn$ptot/1000, TBn$ptot)
  
  unique(TCn$year)
  
  # read in stratification table of the survey Copyrigth Tristan Rouyer
  stratification_scheme <- fread("MEDITS_Strata.csv")
  stratum<-stratification_scheme[stratification_scheme$AREA %in% areacode,]
  stratum<-stratum[stratum$COUNTRY %in% countrycode,]
  
  stratum$AREASTRATA[which(stratum$AREASTRATA==0)] <- as.integer(1)
  
  ## prepare the TA file for the next elaborations##
  meandepth=(TAn$shooting_depth+TAn$hauling_depth)/2
  sqkm=TAn$wing_opening/10000000*TAn$distance
  id2=paste(TAn$country,TAn$area,TAn$year,TAn$haul_number,sep="")
  TAn["strata"]=NA
  TAn=cbind(TAn,meandepth,sqkm,id2)
  if(length(grep('nstrate', names(TAn)))==0){
    TAn$nstrate=TAn$number_of_the_stratum
  }
  TAn$strata[]=TAn$nstrate[]
  
  ## Assigning strata code on meandepth value not on "codestrata"
  if (namesurvey=="MEDITS") {
    for (i in 1:length(TAn$strata))
      if(TAn$meandepth[i]>0 & TAn$meandepth[i] < 51){TAn$strata[i]="A"}else{
        if(TAn$meandepth[i]>=51 & TAn$meandepth[i] < 101){TAn$strata[i]="B"}else{
          if(TAn$meandepth[i]>=101 & TAn$meandepth[i] < 201){TAn$strata[i]="C"}else{
            if(TAn$meandepth[i]>=201 & TAn$meandepth[i] <501){TAn$strata[i]="D"}else{
              TAn$strata[i]="E"}}}}
    jones <- as.data.frame(table(TAn$strata))
    moss <- data.frame(Var1=c("A","B","C","D","E"),Freq1=as.integer(c(0,0,0,0,0)))
    cevert <-right_join(jones,moss)
    cevert <- cevert[order(as.character(cevert$Var1)),]
    cevert$FQ <- rowSums(cevert[,c("Freq", "Freq1")], na.rm=TRUE)
  }else{
    if (namesurvey=="SOLEMON"){
      for (i in 1:length(TAn$strata))
        if(TAn$meandepth[i]>0 & TAn$meandepth[i] < 30){TAn$strata[i]="A"}else{
          if(TAn$meandepth[i]>=30 & TAn$meandepth[i] < 50){TAn$strata[i]="B"}else{
            TAn$strata[i]="C"}}
      jones <- as.data.frame(table(TAn$strata))
      moss <- data.frame(Var1=c("A","B","C","D","E"),Freq1=as.integer(c(0,0,0,0,0)))
      cevert <- right_join(jones,moss)
      cevert <- cevert[order(as.character(cevert$Var1)),]
      cevert$FQ <- rowSums(cevert[,c("Freq", "Freq1")], na.rm=TRUE)
    }else{
      if (countrycode=="BGR"){
        for (i in 1:length(TAn$strata))
          if(TAn$meandepth[i]>15 & TAn$meandepth[i] < 50){TAn$strata[i]="A"}else{
            if(TAn$meandepth[i]>=50 & TAn$meandepth[i] < 75){TAn$strata[i]="B"}else{
              TAn$strata[i]="C"}}
        jones <- as.data.frame(table(TAn$strata))
        moss <- data.frame(Var1=c("A","B","C","D","E"),Freq1=as.integer(c(0,0,0,0,0)))
        cevert <- right_join(jones,moss)
        cevert <- cevert[order(as.character(cevert$Var1)),]
        cevert$FQ <- rowSums(cevert[,c("Freq", "Freq1")], na.rm=TRUE)
      }else{
        for (i in 1:length(TAn$strata))
          if(TAn$meandepth[i]>0 & TAn$meandepth[i] < 30){TAn$strata[i]="A"}else{
            if(TAn$meandepth[i]>=30 & TAn$meandepth[i] < 50){TAn$strata[i]="B"}else{
              if(TAn$meandepth[i]>=50 & TAn$meandepth[i] < 70){TAn$strata[i]="C"}}
            jones <- as.data.frame(table(TAn$strata))
            moss <- data.frame(Var1=c("A","B","C","D","E"),Freq1=as.integer(c(0,0,0,0,0)))
            cevert <- right_join(jones,moss)
            cevert <- cevert[order(as.character(cevert$Var1)),]
            cevert$FQ <- rowSums(cevert[,c("Freq", "Freq1")], na.rm=TRUE)
          }}}}
  unique(TAn$strata)
  table(TAn$strata)
  
  # prepare TB for next elaborations#
  id2=paste(TBn$country,TBn$area,TBn$year,TBn$haul_number,sep="")
  TBn=cbind(TBn,id2)
  # create a new database merging TATB
  TATBn=merge(TAn,TBn,by="id2",all=T)
  TATBn$W_sqkm=TATBn$ptot/TATBn$sqkm/1000
  TATBn$N_sqkm=TATBn$nbtot/TATBn$sqkm
  
  WD=getwd()
  
  # Use a stable output folder, not RStudio's currently active document.
  # The latter may point to another tab, be empty or change while running.
  dir_t <- file.path(output_base_dir, paste0("GSA", gsa))
  output_run_dir <- file.path(dir_t, state, alpha_code)
  dir.create(output_run_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(setwd(WD), add = TRUE)
  setwd(output_run_dir)
  write.table(TAn,file=paste("TA_",state,gsa,gen,spec,".csv"),sep=";",row.names=FALSE)
  write.table(TBn,file=paste("TB_",state,gsa,gen,spec,".csv"),sep=";",row.names=FALSE)
  write.table(TCn,file=paste("TC_",state,gsa,gen,spec,".csv"),sep=";",row.names=FALSE)
  write.table(TATBn,file=paste("TATB_",state,gsa,gen,spec,".csv"),sep=";",row.names=FALSE)
  
  # Checking hauls positions ####
  
  TA_lon=formatC(as.numeric(TAn$shooting_longitude),width=7,format='f',digits=2,flag='0')
  TA_gr_lon=substr(TA_lon, 1, 2)
  TA_gr_lon=as.integer(TA_gr_lon); TA_gr_lon[is.na(TA_gr_lon)] <- 0
  TA_mi_lon=substr(TA_lon, 3, 4)
  TA_mi_lon=as.integer(TA_mi_lon); TA_mi_lon[is.na(TA_mi_lon)] <- 0
  TA_se_lon=substr(TA_lon, 6, 7)
  TA_se_lon=as.integer(TA_se_lon); TA_se_lon[is.na(TA_se_lon)] <- 0
  
  TA_lon_fin=TA_gr_lon+(TA_mi_lon/60)+(TA_se_lon/6000)
  TAn$TA_lon_fin=TA_lon_fin
  
  TA_lat=formatC(as.numeric(TAn$shooting_latitude),width=7,format='f',digits=2,flag='0')
  TA_gr_lat=substr(TA_lat, 1, 2)
  TA_gr_lat=as.integer(TA_gr_lat); TA_gr_lat[is.na(TA_gr_lat)] <- 0
  TA_mi_lat=substr(TA_lat, 3, 4)
  TA_mi_lat=as.integer(TA_mi_lat); TA_mi_lat[is.na(TA_mi_lat)] <- 0
  TA_se_lat=substr(TA_lat, 6, 7)
  TA_se_lat=as.integer(TA_se_lat); TA_se_lat[is.na(TA_se_lat)] <- 0
  
  TA_lat_fin=TA_gr_lat+(TA_mi_lat/60)+(TA_se_lat/6000)
  TAn$TA_lat_fin=TA_lat_fin
  TAn$TA_lon_fin=TA_lon_fin
  
  TAn$TA_lon_fin=ifelse(TAn$shooting_quadrant=="7",-TAn$TA_lon_fin,TAn$TA_lon_fin)
  TA_lon_fin=TAn$TA_lon_fin
  
  xlims <- range(pretty(TA_lon_fin))
  ylims <- range(pretty(TA_lat_fin))
  
  # Exploring hauls time series ####
  hauls_table=aggregate(TAn$haul_number,list(TAn$year,TAn$strata),length)
  names(hauls_table)=c("Year","Stratum","Hauls")
  hauls_table2=hauls_table
  hauls_table_tot=aggregate(Hauls~Year, hauls_table, sum)
  hauls_table_tot$Stratum=rep("Total",nrow(hauls_table_tot))
  hauls_table_tot=hauls_table_tot[,c(1,3,2)]
  hauls_table=rbind(hauls_table,hauls_table_tot)
  
  Hauls_time_series=ggplot(hauls_table, aes(x=Year, y=Hauls,color=Stratum)) +
    geom_line(lwd=0.9) + facet_wrap(~Stratum,scales ="free") +
    theme(legend.position="none")+
    scale_x_continuous(breaks = seq(min(TAn$year),max(TAn$year),by=4))
  
  safe_ggsave(filename=paste(state,"_","GSA_",gsa,"Hauls_time_series.jpeg"),
         width = 10, height = 8, dpi = 150, units = "in",
         plot=Hauls_time_series)
  
  hauls_table3 <- dcast(hauls_table, Year ~ Stratum, sum)
  write.csv(hauls_table3,file=paste(state,"_","GSA_",gsa,"TotNb_Hauls_year.csv"),
            row.names = FALSE)
  
  # Exploring swept area ####
  exploring_sweptarea=ggplot(TAn, aes(x=meandepth, y=sqkm,color=strata)) +
    geom_point(cex=0.9) + facet_wrap(~strata,scales ="free") +
    theme(legend.position="none") + geom_smooth()
  safe_ggsave(filename=paste(state,"_","GSA_",gsa,"exploring_sweptarea.jpeg"),
         width = 10, height = 8, dpi = 150, units = "in",
         plot=exploring_sweptarea)
  
  # Exploring wing opening by meandepth area ####
  exploring_wingopening=ggplot(TAn, aes(x=meandepth, y=wing_opening/10,color=strata)) +
    geom_point(cex=0.9) + facet_wrap(~strata,scales ="free") +
    theme(legend.position="none") + geom_smooth()+
    ggtitle(paste0("GSA",gsa,"_",state)) +
    xlab("Mean_depth(m)") + ylab("Wing_Opening(m)")
  safe_ggsave(filename=paste(state,"_","GSA_",gsa,"exploring_wingopening.jpeg"),
         width = 10, height = 8, dpi = 150, units = "in",
         plot=exploring_wingopening)
  
  # Exploring survey period ####
  ggplot(TAn, aes(x=month, y=year))+
    geom_point(size=2.5,colour=2)+
    scale_y_continuous(breaks = seq(min(TAn$year),max(TAn$year),by=2))+
    theme(axis.text.x = element_text(angle = 0, hjust = 1,size=9))
  
  TA_date=TAn[,c('year', 'month','day')]
  TA_date$data=as.Date(paste0(TA_date$year,"-",TA_date$month,"-",TA_date$day))
  TA_date <- TA_date[order(TA_date$data),]
  TA_date$dayofyear <- as.numeric(format(TA_date$data, "%j"))
  TA_date=TA_date[!duplicated(TA_date$data), ]
  
  ggplot(TA_date, aes(x=month, y=dayofyear))+
    geom_point(size=2.5,colour=2)+
    facet_grid(~year)+
    scale_x_continuous(breaks = seq(min(TA_date$month),max(TA_date$month),by=1),
                       minor_breaks = FALSE)+
    theme(axis.text.x = element_text(angle = 0, hjust = 1,size=9))+
    ggtitle(paste0("GSA",gsa,state)) +
    xlab("Month") + ylab("Day_of_year")
  
  survey_period=ggplot(TA_date, aes(x=month, y=dayofyear))+
    geom_rect(data=NULL,aes(xmin=min(month),xmax=max(month),ymin=0,ymax=90,
                            fill="q1_winter"))+
    geom_rect(data=NULL,aes(xmin=min(month),xmax=max(month),ymin=91,ymax=180,
                            fill="q2_spring"))+
    geom_rect(data=NULL,aes(xmin=min(month),xmax=max(month),ymin=181,ymax=270,
                            fill="q3_summer"))+
    geom_rect(data=NULL,aes(xmin=min(month),xmax=max(month),ymin=271,ymax=365,
                            fill="q4_fall"))+
    scale_fill_manual('Season_quarter',
                      values = c("lightblue","lightgreen","lightyellow","lightsalmon"))+
    geom_point(size=3,colour="red")+
    facet_wrap(~year)+
    ggtitle(paste0("GSA",gsa,state)) +
    xlab("Month") + ylab("Day_of_year")
  
  safe_ggsave(filename=paste(state,"_","GSA_",gsa,"survey_period.jpeg"),
         width = 10, height = 8, dpi = 150, units = "in",
         plot=survey_period)
  
  # Species occurence ########
  tothaulyear=tapply(TATBn$haul_number.x,TATBn$year.x,length)
  TATBn['pos']=1
  TATBn['neg']=0
  TATBn['occurence']=NA
  TATBn$nbtot[is.na(TATBn$nbtot)]=0
  occur=ifelse((TATBn$nbtot>0),
               (TATBn$occurence=TATBn$pos),
               (TATBn$occurence=TATBn$neg))
  TATBn=cbind(TATBn,occur)
  tothaulposyear=tapply(TATBn$occur,TATBn$year.x,sum)
  occurence=tothaulposyear/tothaulyear
  year=sort(unique(TAn$year))
  occurence=as.data.frame(cbind(year,occurence))
  names(occurence)=c("year","occurence")
  write.csv(occurence,file=paste(gen,spec,"_GSA_",gsa,state,"_",
                                 "Occurence.csv",sep=""),row.names=FALSE)
  
  # Abundance Indexes by strata####
  strata=aggregate(TAn$sqkm,list(TAn$year,TAn$strata),sum)
  names(strata)=c("year","stratum","sweptarea")
  weigthbystrata=aggregate(TATBn$ptot,list(TATBn$year.y,TATBn$strata),sum)
  names(weigthbystrata)=c("year","stratum","weigth")
  numberbystrata=aggregate(TATBn$nbtot,list(TATBn$year.y,TATBn$strata),sum)
  names(numberbystrata)=c("year","stratum","number")
  biomass=merge(strata,weigthbystrata,by=c("year","stratum"))
  biomass$kg_km2=((biomass$weigth/biomass$sweptarea)/1000)
  density=merge(strata,numberbystrata,by=c("year","stratum"))
  density$n_km2=((density$number/density$sweptarea))
  
  # Area of each stratum####
  area.shelf1=aggregate(stratum$AREASTRATA,
                        by=list(stratum$CODESTRATA=="A"|
                                  stratum$CODESTRATA=="B"|
                                  stratum$CODESTRATA=="C"),FUN=sum)
  area.str.A=aggregate(stratum$AREASTRATA,by=list(stratum$CODESTRATA=="A"),FUN=sum)
  area.str.A1=area.str.A[2,2]
  area.str.B=aggregate(stratum$AREASTRATA,by=list(stratum$CODESTRATA=="B"),FUN=sum)
  area.str.B1=area.str.B[2,2]
  area.str.C=aggregate(stratum$AREASTRATA,by=list(stratum$CODESTRATA=="C"),FUN=sum)
  area.str.C1=area.str.C[2,2]
  area.str.D=aggregate(stratum$AREASTRATA,by=list(stratum$CODESTRATA=="D"),FUN=sum)
  area.str.D1=area.str.D[2,2]
  area.str.E=aggregate(stratum$AREASTRATA,by=list(stratum$CODESTRATA=="E"),FUN=sum)
  area.str.E1=area.str.E[2,2]
  area.shelf=area.shelf1[2,2]
  area.slope=area.shelf1[1,2]
  area.tot=area.shelf+area.slope
  
  # Weigth coefficients by stratum####
  W1sh=area.str.A1/area.shelf
  W2sh=area.str.B1/area.shelf
  W3sh=area.str.C1/area.shelf
  W4sl=area.str.D1/area.slope
  W5sl=area.str.E1/area.slope
  W1tot=area.str.A1/area.tot
  W2tot=area.str.B1/area.tot
  W3tot=area.str.C1/area.tot
  W4tot=area.str.D1/area.tot
  W5tot=area.str.E1/area.tot
  
  # Indexes by macrostrata (shelf, slope and total)####
  # BIOMASS ####
  #SHELF#
  BIA=subset(biomass,stratum=="A")
  BIA$W1sh=rep(W1sh,nrow(BIA))
  BIA$biomA=BIA$kg_km2*BIA$W1sh
  BIB=subset(biomass,stratum=="B")
  BIB$W2sh=rep(W2sh,nrow(BIB))
  BIB$biomB=BIB$kg_km2*BIB$W2sh
  BIC=subset(biomass,stratum=="C")
  BIC$W3sh=rep(W3sh,nrow(BIC))
  BIC$biomC=BIC$kg_km2*BIC$W3sh
  BIshelf=merge(BIA,BIB,by="year",all=TRUE)
  BIshelf=merge(BIshelf,BIC,by="year",all=TRUE)
  BIshelf$shelf_biomass <- rowSums(BIshelf[c(7,13,19)], na.rm=TRUE)
  #SLOPE#
  BID=subset(biomass,stratum=="D")
  BID$W4sl=rep(W4sl,nrow(BID))
  BID$biomD=BID$kg_km2*BID$W4sl
  BIE=subset(biomass,stratum=="E")
  BIE$W5sl=rep(W5sl,nrow(BIE))
  BIE$biomE=BIE$kg_km2*BIE$W5sl
  BIslope=merge(BID,BIE,by="year",all=TRUE)
  BIslope$slope_biomass <- rowSums(BIslope[c(7,13)], na.rm=TRUE)
  #Total#
  BIA=subset(biomass,stratum=="A")
  BIA$W1tot=rep(W1tot,nrow(BIA))
  BIA$biomA=BIA$kg_km2*BIA$W1tot
  BIB=subset(biomass,stratum=="B")
  BIB$W2tot=rep(W2tot,nrow(BIB))
  BIB$biomB=BIB$kg_km2*BIB$W2tot
  BIC=subset(biomass,stratum=="C")
  BIC$W3tot=rep(W3tot,nrow(BIC))
  BIC$biomC=BIC$kg_km2*BIC$W3tot
  BID=subset(biomass,stratum=="D")
  BID$W4tot=rep(W4tot,nrow(BID))
  BID$biomD=BID$kg_km2*BID$W4tot
  BIE=subset(biomass,stratum=="E")
  BIE$W5tot=rep(W5tot,nrow(BIE))
  BIE$biomE=BIE$kg_km2*BIE$W5tot
  BITot=merge(BIA,BIB,by="year",all=TRUE)
  BITot=merge(BITot,BIC,by="year",all=TRUE)
  BITot=merge(BITot,BID,by="year",all=TRUE)
  BITot=merge(BITot,BIE,by="year",all=TRUE)
  BITot$total_biomass <- rowSums(BITot[c(7,13,19,25,31)], na.rm=TRUE)
  
  BIOMASS=merge(BIshelf,BIslope,by="year",all=TRUE)
  BIOMASS=merge(BIOMASS,BITot,by="year",all=TRUE)
  BIOMASS=BIOMASS[,-c(2:19,21:32,34:63)]
  
  #plotting some graphs
  BMS=reshape2::melt(BIOMASS,id="year")
  safe_ggsave(filename=paste(state,"_","GSA_",gsa,"BIOMASS.jpeg"),
         width = 10, height = 8, dpi = 150, units = "in",
         plot=ggplot(BMS,aes(x=year,y=value,color=variable))+
           geom_line()+
           facet_grid(~variable)+
           ggtitle(paste0(gen,spec,"_","GSA",gsa,"_",state)) +
           xlab("Year") + ylab("Biomass (kg/km2)")+
           theme(axis.text.x = element_text(angle=65, vjust=0.6))+
           theme(legend.position="none"))
  
  # DENSITY ####
  #SHELF#
  DIA=subset(density,stratum=="A")
  DIA$W1sh=rep(W1sh,nrow(DIA))
  DIA$biomA=DIA$n_km2*DIA$W1sh
  DIB=subset(density,stratum=="B")
  DIB$W2sh=rep(W2sh,nrow(DIB))
  DIB$biomB=DIB$n_km2*DIB$W2sh
  DIC=subset(density,stratum=="C")
  DIC$W3sh=rep(W3sh,nrow(DIC))
  DIC$biomC=DIC$n_km2*DIC$W3sh
  DIshelf=merge(DIA,DIB,by="year",all=TRUE)
  DIshelf=merge(DIshelf,DIC,by="year",all=TRUE)
  DIshelf$shelf_density <- rowSums(DIshelf[c(7,13,19)], na.rm=TRUE)
  #SLOPE#
  DID=subset(density,stratum=="D")
  DID$W4sl=rep(W4sl,nrow(DID))
  DID$biomD=DID$n_km2*DID$W4sl
  DIE=subset(density,stratum=="E")
  DIE$W5sl=rep(W5sl,nrow(DIE))
  DIE$biomE=DIE$n_km2*DIE$W5sl
  DIslope=merge(DID,DIE,by="year",all=TRUE)
  DIslope$slope_density <- rowSums(DIslope[c(7,13)], na.rm=TRUE)
  #Total#
  DIA=subset(density,stratum=="A")
  DIA$W1tot=rep(W1tot,nrow(DIA))
  DIA$biomA=DIA$n_km2*DIA$W1tot
  DIB=subset(density,stratum=="B")
  DIB$W2tot=rep(W2tot,nrow(DIB))
  DIB$biomB=DIB$n_km2*DIB$W2tot
  DIC=subset(density,stratum=="C")
  DIC$W3tot=rep(W3tot,nrow(DIC))
  DIC$biomC=DIC$n_km2*DIC$W3tot
  DID=subset(density,stratum=="D")
  DID$W4tot=rep(W4tot,nrow(DID))
  DID$biomD=DID$n_km2*DID$W4tot
  DIE=subset(density,stratum=="E")
  DIE$W5tot=rep(W5tot,nrow(DIE))
  DIE$biomE=DIE$n_km2*DIE$W5tot
  DITot=merge(DIA,DIB,by="year",all=TRUE)
  DITot=merge(DITot,DIC,by="year",all=TRUE)
  DITot=merge(DITot,DID,by="year",all=TRUE)
  DITot=merge(DITot,DIE,by="year",all=TRUE)
  DITot$total_density <- rowSums(DITot[c(7,13,19,25,31)], na.rm=TRUE)
  
  DENSITY=merge(DIshelf,DIslope,by="year",all=TRUE)
  DENSITY=merge(DENSITY,DITot,by="year",all=TRUE)
  DENSITY=DENSITY[,-c(2:19,21:32,34:63)]
  
  DMS=reshape2::melt(DENSITY,id="year")
  safe_ggsave(filename=paste(state,"_","GSA_",gsa,"DENSITY.jpeg"),
         width = 10, height = 8, dpi = 150, units = "in",
         plot=ggplot(DMS,aes(x=year,y=value,color=variable))+
           geom_line()+
           facet_grid(~variable)+
           ggtitle(paste0(gen,spec,"_","GSA",gsa,"_",state)) +
           xlab("Year") + ylab("DENSITY (n/km2)")+
           theme(axis.text.x = element_text(angle=65, vjust=0.6))+
           theme(legend.position="none"))
  
  # VARIABILITY   ####
  # Hauls by year ####
  hauls_table2A=subset(hauls_table,Stratum=="A")
  hauls_table2B=subset(hauls_table,Stratum=="B")
  hauls_table2C=subset(hauls_table,Stratum=="C")
  hauls_table2D=subset(hauls_table,Stratum=="D")
  hauls_table2E=subset(hauls_table,Stratum=="E")
  
  # Swept area by stratum ####
  if(cevert$FQ[1]>0){
    swA=subset(strata,stratum=="A")
  }else{swA=data.frame(year=as.integer(year),stratum="A",
                       sweptarea=as.numeric(0.001))
  }
  if(cevert$FQ[2]>0){
    swB=subset(strata,stratum=="B")
  }else{swB=data.frame(year=as.integer(year),stratum="B",
                       sweptarea=as.numeric(0.001))
  }
  if(cevert$FQ[3]>0){
    swC=subset(strata,stratum=="C")
  }else{swC=data.frame(year=as.integer(year),stratum="C",
                       sweptarea=as.numeric(0.001))
  }
  if(cevert$FQ[4]>0){
    swD=subset(strata,stratum=="D")
  }else{swD=data.frame(year=as.integer(year),stratum="D",
                       sweptarea=as.numeric(0.001))
  }
  if(cevert$FQ[5]>0){
    swE=subset(strata,stratum=="E")
  }else{swE=data.frame(year=as.integer(year),stratum="E",
                       sweptarea=as.numeric(0.001))
  }
  
  # Total area ####
  area.str.A<- sum(stratum[stratum$CODESTRATA %in% "A",7])
  area.str.B<- sum(stratum[stratum$CODESTRATA %in% "B",7])
  area.str.C<- sum(stratum[stratum$CODESTRATA %in% "C",7])
  area.str.D<- sum(stratum[stratum$CODESTRATA %in% "D",7])
  area.str.E<- sum(stratum[stratum$CODESTRATA %in% "E",7])
  
  # Check GSA surface #### 
  tot_shelf=sum(stratum[stratum$CODESTRATA %in% c("A","B","C"),7])
  tot_slope=sum(stratum[stratum$CODESTRATA %in% c("D","E"),7])
  tot_area=sum(stratum[stratum$CODESTRATA %in% c("A","B","C","D","E"),7])
  
  # Correction factor f ####
  swA$f=swA$sweptarea/as.integer(area.str.A)
  swB$f=swB$sweptarea/as.integer(area.str.B)
  swC$f=swC$sweptarea/as.integer(area.str.C)
  swD$f=swD$sweptarea/as.integer(area.str.D)
  swE$f=swE$sweptarea/as.integer(area.str.E)
  SW=rbind(swA,swB,swC,swD,swE)
  SW$id3=paste0(SW$year,SW$stratum,sep="")
  
  # Elaboration ####
  TATBn$id3=paste0(TATBn$year.x,TATBn$strata,sep="")
  noTAwithTB <- TATBn[is.na(TATBn$strata),]
  write.csv(noTAwithTB,file="TBhaul_no_in_TAhaul.csv",row.names = FALSE)
  TATBn <- TATBn[complete.cases(TATBn[ , 12]),]
  biomass$id3=paste0(biomass$year,biomass$stratum,sep="")
  density$id3=paste0(density$year,density$stratum,sep="")
  hauls_table2$id3=paste0(hauls_table2$Year,hauls_table2$Stratum,sep="")
  TATBBI=merge(TATBn,hauls_table2,by="id3",all=TRUE)
  TATBBI=merge(TATBBI,biomass,by="id3",all.y=TRUE)
  TATBDI=merge(TATBn,hauls_table2,by="id3",all=TRUE)
  TATBDI=merge(TATBDI,density,by="id3",all.y=TRUE)
  
  # BIOMASS ####
  TATBBI$var1=TATBBI$sqkm*(TATBBI$W_sqkm-TATBBI$kg_km2)^2
  var1=aggregate(TATBBI$var1,list(TATBBI$year,TATBBI$strata),sum,na.rm=TRUE)
  var1$id3=paste0(var1$Group.1,var1$Group.2,sep="")
  
  var2a=merge(var1,hauls_table2,by="id3",all=TRUE)
  var2a$variance=(1/(var2a$Hauls-1))*var2a$x
  var2a$devst=sqrt(var2a$variance)
  var2a=merge(var2a,SW,by="id3",all=TRUE)
  
  var2a$Stratum <- var2a$stratum
  var2a$Group.2 <- var2a$stratum
  var2a$Group.1 <- var2a$year
  var2a$Year <- var2a$year
  var2a[is.na(var2a)] <- 0
  
  VarIA=subset(var2a,var2a$Stratum=="A")
  VarIB=subset(var2a,var2a$Stratum=="B")
  VarIC=subset(var2a,var2a$Stratum=="C")
  VarID=subset(var2a,var2a$Stratum=="D")
  VarIE=subset(var2a,var2a$Stratum=="E")
  
  VarIA$W1sh=W1sh
  VarIB$W2sh=W2sh
  VarIC$W3sh=W3sh
  VarID$W4sl=W4sl
  VarIE$W5sl=W5sl
  VarIA$W1tot=W1tot
  VarIB$W2tot=W2tot
  VarIC$W3tot=W3tot
  VarID$W4tot=W4tot
  VarIE$W5tot=W5tot
  
  VarIA$Vsh=(VarIA$variance*VarIA$W1sh^2)/(VarIA$sweptarea*(1-VarIA$f))
  VarIB$Vsh=(VarIB$variance*VarIB$W2sh^2)/(VarIB$sweptarea*(1-VarIB$f))
  VarIC$Vsh=(VarIC$variance*VarIC$W3sh^2)/(VarIC$sweptarea*(1-VarIC$f))
  VarID$Vsl=(VarID$variance*VarID$W4sl^2)/(VarID$sweptarea*(1-VarID$f))
  VarIE$Vsl=(VarIE$variance*VarIE$W5sl^2)/(VarIE$sweptarea*(1-VarIE$f))
  
  VarIA$Vt=(VarIA$variance*VarIA$W1tot^2)/(VarIA$sweptarea*(1-VarIA$f))
  VarIB$Vt=(VarIB$variance*VarIB$W2tot^2)/(VarIB$sweptarea*(1-VarIB$f))
  VarIC$Vt=(VarIC$variance*VarIC$W3tot^2)/(VarIC$sweptarea*(1-VarIC$f))
  VarID$Vt=(VarID$variance*VarID$W4tot^2)/(VarID$sweptarea*(1-VarID$f))
  VarIE$Vt=(VarIE$variance*VarIE$W5tot^2)/(VarIE$sweptarea*(1-VarIE$f))
  
  VarIA2=VarIA[,c(10,16)]
  VarIB2=VarIB[,c(10,16)]
  VarIC2=VarIC[,c(10,16)]
  VarID2=VarID[,c(10,16)]
  VarIE2=VarIE[,c(10,16)]
  
  prova=rbind(VarIA2,VarIB2,VarIC2)
  prova1=aggregate(prova$V,list(prova$year),sum,na.rm=TRUE)
  prova1$devst=sqrt(prova1$x)
  names(prova1)=c("year","var_shelf","stdev_shelf")
  prova2=rbind(VarID2,VarIE2)
  prova3=aggregate(prova2$V,list(prova2$year),sum,na.rm=TRUE)
  prova3$devst=sqrt(prova3$x)
  names(prova3)=c("year","var_slope","stdev_slope")
  
  VarIA3=VarIA[,c(10,17)]
  VarIB3=VarIB[,c(10,17)]
  VarIC3=VarIC[,c(10,17)]
  VarID3=VarID[,c(10,17)]
  VarIE3=VarIE[,c(10,17)]
  
  prova4=rbind(VarIA3,VarIB3,VarIC3,VarID3,VarIE3)
  prova5=aggregate(prova4$V,list(prova4$year),sum,na.rm=TRUE)
  prova5$devst=sqrt(prova5$x)
  names(prova5)=c("year","var_tot","stdev_tot")
  
  BIOMASS=merge(BIOMASS,prova1,by="year",all=TRUE)
  BIOMASS=merge(BIOMASS,prova3,by="year",all=TRUE)
  BIOMASS=merge(BIOMASS,prova5,by="year",all=TRUE)
  BIOMASS[is.na(BIOMASS)] <- 0
  yr <- seq(min(year),max(year),1)
  BIOMASS <- right_join(BIOMASS,setnames(as.data.frame(yr),"year"))
  BIOMASS[is.na(BIOMASS)] <- 0
  write.csv(BIOMASS,file=paste0(gen,spec,"_","GSA",gsa,"_",state,"_","BIOMASS.csv"))
  
  safe_ggsave(filename=paste(state,"_","GSA_",gsa,"Total_biomass.jpeg"),
         plot=ggplot(BIOMASS, aes(x=year, y=total_biomass)) + 
           geom_line(linetype = "solid",size=1.25,col="blue") +
           geom_point(colour="red",cex=3)+
           geom_errorbar(aes(ymin=total_biomass-stdev_tot,
                             ymax=total_biomass+stdev_tot),
                         width=.2,colour="red",
                         position=position_dodge(0.05),
                         linetype = "dashed",size=0.75)+
           ylab("kg/km2")+
           ggtitle(paste0(gen,spec,"_","GSA",gsa,"_",state,"Total_biomass")),
         width=7, height = 5)
  
  # DENSITY #####
  TATBDI$var1=TATBDI$sqkm*(TATBDI$N_sqkm-TATBDI$n_km2)^2
  var1a=aggregate(TATBDI$var1,list(TATBDI$year,TATBDI$strata),sum,na.rm=TRUE)
  var1a$id3=paste0(var1a$Group.1,var1a$Group.2,sep="")
  
  var2aa=merge(var1a,hauls_table2,by="id3",all=TRUE)
  var2aa$variance=(1/(var2aa$Hauls-1))*var2aa$x
  var2aa$devst=sqrt(var2aa$variance)
  var2aa=merge(var2aa,SW,by="id3",all=TRUE)
  
  var2aa$Stratum <- var2aa$stratum
  var2aa$Group.2 <- var2aa$stratum
  var2aa$Group.1 <- var2aa$year
  var2aa$Year <- var2aa$year
  var2aa[is.na(var2aa)] <- 0
  
  VarIAa=subset(var2aa,var2aa$Stratum=="A")
  VarIBa=subset(var2aa,var2aa$Stratum=="B")
  VarICa=subset(var2aa,var2aa$Stratum=="C")
  VarIDa=subset(var2aa,var2aa$Stratum=="D")
  VarIEa=subset(var2aa,var2aa$Stratum=="E")
  
  VarIAa$W1sh=W1sh
  VarIBa$W2sh=W2sh
  VarICa$W3sh=W3sh
  VarIDa$W4sl=W4sl
  VarIEa$W5sl=W5sl
  VarIAa$W1tot=W1tot
  VarIBa$W2tot=W2tot
  VarICa$W3tot=W3tot
  VarIDa$W4tot=W4tot
  VarIEa$W5tot=W5tot
  
  VarIAa$Vsh=(VarIAa$variance*VarIAa$W1sh^2)/(VarIAa$sweptarea*(1-VarIAa$f))
  VarIBa$Vsh=(VarIBa$variance*VarIBa$W2sh^2)/(VarIBa$sweptarea*(1-VarIBa$f))
  VarICa$Vsh=(VarICa$variance*VarICa$W3sh^2)/(VarICa$sweptarea*(1-VarICa$f))
  VarIDa$Vsl=(VarIDa$variance*VarIDa$W4sl^2)/(VarIDa$sweptarea*(1-VarIDa$f))
  VarIEa$Vsl=(VarIEa$variance*VarIEa$W5sl^2)/(VarIEa$sweptarea*(1-VarIEa$f))
  
  VarIAa$Vt=(VarIAa$variance*VarIAa$W1tot^2)/(VarIAa$sweptarea*(1-VarIAa$f))
  VarIBa$Vt=(VarIBa$variance*VarIBa$W2tot^2)/(VarIBa$sweptarea*(1-VarIBa$f))
  VarICa$Vt=(VarICa$variance*VarICa$W3tot^2)/(VarICa$sweptarea*(1-VarICa$f))
  VarIDa$Vt=(VarIDa$variance*VarIDa$W4tot^2)/(VarIDa$sweptarea*(1-VarIDa$f))
  VarIEa$Vt=(VarIEa$variance*VarIEa$W5tot^2)/(VarIEa$sweptarea*(1-VarIEa$f))
  
  VarIA2a=VarIAa[,c(10,16)]
  VarIB2a=VarIBa[,c(10,16)]
  VarIC2a=VarICa[,c(10,16)]
  VarID2a=VarIDa[,c(10,16)]
  VarIE2a=VarIEa[,c(10,16)]
  
  provaa=rbind(VarIA2a,VarIB2a,VarIC2a)
  prova1a=aggregate(provaa$V,list(provaa$year),sum,na.rm=TRUE)
  prova1a$devst=sqrt(prova1a$x)
  names(prova1a)=c("year","var_shelf","stdev_shelf")
  prova2a=rbind(VarID2a,VarIE2a)
  prova3a=aggregate(prova2a$V,list(prova2a$year),sum,na.rm=TRUE)
  prova3a$devst=sqrt(prova3a$x)
  names(prova3a)=c("year","var_slope","stdev_slope")
  
  VarIA3a=VarIAa[,c(10,17)]
  VarIB3a=VarIBa[,c(10,17)]
  VarIC3a=VarICa[,c(10,17)]
  VarID3a=VarIDa[,c(10,17)]
  VarIE3a=VarIEa[,c(10,17)]
  
  prova4a=rbind(VarIA3a,VarIB3a,VarIC3a,VarID3a,VarIE3a)
  prova5a=aggregate(prova4a$V,list(prova4a$year),sum,na.rm=TRUE)
  prova5a$devst=sqrt(prova5a$x)
  names(prova5a)=c("year","var_tot","stdev_tot")
  
  DENSITY=merge(DENSITY,prova1a,by="year",all=TRUE)
  DENSITY=merge(DENSITY,prova3a,by="year",all=TRUE)
  DENSITY=merge(DENSITY,prova5a,by="year",all=TRUE)
  DENSITY[is.na(DENSITY)] <- 0
  yr <- seq(min(year),max(year),1)
  DENSITY <- right_join(DENSITY,setnames(as.data.frame(yr),"year"))
  DENSITY[is.na(DENSITY)] <- 0
  
  write.csv(DENSITY,file=paste0(gen,spec,"_","GSA",gsa,"_",state,"_","DENSITY.csv"))
  
  safe_ggsave(filename=paste(state,"_","GSA_",gsa,"Total_density.jpeg"),
         plot=ggplot(DENSITY, aes(x=year, y=total_density)) + 
           geom_line(linetype = "solid",size=1,colour="blue") +
           geom_point(colour="red",cex=3)+
           geom_errorbar(aes(ymin=total_density-stdev_tot,
                             ymax=total_density+stdev_tot),
                         width=.2,colour="red",
                         position=position_dodge(0.05),
                         linetype = "dashed",size=0.75)+
           ylab("n/km2")+
           ggtitle(paste0(gen,spec,"_","GSA",gsa,"_",state,"Total_density")),
         width=7, height = 5)
  # Mean weigth by year as ratio between biomass and density indexes #######
  subTB=as.data.frame(cbind(TBn$year,TBn$haul_number,TBn$ptot,TBn$nbtot))
  names(subTB)=c("year","haul","weigth","number")
  MW=subTB$weigth/subTB$number
  subTB=cbind(subTB,MW)
  
  jpeg('boxplot_outliners.jpg')
  boxplot(subTB$MW~subTB$year,col="red",
          main=paste(gen,spec,sep=""),ylab="MW(g) by haul")
  dev.off()
  
  jpeg('boxplot_no_outliners.jpg')
  boxplot(subTB$MW~subTB$year,col="green",
          main=paste(gen,spec,sep=""),ylab="MW(g) by haul",outline=FALSE)
  dev.off()
  
  meanweight <- aggregate(subTB$weigth,by=list(subTB$year), FUN=sum) /
    aggregate(subTB$number,by=list(subTB$year), FUN=sum)
  Mweight=cbind(sort(unique(subTB$year)),meanweight)
  Mweight=Mweight[,-2]
  names(Mweight)=c("Year","Mean_weight(g)")
  Mweight <- right_join(Mweight,setnames(as.data.frame(yr),"Year"))
  write.csv(Mweight,file=paste0(gen,spec,"_GSA_",gsa,"_",state,"_Mean_weight.csv"),
            row.names=FALSE)
  Mean_weight=ggplot(data=Mweight,aes(x=Year,y=`Mean_weight(g)`))+
    geom_line()+geom_point()
  safe_ggsave(filename=paste(state,"_","GSA_",gsa,"Mean_weight.jpeg"),
         plot=Mean_weight)
  
  # Domain-specific individual mean weight from the original TB weight/number.
  # IMPORTANT: do not alter biomass/density indices or the legacy Mweight table.
  # Join species-specific TB to VALID TA hauls by the same country/GSA/year/haul
  # key used in the MEDITS engine. Each TB record contributes once only.
  medits_domain_weight <- function(ta, tb, historical_mw, observed_years) {
    result <- data.frame(Year=sort(unique(as.numeric(as.character(observed_years)))),
                         stringsAsFactors=FALSE)
    result <- result[is.finite(result$Year),,drop=FALSE]
    for (dom in c("total","shelf","slope")) {
      result[[paste0(dom,"_mean_weight")]] <- NA_real_
      result[[paste0(dom,"_N_individuals")]] <- NA_real_
      result[[paste0(dom,"_N_hauls")]] <- NA_integer_
    }
    if (is.data.frame(ta) && is.data.frame(tb) && nrow(ta) && nrow(tb)) {
      good_ta <- !is.na(ta$id2) & nzchar(as.character(ta$id2)) &
                 as.character(ta$strata) %in% c("A","B","C","D","E")
      ta_key <- as.character(ta$id2[good_ta])
      # A repeated TA haul ID would make the haul-to-depth assignment ambiguous.
      if (anyDuplicated(ta_key)) stop("Mean weight domain: duplicate valid TA haul keys")
      loc <- match(as.character(tb$id2),ta_key)
      stratum <- rep(NA_character_,nrow(tb))
      stratum[!is.na(loc)] <- as.character(ta$strata[good_ta])[loc[!is.na(loc)]]
      year <- suppressWarnings(as.numeric(as.character(tb$year)))
      wt <- suppressWarnings(as.numeric(as.character(tb$ptot)))
      nn <- suppressWarnings(as.numeric(as.character(tb$nbtot)))
      for (dom in c("total","shelf","slope")) {
        strata <- switch(dom,total=c("A","B","C","D","E"),
                             shelf=c("A","B","C"),slope=c("D","E"))
        valid <- is.finite(year) & stratum %in% strata & is.finite(wt) &
                 is.finite(nn) & wt >= 0 & nn >= 0
        if (!any(valid)) next
        vals <- data.frame(Year=year[valid],Weight=wt[valid],Number=nn[valid],
                           Haul=as.character(tb$id2[valid]))
        groups <- split(vals,vals$Year)
        for (yy in names(groups)) {
          z <- groups[[yy]]
          ix <- match(as.numeric(yy),result$Year)
          if (is.na(ix)) next
          sw <- sum(z$Weight)
          sn <- sum(z$Number)
          result[[paste0(dom,"_mean_weight")]][ix] <- if(sn > 0) sw/sn else NA_real_
          result[[paste0(dom,"_N_individuals")]][ix] <- sn
          result[[paste0(dom,"_N_hauls")]][ix] <- length(unique(z$Haul))
        }
      }
    }
    # Historical total is retained as an explicit comparison, not silently
    # substituted for the valid-haul-linked estimate.
    result$historical_total_mean_weight <- NA_real_
    if(is.data.frame(historical_mw) &&
       all(c("Year","Mean_weight(g)") %in% names(historical_mw))) {
      ix <- match(result$Year,as.numeric(as.character(historical_mw$Year)))
      result$historical_total_mean_weight <- as.numeric(as.character(
        historical_mw[["Mean_weight(g)"]][ix]))
    }
    result$total_minus_historical_g <- result$total_mean_weight -
      result$historical_total_mean_weight
    result
  }
  Mweight_domain <- medits_domain_weight(TAn,TBn,Mweight,yr)
  write.csv(Mweight_domain,
            file=paste0(gen,spec,"_GSA_",gsa,"_",state,"_Mean_weight_by_domain.csv"),
            row.names=FALSE)

  # Standardized LFDs by km2 ####
  raise=TCn$pfrac/TCn$pechan
  nblonraise=TCn$nblon*raise
  TCn=cbind(TCn,raise,nblonraise)
  
  id2=paste(TCn$country,TCn$area,TCn$year,TCn$haul_number,sep="")
  TCn=cbind(TCn,id2)
  TATCn=merge(TAn,TCn,by=c("id2","year","haul_number"),all=TRUE)
  
  TATCnn <- TATCn
  TATCnn$id3=paste0(TATCn$year.x,TATCn$strata,sep="")
  noTAwithTC <- TATCnn[is.na(TATCnn$strata),]
  write.csv(noTAwithTC,file="TChaul_no_in_TAhaul.csv",row.names = FALSE)
  
  # Check TC vs TB ####
  tempwgB <- TBn%>% group_by(country,area,year,haul_number)%>%
    summarize(totwgB=sum(ptot), .groups="drop")
  tempnbB <- TBn%>% group_by(country,area,year,haul_number)%>%
    summarize(totnbB=sum(nbtot), .groups="drop")
  tempTB <- merge(tempwgB,tempnbB,by=c("country","area","year","haul_number"),all=TRUE)
  
  tempwgC <- TCn%>% group_by(country,area,year,haul_number)%>%
    summarize(totwgC=sum(unique(pfrac)), .groups="drop")
  tempnbC <- TCn%>% group_by(country,area,year,haul_number)%>%
    summarize(totnbC=sum(nblonraise), .groups="drop")
  tempTC <- merge(tempwgC,tempnbC,by=c("country","area","year","haul_number"),all=TRUE)
  
  TBTCcheck <- merge(tempTB,tempTC,by=c("country", "area", "year", "haul_number"))
  TBTCcheck$wgratio <- TBTCcheck$totwgB/TBTCcheck$totwgC
  TBTCcheck$nbratio <- TBTCcheck$totnbB/TBTCcheck$totnbC
  
  write.csv(TBTCcheck[which(TBTCcheck$wgratio != 1 | TBTCcheck$nbratio != 1),],
            file=paste0("TBTCtoCheck_",state,gsa,gen,spec,".csv"),
            row.names = FALSE)
  
  # TOTAL LFD ####
  TATCn1=TATCn[,c("year","haul_number","strata","sqkm","length_class","nblonraise")]
  st=TATCn1$nblonraise/TATCn1$sqkm
  TATC1=cbind(TATCn1,st)
  write.table(TATCn,file=paste("TATC_",state,gsa,gen,spec,".csv"),
              sep=";",row.names=FALSE)
  
  lf=data.frame(aggregate(TATC1$nblonraise,
                          list(TATC1$length_class,TATC1$year,TATC1$strata),sum))
  names(lf)=c("LC","Year","Stratum","Value")
  lf$id3=paste0(lf$Year,lf$Stratum,sep="")
  lfst=merge(lf,SW,by="id3",all=TRUE)
  lfst$st=lfst$Value/lfst$sweptarea
  
  lfA=subset(lfst,lfst$stratum=="A")
  lfA$Wtot=W1tot; lfA[is.na(lfA)]=0; lfA$lfd=(lfA$Wtot*lfA$st)
  lfB=subset(lfst,lfst$stratum=="B")
  lfB$Wtot=W2tot; lfB[is.na(lfB)]=0; lfB$lfd=(lfB$Wtot*lfB$st)
  lfC=subset(lfst,lfst$stratum=="C")
  lfC$Wtot=W3tot; lfC[is.na(lfC)]=0; lfC$lfd=(lfC$Wtot*lfC$st)
  lfD=subset(lfst,lfst$stratum=="D")
  lfD$Wtot=W4tot; lfD[is.na(lfD)]=0; lfD$lfd=(lfD$Wtot*lfD$st)
  lfE=subset(lfst,lfst$stratum=="E")
  lfE$Wtot=W5tot; lfE[is.na(lfE)]=0; lfE$lfd=(lfE$Wtot*lfE$st)
  
  LFD=rbind(lfA,lfB,lfC,lfD,lfE)
  LFD_shelf <- medits_lfd_for_strata(LFD,c("A","B","C"))
  LFD_slope <- medits_lfd_for_strata(LFD,c("D","E"))
  LFD=data.frame(aggregate(LFD$lfd,list(LFD$LC,LFD$year),sum))
  names(LFD)=c("Length","Year","Frequency")
  LFD=LFD[,c(2,1,3)]
  LFD=LFD[!(LFD$Length==0),]
  MAX <- mround(max(LFD$Length), lfstep)
  
  # tutti gli anni del survey (da TA), anche se nel TC non c'Ã¨ nulla
  len_seq <- seq(lfstep, MAX, lfstep)
  lclasses <- expand.grid(
    Length = len_seq,
    Year   = year            # <- usa il vettore "year" di TAn
  )
  
  tempLFD <- merge(lclasses,LFD,by=c("Length","Year"),all.x=TRUE)
  errorLFD <- LFD[!LFD$Length%in%lclasses$Length,]
  write.csv(errorLFD,file=paste0("ErrorLFDtot_",state,gsa,gen,spec,".csv"),
            row.names = FALSE)
  tempLFD[is.na(tempLFD)] <- 0
  tempLFD <- tempLFD[order(tempLFD$Year,tempLFD$Length),]
  
  safe_ggsave(filename=paste(state,"_","GSA_",gsa,"LFDTOT_10-800m.jpeg"),
         width = 10, height = 8, dpi = 150, units = "in",
         plot=ggplot(tempLFD, aes(y=Frequency, x=Length,col=Year))+
           geom_bar(stat= "identity")+
           facet_wrap(~Year)+
           ggtitle(paste(gen,spec,"LFDs_10-800m_GSA",gsa,state))+
           xlab("Length")+ylab("n/km2")+
           theme(legend.position="none"))
  
  jpeg(file=paste0(gen,spec,"_GSA_",gsa,"_",state,"_","WireframeTotal.jpg"),
       width = 350, height = 350)
  theseCol=rainbow(length(tempLFD$Length))
  wireframe(Frequency~Year*Length, data = tempLFD,drape=TRUE,
            col.regions=theseCol,zlab="n/km2", xlab="year", ylab="length",
            scales = list(arrows = FALSE,tick.number = 10,col="black",
                          font= 1, tck=0.8,rot=90),
            colorkey=TRUE, ylim = rev(range(tempLFD$Length)),
            screen = list(z = 28, x = -70, y = 3),
            main=paste(gen,spec,sep="")) 
  dev.off()
  write.csv(tempLFD,file=paste0("LFDTOT_10-800m.csv"),row.names = FALSE)
  
  meanlength <- aggregate(tempLFD$Length*tempLFD$Frequency,
                          by = list(tempLFD$Year), FUN = sum) /
    aggregate(tempLFD$Frequency,
              by = list(tempLFD$Year), FUN = sum)
  
  ML <- cbind(unique(tempLFD$Year), meanlength)
  ML <- ML[,-2]
  names(ML) <- c("Year","Mean_Length_mm")   # <<< NOME PULITO
  
  write.csv(
    ML,
    file = paste0(gen,spec,"_GSA_",gsa,"_",state,"_","Mean_Length.csv"),
    row.names = FALSE
  )
  
  Mean_Length <- ggplot(ML, aes(x = Year, y = Mean_Length_mm)) +   # <<< usa il nuovo nome
    geom_line() +
    geom_point() +
    ggtitle(paste0(gen, spec, "_GSA_", gsa, "_", state)) +
    ylab("Mean length (mm)")
  
  ggsave(
    filename = paste(state,"_","GSA_",gsa,"Mean_Length.jpeg"),
    plot = Mean_Length
  )
  
  
  # FEMALE ####

  if (!any(TATCn$sex == "F", na.rm = TRUE)) {
    # nessuna femmina nei dati -> LFD F vuota ma oggetto definito
    LFDF     <- data.frame(Year = integer(0),
                           Length = numeric(0),
                           Frequency = numeric(0))
    tempLFDF <- LFDF
    LFDF_shelf <- LFDF
    LFDF_slope <- LFDF
  } else {
    TATCnF <- subset(TATCn, TATCn$sex=="F")
    TATCn0 <- TATCnF[,c("year","haul_number","strata","sqkm","length_class","nblonraise")]
    st     <- TATCn0$nblonraise/TATCn0$sqkm
    TATC0  <- cbind(TATCn0,st)
    lff    <- data.frame(aggregate(TATC0$nblonraise,
                                   list(TATC0$length_class,TATC0$year,TATC0$strata),sum))
    names(lff)=c("LC","Year","Stratum","Value")
    lff$id3=paste0(lff$Year,lff$Stratum,sep="")
    lfstf=merge(lff,SW,by="id3",all=TRUE)
    lfstf$st=lfstf$Value/lfstf$sweptarea
    
    lfAf=subset(lfstf,lfstf$stratum=="A")
    lfAf$Wtot=W1tot; lfAf[is.na(lfAf)]=0; lfAf$lfd=(lfAf$Wtot*lfAf$st)
    lfBf=subset(lfstf,lfstf$stratum=="B")
    lfBf$Wtot=W2tot; lfBf[is.na(lfBf)]=0; lfBf$lfd=(lfBf$Wtot*lfBf$st)
    lfCf=subset(lfstf,lfstf$stratum=="C")
    lfCf$Wtot=W3tot; lfCf[is.na(lfCf)]=0; lfCf$lfd=(lfCf$Wtot*lfCf$st)
    lfDf=subset(lfstf,lfstf$stratum=="D")
    lfDf$Wtot=W4tot; lfDf[is.na(lfDf)]=0; lfDf$lfd=(lfDf$Wtot*lfDf$st)
    lfEf=subset(lfstf,lfstf$stratum=="E")
    lfEf$Wtot=W5tot; lfEf[is.na(lfEf)]=0; lfEf$lfd=(lfEf$Wtot*lfEf$st)
    
    LFDF=rbind(lfAf,lfBf,lfCf,lfDf,lfEf)
    LFDF_shelf <- medits_lfd_for_strata(LFDF,c("A","B","C"))
    LFDF_slope <- medits_lfd_for_strata(LFDF,c("D","E"))
    LFDF=data.frame(aggregate(LFDF$lfd,list(LFDF$LC,LFDF$year),sum))
    names(LFDF)=c("Length","Year","Frequency")
    LFDF=LFDF[,c(2,1,3)]
    LFDF=LFDF[!(LFDF$Length==0),]
    MAXF <- mround(max(LFDF$Length), lfstep)
    
    len_seqF  <- seq(lfstep, MAXF, lfstep)
    lclassesF <- expand.grid(
      Length = len_seqF,
      Year   = year
    )
    
    tempLFDF <- merge(lclassesF,LFDF,by=c("Length","Year"),all.x=TRUE)
    tempLFDF[is.na(tempLFDF)] <- 0
    tempLFDF <- tempLFDF[order(tempLFDF$Year,tempLFDF$Length),]
    
    safe_ggsave(filename=paste(state,"_","GSA_",gsa,"FEMALE_LFDs_10-800m.jpeg"),
                width = 10, height = 8, dpi = 150, units = "in",
                plot=ggplot(tempLFDF, aes(y=Frequency, x=Length,col=Year))+
                  geom_bar(stat= "identity")+
                  facet_wrap(~Year)+
                  ggtitle(paste(gen,spec,"FEMALE_LFDs_10-800m_GSA",gsa,state))+
                  xlab("Length")+ylab("n/km2")+
                  theme(legend.position="none"))
    
    jpeg(file=paste0(gen,spec,"_GSA_",gsa,"_",state,"_","WireframeFemale.jpg"),
         width = 350, height = 350)
    theseCol=rainbow(length(tempLFDF$Length))
    wireframe(Frequency~Year*Length, data = tempLFDF,drape=TRUE,
              col.regions=theseCol,zlab="n/km2", xlab="year", ylab="length",
              scales = list(arrows = FALSE,tick.number = 10,col="black",
                            font= 1, tck=0.8,rot=90),
              colorkey=TRUE, ylim = rev(range(tempLFDF$Length)),
              screen = list(z = 28, x = -70, y = 3),
              main=paste(gen,spec,"_FEMALE",sep="")) 
    dev.off()
    write.csv(tempLFDF,file=paste0("LFDFEM_10-800m.csv"),row.names = FALSE)
  }
  
  # MALE ####
  # MALE ####
  if (!any(TATCn$sex == "M", na.rm = TRUE)) {
    LFDM     <- data.frame(Year = integer(0),
                           Length = numeric(0),
                           Frequency = numeric(0))
    tempLFDM <- LFDM
    LFDM_shelf <- LFDM
    LFDM_slope <- LFDM
  } else {
    TATCnM=subset(TATCn,TATCn$sex=="M")
    TATCn2=TATCnM[,c("year","haul_number","strata","sqkm","length_class","nblonraise")]
    st=TATCn2$nblonraise/TATCn2$sqkm
    TATC2=cbind(TATCn2,st)
    lfm=data.frame(aggregate(TATC2$nblonraise,
                             list(TATC2$length_class,TATC2$year,TATC2$strata),sum))
    names(lfm)=c("LC","Year","Stratum","Value")
    lfm$id3=paste0(lfm$Year,lfm$Stratum,sep="")
    lfstm=merge(lfm,SW,by="id3",all=TRUE)
    lfstm$st=lfstm$Value/lfstm$sweptarea
    
    lfAm=subset(lfstm,lfstm$stratum=="A")
    lfAm$Wtot=W1tot; lfAm[is.na(lfAm)]=0; lfAm$lfd=(lfAm$Wtot*lfAm$st)
    lfBm=subset(lfstm,lfstm$stratum=="B")
    lfBm$Wtot=W2tot; lfBm[is.na(lfBm)]=0; lfBm$lfd=(lfBm$Wtot*lfBm$st)
    lfCm=subset(lfstm,lfstm$stratum=="C")
    lfCm$Wtot=W3tot; lfCm[is.na(lfCm)]=0; lfCm$lfd=(lfCm$Wtot*lfCm$st)
    lfDm=subset(lfstm,lfstm$stratum=="D")
    lfDm$Wtot=W4tot; lfDm[is.na(lfDm)]=0; lfDm$lfd=(lfDm$Wtot*lfDm$st)
    lfEm=subset(lfstm,lfstm$stratum=="E")
    lfEm$Wtot=W5tot; lfEm[is.na(lfEm)]=0; lfEm$lfd=(lfEm$Wtot*lfEm$st)
    
    LFDM=rbind(lfAm,lfBm,lfCm,lfDm,lfEm)
    LFDM_shelf <- medits_lfd_for_strata(LFDM,c("A","B","C"))
    LFDM_slope <- medits_lfd_for_strata(LFDM,c("D","E"))
    LFDM=data.frame(aggregate(LFDM$lfd,list(LFDM$LC,LFDM$year),sum))
    names(LFDM)=c("Length","Year","Frequency")
    LFDM=LFDM[,c(2,1,3)]
    LFDM=LFDM[!(LFDM$Length==0),]
    
    MAXM <- mround(max(LFDM$Length), lfstep)
    
    len_seqM  <- seq(lfstep, MAXM, lfstep)
    lclassesM <- expand.grid(
      Length = len_seqM,
      Year   = year
    )
    
    tempLFDM <- merge(lclassesM,LFDM,by=c("Length","Year"),all.x=TRUE)
    tempLFDM[is.na(tempLFDM)] <- 0
    tempLFDM <- tempLFDM[order(tempLFDM$Year,tempLFDM$Length),]
    
    safe_ggsave(filename=paste(state,"_","GSA_",gsa,"MALE_LFDs_10-800m.jpeg"),
                width = 10, height = 8, dpi = 150, units = "in",
                plot=ggplot(tempLFDM, aes(y=Frequency, x=Length,col=Year))+
                  geom_bar(stat= "identity")+
                  facet_wrap(~Year)+
                  ggtitle(paste(gen,spec,"MALE_LFDs_10-800m_GSA",gsa,state))+
                  xlab("Length")+ylab("n/km2")+
                  theme(legend.position="none"))
    
    jpeg(file=paste0(gen,spec,"_GSA_",gsa,"_",state,"_","WireframeMale.jpg"),
         width = 350, height = 350)
    theseCol=rainbow(length(tempLFDM$Length))
    wireframe(Frequency~Year*Length, data = tempLFDM,drape=TRUE,
              col.regions=theseCol,zlab="n/km2", xlab="year", ylab="length",
              scales = list(arrows = FALSE,tick.number = 10,col="black",
                            font= 1, tck=0.8,rot=90),
              colorkey=TRUE, ylim = rev(range(tempLFDM$Length)),
              screen = list(z = 28, x = -70, y = 3),
              main=paste(gen,spec,"_MALE",sep=""))
    dev.off()
    write.csv(tempLFDM,file=paste0("LFDMALE_10-800m.csv"),row.names = FALSE)
  }
  
  
  LFDU <- medits_lfd_undetermined(LFD,LFDF,LFDM)
  LFDU_shelf <- medits_lfd_undetermined(LFD_shelf,LFDF_shelf,LFDM_shelf)
  LFDU_slope <- medits_lfd_undetermined(LFD_slope,LFDF_slope,LFDM_slope)
  LFDSR_shelf <- medits_lfd_sex_ratio(LFDF_shelf,LFDM_shelf)
  LFDSR_slope <- medits_lfd_sex_ratio(LFDF_slope,LFDM_slope)

  # Sex_ratio ####
  if (nrow(LFDF) == 0 || nrow(LFDM) == 0) {
    # niente dati F o M -> nessuna sex ratio calcolabile
    LFDSR <- data.frame()
  } else {
    LFDSR <- merge(LFDF, LFDM, by = c("Year","Length"), all = TRUE)
    
    ## prendo solo le lunghezze non NA
    len_non_na <- LFDSR$Length[!is.na(LFDSR$Length)]
    
    if (!length(len_non_na)) {
      # tutte le lunghezze sono NA -> nessuna sex ratio calcolabile
      LFDSR <- data.frame()
    } else {
      # griglia completa anni × lunghezze usando SOLO lunghezze finite
      len_all <- seq(lfstep, max(len_non_na), lfstep)
      
      grid_SR <- expand.grid(
        Year   = year,    # vettore year definito più sopra in funzione
        Length = len_all
      )
      
      LFDSR <- merge(grid_SR, LFDSR, by = c("Year","Length"), all.x = TRUE)
      
      # metto a zero le frequenze mancanti (F e M)
      if (!"Frequency.x" %in% names(LFDSR)) LFDSR$Frequency.x <- 0
      if (!"Frequency.y" %in% names(LFDSR)) LFDSR$Frequency.y <- 0
      
      LFDSR$Frequency.x[is.na(LFDSR$Frequency.x)] <- 0
      LFDSR$Frequency.y[is.na(LFDSR$Frequency.y)] <- 0
      
      tot <- LFDSR$Frequency.x + LFDSR$Frequency.y
      
      LFDSR$SR <- ifelse(
        tot == 0,
        NA_real_,
        LFDSR$Frequency.x / (LFDSR$Frequency.x + LFDSR$Frequency.y)
      )
    }
  }
  
  ## Modelling sex ratio vector ##
  if (is.null(LFDSR) || nrow(LFDSR) == 0 || !"Length" %in% names(LFDSR)) {
    # Nessun dato F/M -> nessun vettore di SR stimabile
    SR_vec <- data.frame(Length = numeric(0), sr = numeric(0))
  } else {
    ## uso solo righe con Length e SR non NA
    LFDSR_mod <- LFDSR[!is.na(LFDSR$Length) & !is.na(LFDSR$SR), ]
    
    if (nrow(LFDSR_mod) == 0) {
      # ci sono righe in LFDSR ma tutte SR NA -> niente modello
      SR_vec <- data.frame(Length = numeric(0), sr = numeric(0))
    } else {
      # plot "interno" (puoi anche toglierlo, non è usato)
      ggplot(LFDSR_mod, aes(x = Length, y = SR)) +
        geom_point() +
        stat_smooth(method = "loess", fullrange = TRUE, span = 0.6) +
        ggtitle("SexRatio") + xlab("Length") + ylab("FF/FF+MM") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.6)) +
        scale_x_continuous(
          breaks = seq(0, max(LFDSR_mod$Length, na.rm = TRUE), by = 10)
        ) +
        theme(legend.position = "none")
      
      ## modello loess solo sui dati "buoni"
      smooth_vals <- predict(
        loess(SR ~ Length, data = LFDSR_mod),
        newdata = data.frame(
          Length = sort(unique(LFDSR_mod$Length))
        )
      )
      
      smooth_vals[which(smooth_vals < 0)] <- 0
      smooth_vals[which(smooth_vals > 1)] <- 1
      
      SR_vec <- data.frame(
        Length = sort(unique(LFDSR_mod$Length)),
        sr     = smooth_vals
      )
      
      SR_vec <- SR_vec[order(SR_vec$Length), ]
      SR_vec <- SR_vec[!SR_vec$Length %in% 37, ]   # tua riga originale
      
      # griglia completa di lunghezze basata solo su valori non NA
      len_non_na2 <- SR_vec$Length[!is.na(SR_vec$Length)]
      if (!length(len_non_na2)) {
        SR_vec <- data.frame(Length = numeric(0), sr = numeric(0))
      } else {
        max_len <- max(len_non_na2)
        lenn    <- data.frame(Length = seq(lfstep, max_len, lfstep))
        
        SR_vec <- dplyr::full_join(lenn, SR_vec, by = "Length")
        
        pos <- which(SR_vec$sr > 0)[1] - 1
        if (!is.na(pos) && pos > 0) {
          SR_vec$sr[1:pos] <- 0.5
        }
        SR_vec$sr[is.na(SR_vec$sr)] <- max_SR
      }
    }
  }
  
  
  
  # The overall run's sex-ratio curve is available to the per-GSA run.
  # Preserve GSA-specific observed/raised LFDs; share only the curve.
  if (!is.null(global_sex_ratio) &&
      is.data.frame(global_sex_ratio) &&
      all(c("Length", "sr") %in% names(global_sex_ratio)) &&
      nrow(global_sex_ratio)) {
    SR_vec <- global_sex_ratio[, c("Length", "sr"), drop = FALSE]
  }

  ## Yearly density and biomass per hauls bubble plots ####
  metallica=merge(TAn,TATBn,by="id2")
  abundance=cbind.data.frame(metallica$country,metallica$year,
                             metallica$TA_lat_fin,metallica$TA_lon_fin,
                             metallica$meandepth.x,metallica$sqkm.x,
                             metallica$ptot,metallica$nbtot)
  colnames(abundance)=c("COUNTRY","YEAR","LAT","LONG","DEPTH",
                        "SQKM","WEIGHT","NUMBER")
  abundance[is.na(abundance)] = 0
  abundance$BIOMASS=(abundance$WEIGHT/1000)/abundance$SQKM
  abundance$DENSITY=abundance$NUMBER/abundance$SQKM
  
  ## ---------------------------------------------------------------
  ##  GAM: trend standardizzato con stesso "design" di cale per anno
  ##  (usa tutte le cale reali per stimare il modello)
  ## ---------------------------------------------------------------
  
  # BIOMASS_gam_index <- NULL
  # DENSITY_gam_index <- NULL
  # mod_biom_gam      <- NULL
  # mod_dens_gam      <- NULL
  # 
  # # dati a livello di cala per il GAM
  # gam_base <- abundance %>%
  #   dplyr::filter(!is.na(BIOMASS),
  #                 !is.na(DENSITY),
  #                 SQKM > 0,
  #                 !is.na(YEAR),
  #                 !is.na(LAT),
  #                 !is.na(LONG),
  #                 !is.na(DEPTH))
  # 
  # if (nrow(gam_base) > 0) {
  #   
  #   ## "design" di riferimento: stesse posizioni (LAT, LONG, DEPTH) per tutti gli anni
  #   design_ref <- gam_base %>%
  #     dplyr::select(LAT, LONG, DEPTH) %>%
  #     dplyr::distinct()
  #   hauls_per_year <- gam_base %>%
  #     dplyr::group_by(YEAR) %>%
  #     dplyr::summarise(n_hauls = dplyr::n(), .groups = "drop")
  #   
  #   cat("GAM - numero di cale reali per anno:\n")
  #   print(hauls_per_year)
  #   # se è enorme puoi ridurre:
  #   # design_ref <- design_ref %>% dplyr::sample_n(size = min(500, nrow(design_ref)))
  #   n_ref <- nrow(design_ref)
  #   cat("GAM - numero di cale standard (design_ref):", n_ref, "\n")
  #   ## preparo dati con variabile anno numerica
  #   dat_biom <- gam_base %>%
  #     dplyr::mutate(YEAR_num = as.numeric(YEAR))
  #   
  #   dat_dens <- gam_base %>%
  #     dplyr::mutate(YEAR_num = as.numeric(YEAR))
  #   
  #   ## GAM per BIOMASS (kg/km2) - usa tutte le cale reali
  #   mod_biom_gam <- mgcv::gam(
  #     BIOMASS ~ s(YEAR_num, k = 10) +
  #       s(LONG, LAT, k = 100) +
  #       s(DEPTH, k = 7),
  #     data   = dat_biom,
  #     family = mgcv::tw(link = "log"),  # Tweedie: continua con zeri
  #     method = "REML"
  #   )
  #   
  #   ## GAM per DENSITY (n/km2)
  #   mod_dens_gam <- mgcv::gam(
  #     DENSITY ~ s(YEAR_num, k = 10) +
  #       s(LONG, LAT, k = 100) +
  #       s(DEPTH, k = 7),
  #     data   = dat_dens,
  #     family = mgcv::tw(link = "log"),
  #     method = "REML"
  #   )
  #   
  #   ## Predizione di indici annuali standardizzati
  #   years_vec <- sort(unique(gam_base$YEAR))
  #   
  #   # costruisco la stessa griglia di "pseudo-cale" per tutti gli anni
  #   new_all <- do.call(
  #     rbind,
  #     lapply(years_vec, function(yy) {
  #       data.frame(
  #         YEAR      = yy,
  #         YEAR_num  = as.numeric(yy),
  #         LAT       = design_ref$LAT,
  #         LONG      = design_ref$LONG,
  #         DEPTH     = design_ref$DEPTH
  #       )
  #     })
  #   )
  #   
  #   # uso la stessa griglia sia per biomassa che per densità
  #   pred_biom <- predict(mod_biom_gam, newdata = new_all, type = "response")
  #   pred_dens <- predict(mod_dens_gam, newdata = new_all, type = "response")
  #   
  #   new_all$pred_biom <- pred_biom
  #   new_all$pred_dens <- pred_dens
  #   
  #   # indice annuale = media delle predizioni sulla griglia standard
  #   n_ref <- nrow(design_ref)  # numero di pseudo-cale nel design standard
  #   
  #   BIOMASS_gam_index <- new_all %>%
  #     dplyr::group_by(YEAR) %>%
  #     dplyr::summarise(
  #       index_biom_gam = mean(pred_biom, na.rm = TRUE),
  #       sd_biom_gam    = sd(pred_biom,  na.rm = TRUE) / sqrt(n_ref),  # <<< SE dell'indice
  #       .groups = "drop"
  #     )
  #  
  #   DENSITY_gam_index <- new_all %>%
  #     dplyr::group_by(YEAR) %>%
  #     dplyr::summarise(
  #       index_dens_gam = mean(pred_dens, na.rm = TRUE),
  #       sd_dens_gam    = sd(pred_dens,  na.rm = TRUE) / sqrt(n_ref),
  #       .groups = "drop")
  #     
  # }
  # 
  #
  # DENSITY bubbles ##
  DIbubble=subset(abundance[!is.na(abundance$DENSITY),])
  DIbubble$levels_median=NA
  DIbubble=subset(abundance[abundance$DENSITY>0,])
  plot_list5 = list()
  
  for (i in min(DIbubble$YEAR):max(DIbubble$YEAR)){
    DIbubble2 <- DIbubble[DIbubble$YEAR==i,]
    DIbubble2medianlow=quantile(DIbubble2$DENSITY,prob=(0.25))
    DIbubble2medianup=quantile(DIbubble2$DENSITY,prob=(0.75))
    DIbubble2$levels_median <-ifelse(DIbubble2$DENSITY <= DIbubble2medianlow,"low",
                                     ifelse(DIbubble2$DENSITY > DIbubble2medianlow  &
                                              DIbubble2$DENSITY <= DIbubble2medianup,
                                            "medium","high"))
    
    plot_list5[[i]]=ggplot(data = DIbubble2,
                           aes(x =LONG, y= LAT,cex=DENSITY,col=levels_median))+
      borders("worldHires",fill="lightblue",colour="black",
              xlim=xlims,ylim=ylims)+
      geom_point()+
      scale_colour_manual(values=c("low"="yellow",
                                   "medium"="orange",
                                   "high" ="red"))+
      coord_quickmap(xlim = xlims,ylim=ylims)+
      labs(x="Longitude", y="Latitude", cex="N/sqkm", col="Level")+
      ggtitle(paste0("Density index_",state,"_","GSA",gsa," ",
                     gen,spec," ",i))
  }
  
  for (i in min(DIbubble$YEAR):max(DIbubble$YEAR)) {
    file_name = paste0("DENSITY_BUBBLE_YR_",gen,spec,"_GSA_",gsa,"_",
                       state,"_",i, ".tiff")
    tiff(file_name, width = 800, height = 600, units = "px")
    print(plot_list5[[i]])
    dev.off()
  }
  
  # Expose the already-computed per-haul fields for the Shiny maps.
  # IMPORTANT: these are observed per-haul values, NOT area-stratified annual indices.
  # Keep valid zero catches where present; never replace missing hauls with zero.
  spatial_ok <- is.finite(suppressWarnings(as.numeric(abundance$SQKM))) &
    suppressWarnings(as.numeric(abundance$SQKM)) > 0 &
    is.finite(suppressWarnings(as.numeric(abundance$LAT))) &
    is.finite(suppressWarnings(as.numeric(abundance$LONG))) &
    is.finite(suppressWarnings(as.numeric(abundance$YEAR)))
  medits_spatial <- abundance[spatial_ok, , drop=FALSE]
  BIbubble_map <- medits_spatial[is.finite(medits_spatial$BIOMASS) &
                                 medits_spatial$BIOMASS >= 0, , drop=FALSE]
  DIbubble_map <- medits_spatial[is.finite(medits_spatial$DENSITY) &
                                 medits_spatial$DENSITY >= 0, , drop=FALSE]

  ## BIOMASS bubbles ##
  BIbubble=subset(abundance[!is.na(abundance$BIOMASS),])
  BIbubble$levels_median=NA
  BIbubble=subset(abundance[abundance$BIOMASS>0,])
  
  plot_list6 = list()
  
  for (i in min(BIbubble$YEAR):max(BIbubble$YEAR)){
    BIbubble2 <- BIbubble[BIbubble$YEAR==i,]
    BIbubble2medianlow=quantile(BIbubble2$BIOMASS,prob=(0.25))
    BIbubble2medianup=quantile(BIbubble2$BIOMASS,prob=(0.75))
    BIbubble2$levels_median <-ifelse(BIbubble2$BIOMASS <= BIbubble2medianlow,"low",
                                     ifelse(BIbubble2$BIOMASS > BIbubble2medianlow  &
                                              BIbubble2$BIOMASS <= BIbubble2medianup,
                                            "medium","high"))
    
    plot_list6[[i]]=ggplot(data = BIbubble2,
                           aes(x =LONG, y= LAT,cex=BIOMASS,col=levels_median))+
      borders("worldHires",fill="lightblue",colour="black",
              xlim=xlims,ylim=ylims)+
      geom_point()+
      scale_colour_manual(values=c("low"="yellow",
                                   "medium"="orange",
                                   "high" ="red"))+
      coord_quickmap(xlim = xlims,ylim=ylims)+
      labs(x="Longitude", y="Latitude", cex="kg/sqkm", col="Level")+
      ggtitle(paste0("Biomass index_",state,"_","GSA",gsa," ",
                     gen,spec," ",i))
  }
  
  for (i in min(BIbubble$YEAR):max(BIbubble$YEAR)) {
    file_name = paste0("BIOMASS_BUBBLE_YR_",gen,spec,"_GSA_",gsa,"_",
                       state,"_",i, ".tiff")
    tiff(file_name, width = 800, height = 600, units = "px")
    print(plot_list6[[i]])
    dev.off()
  }
 
  setwd(WD)
  ######### END OF SCRIPT ###############
  
  ## -----------------------------------------------------------------
  ## ALLA FINE DELLA FUNZIONE RESTITUIAMO GLI OGGETTI CHE SERVONO
  ## AI GRAFICI STANDARD
  ## -----------------------------------------------------------------
  out <- list(
    TAn          = TAn,
    BIbubble     = BIbubble,
    BIbubble_map = BIbubble_map,
    DIbubble_map = DIbubble_map,
    BIOMASS      = BIOMASS,
    DENSITY      = DENSITY,
    # BIOMASS_gam_index = BIOMASS_gam_index,  # <<< NUOVO
    # DENSITY_gam_index = DENSITY_gam_index,  # <<< NUOVO
    # mod_biom_gam      = mod_biom_gam,      # <<< NUOVO (modello completo)
    # mod_dens_gam      = mod_dens_gam,      # <<< NUOVO
    Mweight      = Mweight,
    Mweight_domain = Mweight_domain,
    ML           = ML,
    LFD          = LFD,
    LFD_shelf    = LFD_shelf,
    LFD_slope    = LFD_slope,
    LFDF         = LFDF,
    LFDF_shelf   = LFDF_shelf,
    LFDF_slope   = LFDF_slope,
    LFDM         = LFDM,
    LFDM_shelf   = LFDM_shelf,
    LFDM_slope   = LFDM_slope,
    LFDU         = LFDU,
    LFDU_shelf   = LFDU_shelf,
    LFDU_slope   = LFDU_slope,
    LFDSR        = LFDSR,
    LFDSR_shelf  = LFDSR_shelf,
    LFDSR_slope  = LFDSR_slope,
    SR_vec       = SR_vec,
    gsa          = gsa,
    state        = state,
    gen          = gen,
    spec         = spec,
    lfstep       = lfstep,
    max_SR       = max_SR,
    TCn          = TCn,
    output_run_dir = output_run_dir
  )
  
  
  return(out)
}