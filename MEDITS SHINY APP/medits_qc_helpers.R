# MEDITS data-quality helpers – original-unit, record-level, non-destructive.
# No automatic exclusions, no FE changes, and no assumption that maturity = age.

medits_qc_norm <- function(x) toupper(gsub('[^A-Za-z0-9]', '', as.character(x)))

medits_qc_catalog <- function(data) {
  targets <- list(
    TC=list(
      LENGTH=c('LENGTH_CLASS','LENGTHCLASS','LENGTH_MM','LENGTH','LEN','LENGTH_CM'),
      WEIGHT=c('WEIGHT_OF_THE_FRACTION','WEIGHT_OF_THE_SAMPLE_MEASURED','PFRAC','PECHAN',
               'INDIVIDUAL_WEIGHT','WEIGHT_G','IND_WEIGHT'),
      AGE=c('AGE','AGE_YEARS','AGE_IN_YEARS','AGE_GROUP','AGE_CLASS','ETA','ETÀ')),
    TB=list(
      WEIGHT=c('TOTAL_WEIGHT_IN_THE_HAUL','PTOT','TOTAL_WEIGHT','TOTAL_CATCH_WEIGHT')),
    AGE=list(
      AGE=c('AGE','AGE_YEARS','AGE_IN_YEARS','AGE_GROUP','AGE_CLASS','ETA','ETÀ'),
      LENGTH=c('LENGTH_CLASS','LENGTH_MM','LENGTH','LENGTH_CM','LEN'),
      WEIGHT=c('INDIVIDUAL_WEIGHT','WEIGHT_G','WEIGHT','WEIGHT_GRAMS','IND_WEIGHT')))
  ans <- list()
  for (source in intersect(names(targets),names(data))) {
    d <- data[[source]]
    if (is.null(d) || !is.data.frame(d) || !nrow(d)) next
    nms <- names(d); norm <- medits_qc_norm(nms)
    used <- character()
    for (typ in names(targets[[source]])) {
      for (target in targets[[source]][[typ]]) {
        idx <- which(norm==medits_qc_norm(target))
        if (length(idx)) {
          field <- nms[idx[1]]
          if (!(field %in% used)) {
            ans[[length(ans)+1L]] <- data.frame(Source=source,Field=field,Type=typ,
                                                  stringsAsFactors=FALSE)
            used <- c(used,field)
          }
        }
      }
    }
  }
  if (!length(ans)) return(data.frame(Source=character(),Field=character(),Type=character(),
                                      stringsAsFactors=FALSE))
  do.call(rbind,ans)
}

medits_qc_meta <- function(d, choices) {
  ix <- match(medits_qc_norm(choices),medits_qc_norm(names(d)))
  ix <- ix[!is.na(ix)]
  if (!length(ix)) return(rep('',nrow(d)))
  as.character(d[[ix[1]]])
}

medits_qc_scan_data <- function(data, max_length, max_weight, max_age,
                                include_missing=FALSE, selected_type='ALL',
                                selected_rows=NULL) {
  catalog <- medits_qc_catalog(data)
  stats <- list(); found <- list()
  for (j in seq_len(nrow(catalog))) {
    src <- catalog$Source[j]; field <- catalog$Field[j]; typ <- catalog$Type[j]
    d <- data[[src]]
    # Preserve the original source-row numbers used by the correction register.
    row_ids <- if (!is.null(selected_rows) && src %in% names(selected_rows))
      as.integer(selected_rows[[src]]) else seq_len(nrow(d))
    row_ids <- row_ids[!is.na(row_ids) & row_ids >= 1L & row_ids <= nrow(d)]
    d <- d[row_ids,,drop=FALSE]
    v <- suppressWarnings(as.numeric(gsub(',', '.', as.character(d[[field]]),fixed=TRUE)))
    missing <- !is.finite(v)
    limit <- switch(typ,LENGTH=max_length,WEIGHT=max_weight,AGE=max_age)
    invalid <- if (typ=='LENGTH') is.finite(v) & v<=0 else is.finite(v) & v<0
    above <- is.finite(v) & v>limit
    flagged <- invalid | above | (isTRUE(include_missing) & missing)
    stats[[length(stats)+1L]] <- data.frame(
      Source=src,Field=field,Type=typ,Rows=nrow(d),
      Non_numeric_or_missing=sum(missing),Invalid_or_negative=sum(invalid),
      Above_limit=sum(above),Flagged=sum(flagged),
      Min=if(any(is.finite(v))) min(v[is.finite(v)]) else NA_real_,
      Max=if(any(is.finite(v))) max(v[is.finite(v)]) else NA_real_,
      Manual_limit=limit,stringsAsFactors=FALSE)
    if (selected_type!='ALL' && typ!=selected_type) next
    ix <- which(flagged)
    if (!length(ix)) next
    reason <- ifelse(missing[ix],'missing or nonnumeric',
                     ifelse(invalid[ix],if(typ=='LENGTH') 'nonpositive length' else 'negative value',
                            paste0('above manual limit: ',limit)))
    found[[length(found)+1L]] <- data.frame(
      Source=src,Type=typ,Row=row_ids[ix],Field=field,Value=v[ix],Reason=reason,
      Country=medits_qc_meta(d,c('COUNTRY'))[ix],
      GSA=medits_qc_meta(d,c('AREA','GSA'))[ix],
      Year=medits_qc_meta(d,c('YEAR'))[ix],
      Haul=medits_qc_meta(d,c('HAUL_NUMBER','HAUL','HAULID'))[ix],
      Species=medits_qc_meta(d,c('SPECIES','SPECIES_CODE','SPP'))[ix],
      Sex=medits_qc_meta(d,c('SEX'))[ix],
      Length=medits_qc_meta(d,c('LENGTH_CLASS','LENGTH','LENGTH_MM'))[ix],
      Age=medits_qc_meta(d,c('AGE','AGE_YEARS','AGE_CLASS','AGE_GROUP','ETA'))[ix],
      stringsAsFactors=FALSE)
  }
  cols <- c('ID','Source','Type','Row','Field','Value','Reason','Country','GSA',
            'Year','Haul','Species','Sex','Length','Age')
  candidates <- if (length(found)) do.call(rbind,found) else
    data.frame(Source=character(),Type=character(),Row=integer(),Field=character(),
               Value=numeric(),Reason=character(),Country=character(),GSA=character(),
               Year=character(),Haul=character(),Species=character(),Sex=character(),
               Length=character(),Age=character(),stringsAsFactors=FALSE)
  candidates$ID <- seq_len(nrow(candidates))
  candidates <- candidates[,cols,drop=FALSE]
  summary <- if (length(stats)) do.call(rbind,stats) else data.frame()
  list(candidates=candidates,fields=summary,catalog=catalog)
}

medits_qc_subset_candidates <- function(df,query='') {
  if (!nrow(df) || is.null(query) || !nzchar(trimws(query))) return(df)
  hay <- apply(df[,intersect(c('Source','Type','Row','Field','Value','Country','GSA',
                               'Year','Haul','Species','Sex','Reason'),names(df)),drop=FALSE],
               1, paste,collapse=' | ')
  df[grepl(trimws(query),hay,ignore.case=TRUE,fixed=TRUE),,drop=FALSE]
}


# Strict species selection follows the same lookup used by run_medits_from_raw.
# Data are NOT subsetted in memory for analysis: these are row positions for QC only.
medits_qc_species_lookup <- function(alpha, species_file='Sp_Medits.csv',
                                     asfis_file='ASFIS_2020.csv') {
  alpha <- toupper(trimws(as.character(alpha[1])))
  if (!nzchar(alpha) || is.na(alpha)) stop('Select a valid 3-alpha species code before biological review.')
  if (!file.exists(species_file) || !file.exists(asfis_file))
    stop('Species lookup missing: place Sp_Medits.csv and ASFIS_2020.csv in the application directory.')
  mdts <- data.table::fread(species_file,showProgress=FALSE)
  asfis <- data.table::fread(asfis_file,showProgress=FALSE)
  if (ncol(mdts)<7L || ncol(asfis)<4L) stop('Invalid species lookup CSV column structure.')
  mdts <- mdts[,c(2,3,7),with=FALSE]
  asfis <- asfis[,c(3,4),with=FALSE]
  names(mdts)[2] <- 'Scientific_name'
  lookup <- merge(mdts,asfis,by='Scientific_name',all=TRUE)
  lookup <- lookup[!is.na(lookup$code) & !is.na(lookup$`3A_CODE`),,drop=FALSE]
  matched <- lookup[toupper(trimws(as.character(lookup$`3A_CODE`)))==alpha,,drop=FALSE]
  codes <- as.character(matched$code)
  pairs <- unique(data.frame(GENUS=substr(codes,1L,4L),
                             SPECIES=substr(codes,5L,7L),stringsAsFactors=FALSE))
  pairs <- pairs[nzchar(pairs$GENUS) & nzchar(pairs$SPECIES),,drop=FALSE]
  if (!nrow(pairs)) stop('No MEDITS GENUS/SPECIES mapping found for species: ',alpha)
  pairs
}

medits_qc_species_rows <- function(data, alpha, gsa=NULL, country=NULL,
                                   species_file='Sp_Medits.csv',
                                   asfis_file='ASFIS_2020.csv') {
  pairs <- medits_qc_species_lookup(alpha,species_file,asfis_file)
  alpha <- toupper(trimws(as.character(alpha[1])))
  norm <- function(v) toupper(trimws(as.character(v)))
  gsas <- norm(gsa); gsas <- gsas[!is.na(gsas) & nzchar(gsas)]
  countries <- norm(country); countries <- countries[!is.na(countries) & nzchar(countries)]
  result <- list()
  notices <- character()
  for (src in intersect(c('TB','TC','AGE'),names(data))) {
    d <- data[[src]]
    if (!is.data.frame(d)) next
    n <- nrow(d)
    keep <- rep(TRUE,n)
    genus <- medits_qc_meta(d,c('GENUS'))
    species <- medits_qc_meta(d,c('SPECIES'))
    if ('GENUS' %in% toupper(names(d)) && 'SPECIES' %in% toupper(names(d))) {
      keep <- keep & (paste(norm(genus),norm(species),sep='|') %in%
           paste(norm(pairs$GENUS),norm(pairs$SPECIES),sep='|'))
    } else {
      code <- medits_qc_meta(d,c('SPECIES_CODE','3A_CODE','ASFIS_CODE','SPP'))
      # Some age-reading tables use SPECIES for the 3-alpha code, not the taxon abbreviation.
      if (all(!nzchar(code) | is.na(code)) && src=='AGE')
        code <- medits_qc_meta(d,c('SPECIES'))
      if (all(!nzchar(code) | is.na(code))) {
        keep[] <- FALSE
        notices <- c(notices,paste0(src,': no species identifier found. Skipped to prevent mixing species.'))
      } else keep <- keep & norm(code)==alpha
    }
    if(length(gsas)) {
      area <- medits_qc_meta(d,c('AREA','GSA'))
      if (all(!nzchar(area) | is.na(area))) {
        keep[] <- FALSE
        notices <- c(notices,paste0(src,': missing GSA/AREA; excluded from scoped review.'))
      } else keep <- keep & norm(area) %in% gsas
    }
    if(length(countries)) {
      cty <- medits_qc_meta(d,c('COUNTRY'))
      if (all(!nzchar(cty) | is.na(cty))) {
        keep[] <- FALSE
        notices <- c(notices,paste0(src,': missing COUNTRY; excluded from scoped review.'))
      } else keep <- keep & norm(cty) %in% countries
    }
    keep[is.na(keep)] <- FALSE
    result[[src]] <- which(keep)
  }
  attr(result,'notices') <- unique(notices)
  result
}