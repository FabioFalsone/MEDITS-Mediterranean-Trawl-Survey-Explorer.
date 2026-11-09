## Medits_app1.R

options(shiny.maxRequestSize = 400 * 1024^2)  # 400 MB

library(shiny)
library(data.table)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(reshape2)
library(FLCore)
library(ggplotFL)
library(dplyr)
library(maps)
library(mapdata)
source("Medits_function_GSAMedits_gsa.R")  # definisce run_medits_from_raw()
source("medits_report_helpers.R")  # report PDF: figure, indice, tabelle e descrizioni numeriche
source("medits_qc_helpers.R")       # record-level QC lengths, weights, optional ages
source("medits_length_stats.R")       # robust LFD summaries, safely handles one-bin years
source("medits_display_helpers.R")   # robust map and LFD plotting (no index recalculation)
source("medits_legacy_lfd_plot.R") # real original plotLengths() LFD layout (v4)

# Fail fast if an old/incomplete functions file was copied into this folder.
# The same helpers are used by Shiny maps, Survey period, and the PDF report.
medits_required_survey_helpers <- c(
  "medits_coast_segments", "medits_survey_unique_hauls",
  "medits_survey_period_summary", "medits_make_survey_period_plot"
)
medits_missing_survey_helpers <- medits_required_survey_helpers[
  !vapply(medits_required_survey_helpers, exists, logical(1),
          mode="function", inherits=TRUE)
]
if (length(medits_missing_survey_helpers)) {
  stop("Missing MEDITS survey/map functions: ",
       paste(medits_missing_survey_helpers, collapse=", "),
       ". Replace Medits_function_GSAMedits_gsa.R with the file supplied ",
       "in this FIX package, then restart RStudio/Shiny.")
}
message("MEDITS v5 OCEAN THEME: survey, map and biological review helpers loaded successfully")

## <<< NUOVO: salvo la cartella della app >>>
app_dir <- normalizePath(".")


ui <- fluidPage(
  tags$head(tags$link(rel="stylesheet",type="text/css",href="medits_theme.css")),
  tags$div(class="medits-hero",
    tags$div(class="eyebrow", icon("ship"), "MARINE RESEARCH  /  SUSTAINABLE FISHERIES  /  MEDITS"),
    tags$h1("MEDITS · Mediterranean Trawl Survey Explorer"),
    tags$p("This application uses the R scripts developed by Mannini et al.2020 for the analysis and estimation of MEDITS survey indicators.")),
  sidebarLayout(
    sidebarPanel(
      fileInput("ta_file", "Upload TA.csv", accept = ".csv"),
      fileInput("tb_file", "Upload TB.csv", accept = ".csv"),
      fileInput("tc_file", "Upload TC.csv", accept = ".csv"),
      fileInput("age_file", "Age readings CSV (optional; quality review only)", accept = ".csv"),
      tags$div(class="medits-sidebar-title", "Analysis settings"),
      # SPECIE: 3-alpha da input
      textInput("species", "Species (3-alpha, e.g. HKE)", value = "HKE"),
      uiOutput("gsa_ui"),
      uiOutput("country_ui"),
      uiOutput("display_gsa_ui"),
      tags$small("Choose AGGREGATED or an individual GSA: the selection changes every existing trend, figure, CSV and PDF."),

      
      ## Scelta dello strato SOLO per la visualizzazione dei plot
      radioButtons(
        "stratum_view",
        "Bathymetric domain (plots, LFD and PDF)",
        choices = c(
          "Total (shelf + slope)"   = "total",
          "Shelf only (10-200 m)"   = "shelf",
          "Slope only (200-800 m)" = "slope"
        ),
        selected = "total"
      ),
      
      textInput("survey", "Survey name", value = "MEDITS"),
      tags$hr(),
      tags$div(class="medits-sidebar-title", "Review and analysis"),
      actionButton("open_bio_qc", tagList(icon("search")," Review length, age and weight"),
                   class="btn-default medits-full-btn"),
      tags$p(class="medits-sidebar-note", "Manual screening of length, weight and age values only. No observation is removed automatically."),
      actionButton("run_btn", tagList(icon("play-circle"), " Run MEDITS analysis"),
                   class="btn-primary medits-run-btn medits-full-btn"),
      tags$div(class="medits-mini-status", tags$span(class="medits-status-dot"),
               textOutput("medits_run_status",inline=TRUE)),
      tags$hr(),
      checkboxInput("report_include_interpretation",
                    "Add automatic numerical descriptions", value = TRUE),
      checkboxInput("report_include_detailed_tables",
                    "Append detailed LFD and sex-ratio matrices (compact tables)", value = FALSE),
      downloadButton("download_full_report", "Download complete PDF report"),
      tags$hr(),
      downloadButton("download_all_csv_jpeg", "Download ALL CSV + JPEG + report figures (ZIP)"),
      br(), br(),
      downloadButton("download_lfd_gsa", "Download LFD by individual GSA (ZIP)")
    ),
    
    mainPanel(
      tabsetPanel(id="main_medits_tabs",
        tabPanel("Survey", value="survey", icon=icon("ship"),
          tags$div(class="medits-group-header",
            tags$div(class="medits-kicker","MARINE DATA · MEDITS"),
            tags$h2("Survey"),
            tags$p("Survey timing, sampling effort and haul positions · a transparent view of the data collection footprint.")),
          tabsetPanel(id="survey_subtabs",type="pills",
          tabPanel(
          "Survey period",
          uiOutput("survey_plot_ui"),
          br(),
          downloadButton("dl_surveyp_png", "Download PNG"),
          downloadButton("dl_surveyp_csv", "Download CSV")
        ),
          tabPanel(
          "Haul map by Year",
          br(),
          checkboxInput("haul_show_all_years", "Show all survey years (comparable panels)", value = TRUE),
          conditionalPanel("!input.haul_show_all_years", uiOutput("year_map_ui")),
          plotOutput("map_hauls", height = "880px"),
          br(),
          downloadButton("dl_haulmap_png", "Download PNG"),
          downloadButton("dl_haulmap_csv", "Download CSV")
        ),
          tabPanel(
          "Hauls by stratum",
          plotOutput("hauls_by_stratum", height = "500px"),
          br(),
          downloadButton("dl_haulstrat_csv", "Download CSV")
        )
          )
        ),
        tabPanel("Indices & maps", value="indices_and_maps", icon=icon("line-chart"),
          tags$div(class="medits-group-header",
            tags$div(class="medits-kicker","MARINE DATA · MEDITS"),
            tags$h2("Indices & maps"),
            tags$p("Biomass and density from stratified MEDITS analyses · aggregated and individual GSA comparisons.")),
          tabsetPanel(id="indices_subtabs",type="pills",
          tabPanel(
          "Biomass / Density / Mean W",
          h4("Biomass"),
          plotOutput("biom_trend",  height = "300px"),
          downloadButton("dl_biomtrend_png", "PNG"),
          downloadButton("dl_biomtrend_csv", "CSV"),
          br(), br(),
          
          h4("Density"),
          plotOutput("dens_trend",  height = "300px"),
          downloadButton("dl_denstrend_png", "PNG"),
          downloadButton("dl_denstrend_csv", "CSV"),
          br(), br(),
          
          h4("Mean weight"),
          plotOutput("mean_weight", height = "300px"),
          downloadButton("dl_meanw_png", "PNG"),
          downloadButton("dl_meanw_csv", "CSV")
        ),
          tabPanel("Compare GSA trends",
          tags$h3("Aggregate and individual GSA trends"),
          tags$p("Biomass and density follow the selected bathymetric domain. Mean weight is an observed annual measure, not macro-stratum-specific. Missing GSA/year values remain absent, never imputed."),
          plotOutput("gsa_comparison_plot",height="850px"),
          downloadButton("dl_gsa_comparison_png","Download comparison figure (PNG)"),
          downloadButton("dl_gsa_comparison_csv","Download trend data: aggregate + each GSA (CSV)"),
          tags$h4("Status of individual GSA calculations"),
          tableOutput("gsa_run_table")),
          tabPanel(
          "Biomass map (log)",
          plotOutput("map_biomass", height = "1050px", width = "100%"),
          br(),
          downloadButton("dl_biomap_png", "Download PNG"),
          downloadButton("dl_biomap_csv", "Download CSV")
        ),
          tabPanel(
          "Density map (log)",
          plotOutput("map_density", height = "1050px", width = "100%"),
          br(),
          downloadButton("dl_denmap_png", "Download PNG"),
          downloadButton("dl_denmap_csv", "Download CSV")
        )
          )
        ),
        tabPanel("Length & sex", value="length_and_sex", icon=icon("leaf"),
          tags$div(class="medits-group-header",
            tags$div(class="medits-kicker","MARINE DATA · MEDITS"),
            tags$h2("Length & sex"),
            tags$p("Classic ggplotFL plotLengths() figures (two year blocks, blue mean and red L95), with observed and LOESS-corrected sex frequencies.")),
          radioButtons("lfd_sex_mode", "Female / male frequency display",
            choices=c("LOESS corrected (includes allocated undetermined)"="corrected",
                      "Observed F/M only (no allocation)"="observed"),
            selected="corrected",inline=TRUE),
          uiOutput("lfd_allocation_status"),
          tabsetPanel(id="structure_subtabs",type="pills",
          tabPanel(
          "Mean length",
          plotOutput("mean_length", height = "500px"),
          br(),
          downloadButton("dl_meanlen_png", "Download PNG"),
          downloadButton("dl_meanlen_csv", "Download CSV")
        ),
          tabPanel(
          "LFD total (all sexes)",
          plotOutput("lfd_blocks_total", height = "850px"),
          br(),
          downloadButton("dl_lfd_tot_png", "Download PNG"),
          downloadButton("dl_lfd_tot_csv", "Download CSV")
        ),
          tabPanel(
          "LFD females",
          p("Corrected frequencies allocate undetermined individuals using the single pooled LOESS sex-ratio curve. Select Observed above to view source F only."),
          plotOutput("lfd_blocks_female", height = "850px"),
          br(),
          downloadButton("dl_lfd_f_png", "Download PNG"),
          downloadButton("dl_lfd_f_csv", "Download CSV")
        ),
          tabPanel(
          "LFD males",
          p("Corrected frequencies allocate undetermined individuals using the same pooled LOESS curve. Select Observed above to view source M only."),
          plotOutput("lfd_blocks_male", height = "850px"),
          br(),
          downloadButton("dl_lfd_m_png", "Download PNG"),
          downloadButton("dl_lfd_m_csv", "Download CSV")
        ),
          tabPanel(
          "LFD undetermined",
          p("Observed unassigned standardised frequencies before sex allocation. This LFD is retained unchanged for auditing."),
          plotOutput("lfd_undet_plot", height = "650px"),
          downloadButton("dl_lfd_u_csv","Download CSV"),
          downloadButton("dl_lfd_u_png","Download PNG")
        ),
          tabPanel(
          "Age–length keys (ALK)",
          tags$div(class="medits-info-note",
             "Upload an optional age-readings CSV in the sidebar (template included in the ZIP). " ,
             "ALKs are calculated from aged observations; they are NOT standardised survey LFDs and are never used to alter the MEDITS model."),
          uiOutput("alk_data_status"),
          selectInput("alk_sex", "ALK group", choices=c("All sexes"="ALL", "Female"="F", "Male"="M"),
                      selected="ALL"),
          uiOutput("alk_year_ui"),
          plotOutput("alk_plot",height="650px"),
          downloadButton("dl_alk_png","Download ALK plot (PNG)"),
          downloadButton("dl_alk_csv","Download ALK probabilities (CSV)"),
          tags$p("ALK: P(Age | Length); probabilities sum to 1 within each length class. The sample count at each length is included in the CSV.")),
          tabPanel(
          "Sample size",
          h4("Total sample size per year (all sexes)"),
          plotOutput("nsample_total",  height = "300px"),
          downloadButton("dl_nsample_tot_png", "PNG"),
          downloadButton("dl_nsample_tot_csv", "CSV"),
          br(), br(),
          
          h4("Female sample size per year"),
          plotOutput("nsample_female", height = "300px"),
          downloadButton("dl_nsample_f_png", "PNG"),
          downloadButton("dl_nsample_f_csv", "CSV"),
          br(), br(),
          
          h4("Male sample size per year"),
          plotOutput("nsample_male",   height = "300px"),
          downloadButton("dl_nsample_m_png", "PNG"),
          downloadButton("dl_nsample_m_csv", "CSV")
        ),
          tabPanel(
          "Sex ratio",
          plotOutput("sex_ratio", height = "800px"),
          br(),
          downloadButton("dl_sr_png", "Download PNG"),
          downloadButton("dl_sr_csv", "Download CSV")
        )
          )
        )
      )
    )
  )
)    

# Review and correction UX follows the commercial CampBiol application:
# manual thresholds, single-record confirmation, searchable audit, no automatic deletion.
medits_bio_review_modal <- function() {
  shiny::modalDialog(
    title=tags$div(class="medits-modal-heading",icon("search"),
                   " Biological review · length, weight and age"),
    size="l", easyClose=TRUE,
    footer=tagList(downloadButton("qc_changes_download","Export review log (CSV)"),
                   modalButton("Close")),
    tags$p(class="medits-page-subtitle",
      "Review is restricted to the selected species, GSA(s), and country/countries. Thresholds are manual; no record is automatically removed or modified."),
    uiOutput("qc_kpis"),
    tabsetPanel(id="bio_review_tabs",
      tabPanel("1 · Screening",
        tags$div(class="medits-row",
          numericInput("qc_length_max","Maximum length (original TC units)",value=1500,min=1),
          numericInput("qc_age_max","Maximum age (years)",value=20,min=0),
          numericInput("qc_weight_max","Maximum weight (original TB/TC units)",value=1000000,min=1)),
        tags$div(class="medits-row",
          selectInput("qc_type","Variable",choices=c("All"="ALL","Length"="LENGTH",
            "Weight"="WEIGHT","Age"="AGE"),selected="ALL"),
          checkboxInput("qc_include_missing","Include missing or non-numeric values",value=FALSE)),
        actionButton("qc_scan","Scan observations",class="btn-primary"),
        verbatimTextOutput("qc_status"),
        tableOutput("qc_fields_table")),
      tabPanel("2 · Review and correct",
        tags$p(class="medits-page-subtitle","Search a flagged observation, inspect the original record and confirm a single targeted correction."),
        textInput("qc_search","Search by GSA, year, haul, species or field",value=""),
        tableOutput("qc_candidates_table"),
        selectizeInput("qc_candidate_id","Record to review",choices=NULL,
          options=list(maxOptions=100,placeholder="Search a value or haul"),width="100%"),
        tableOutput("qc_record_detail"),
        numericInput("qc_replacement","Corrected value (original units)",value=NA_real_),
        tags$div(class="medits-action-row",
          actionButton("qc_apply","Apply correction",class="btn-primary"),
          actionButton("qc_undo_selected","Restore selected record"))),
      tabPanel("3 · Log and exports",
        tableOutput("qc_changes_table"),
        tags$div(class="medits-action-row",
          downloadButton("qc_corrected_download","Download corrected inputs (ZIP)"),
          actionButton("qc_undo","Clear all corrections")),
        fileInput("qc_import","Import previously saved corrections (CSV)",accept=".csv"),
        tags$p(class="medits-page-subtitle",
          "Corrections are applied only to the in-memory copy and to the next Run. If age is supplied separately, it is not used in the LFD calculations."))
    )
  )
}

# GSA comparison is built from the actual, independently estimated GSA runs.
# Aggregated stock indices are not a simple sum/mean of individual GSA series.
# The selected mean-weight table uses TB weight/number for the chosen domain.
# It is an individual-weight statistic, NOT biomass-index / density-index.
medits_mean_weight_view <- function(a, domain='total') {
  domain <- if(domain %in% c('total','shelf','slope')) domain else 'total'
  d <- a$Mweight_domain
  col <- paste0(domain,'_mean_weight')
  ncol <- paste0(domain,'_N_individuals')
  hcol <- paste0(domain,'_N_hauls')
  if(!is.data.frame(d) || !all(c('Year',col,ncol,hcol) %in% names(d))) {
    return(data.frame(Year=numeric(),`Mean_weight(g)`=numeric(),
                      N_individuals=numeric(),N_hauls=integer(),
                      Domain=character(),check.names=FALSE))
  }
  x <- data.frame(Year=as.numeric(as.character(d$Year)),
                  `Mean_weight(g)`=as.numeric(as.character(d[[col]])),
                  N_individuals=as.numeric(as.character(d[[ncol]])),
                  N_hauls=as.integer(d[[hcol]]),
                  Domain=domain,check.names=FALSE)
  x[order(x$Year),,drop=FALSE]
}

medits_collect_gsa_trends <- function(runs, domain='total') {
  if(is.null(runs) || !length(runs)) return(data.frame())
  bc <- switch(domain,total='total_biomass',shelf='shelf_biomass',slope='slope_biomass')
  dc <- switch(domain,total='total_density',shelf='shelf_density',slope='slope_density')
  pieces <- list()
  add <- function(x, col, indicator, gsa, dom, yearcol) {
    if(!is.data.frame(x) || !nrow(x) || !all(c(col,yearcol) %in% names(x))) return(NULL)
    yr <- suppressWarnings(as.numeric(as.character(x[[yearcol]])))
    val <- suppressWarnings(as.numeric(as.character(x[[col]])))
    valid <- is.finite(yr) & is.finite(val)
    if(!any(valid)) return(NULL)
    data.frame(GSA=rep(gsa,sum(valid)),Year=yr[valid],
               Indicator=rep(indicator,sum(valid)),Domain=rep(dom,sum(valid)),
               Value=val[valid],stringsAsFactors=FALSE)
  }
  for(nm in names(runs)) {
    a <- runs[[nm]]
    tag <- if(identical(nm,'AGGREGATED')) 'AGGREGATED' else sub('^GSA_','',nm)
    pieces[[length(pieces)+1L]] <- add(a$BIOMASS,bc,'Biomass (kg/km2)',tag,domain,'year')
    pieces[[length(pieces)+1L]] <- add(a$DENSITY,dc,'Density (n/km2)',tag,domain,'year')
    pieces[[length(pieces)+1L]] <- add(medits_mean_weight_view(a,domain),'Mean_weight(g)','Mean weight (g)',tag,domain,'Year')
  }
  pieces <- Filter(Negate(is.null),pieces)
  if(!length(pieces)) return(data.frame())
  do.call(rbind,pieces)
}

server <- function(input, output, session){
  
  observeEvent(input$open_bio_qc, {
    showModal(medits_bio_review_modal())
  }, ignoreInit=TRUE)

  ## helper: etichetta GSA combinata (es. "16+19")
  gsa_label <- function(gsa_vec) {
    gsa_vec <- unique(as.character(gsa_vec))
    if (length(gsa_vec) > 1) {
      paste(gsa_vec, collapse = "+")
    } else {
      gsa_vec
    }
  }
  
  ## helper: macro stratum selezionato nei plot
  get_view <- function() {
    v <- input$stratum_view
    if (is.null(v)) "total" else v
  }
  
  # SELECTOR: Total / Shelf / Slope (all standardized LFD objects)
  select_lfd <- function(a, view, sex='C') {
    key <- switch(sex, C='LFD',F='LFDF',M='LFDM',U='LFDU')
    if (view %in% c('shelf','slope')) key <- paste0(key,'_',view)
    d <- a[[key]]
    if (!is.data.frame(d)) stop(paste('Missing LFD output:', key))
    d
  }
  select_sr <- function(a,view) {
    key <- if(view %in% c('shelf','slope')) paste0('LFDSR_',view) else 'LFDSR'
    d <- a[[key]]
    if (!is.data.frame(d)) stop(paste('Missing sex ratio output:',key))
    d
  }

  ## helper: prende una colonna cercando sia MAIUSCOLO che minuscolo
  get_col <- function(df, upper, lower) {
    if (upper %in% names(df)) {
      df[[upper]]
    } else if (lower %in% names(df)) {
      df[[lower]]
    } else {
      NULL
    }
  }
  
  ## --- quick read of TA to populate GSA / country selectors ----------
  ta_head <- reactive({
    req(input$ta_file)
    fread(input$ta_file$datapath)
  })
  
  output$gsa_ui <- renderUI({
    ta <- ta_head()
    areas <- sort(unique(get_col(ta, "AREA", "area")))
    selectInput(
      "gsa", "GSA (AREA)",
      choices  = areas,
      selected = areas,
      multiple = TRUE
    )
  })
  
  output$country_ui <- renderUI({
    ta <- ta_head()
    cnt <- sort(unique(get_col(ta, "COUNTRY", "country")))
    selectInput(
      "country", "Country",
      choices  = cnt,
      selected = cnt,
      multiple = TRUE
    )
  })
  
  ## ----- MEDITS record-level QC: lengths, weights, optional ages ----------
  raw_input <- reactive({
    req(input$ta_file,input$tb_file,input$tc_file)
    files <- list(TA=input$ta_file$datapath,TB=input$tb_file$datapath,
                  TC=input$tc_file$datapath)
    if (!is.null(input$age_file) && nrow(input$age_file))
      files$AGE <- input$age_file$datapath
    lapply(files,function(path) as.data.frame(data.table::fread(path)))
  })
  qc_empty_register <- function() data.frame(
    Source=character(),Row=integer(),Field=character(),Type=character(),
    Original=numeric(),Corrected=numeric(),Country=character(),GSA=character(),
    Year=character(),Haul=character(),Species=character(),
    stringsAsFactors=FALSE)
  corrections <- reactiveVal(qc_empty_register())
  observeEvent(list(input$ta_file,input$tb_file,input$tc_file,input$age_file),{
    corrections(qc_empty_register())
  },ignoreInit=TRUE)
  qc_data <- reactive({
    a <- raw_input()
    changes <- corrections()
    for (i in seq_len(nrow(changes))) {
      c <- changes[i,,drop=FALSE]
      if (c$Source %in% names(a) && c$Field %in% names(a[[c$Source]]) &&
          c$Row>=1L && c$Row<=nrow(a[[c$Source]]))
        a[[c$Source]][c$Row,c$Field] <- c$Corrected
    }
    a
  })
  # QC is species- and geography-scoped, independently of the selected figure.
  # Original row positions are kept, so applying a correction remains safe.
  qc_selected_rows <- reactive({
    req(input$species,input$gsa,input$country)
    medits_qc_species_rows(qc_data(),alpha=input$species,
       gsa=input$gsa,country=input$country,
       species_file=file.path(app_dir,'Sp_Medits.csv'),
       asfis_file=file.path(app_dir,'ASFIS_2020.csv'))
  })
  qc_all_results <- reactive({
    req(input$qc_scan > 0)
    medits_qc_scan_data(qc_data(),max_length=input$qc_length_max,
      max_weight=input$qc_weight_max,max_age=input$qc_age_max,
      include_missing=isTRUE(input$qc_include_missing),selected_type='ALL',
      selected_rows=qc_selected_rows())
  })
  qc_candidates <- reactive({
    z <- qc_all_results()$candidates
    if (!identical(input$qc_type,'ALL')) z <- z[z$Type==input$qc_type,,drop=FALSE]
    z
  })
  qc_visible <- reactive(medits_qc_subset_candidates(qc_candidates(),input$qc_search))
  qc_kpi_card <- function(value,label) tags$div(class='medits-kpi',
                        tags$div(class='value',as.character(value)),
                        tags$div(class='label',label))
  output$qc_kpis <- renderUI({
    n_corrected <- nrow(corrections())
    if (is.null(input$qc_scan) || input$qc_scan==0)
      return(tags$div(class='medits-kpis',
        qc_kpi_card('—','Length flags'),qc_kpi_card('—','Weight flags'),
        qc_kpi_card('—','Age flags'),qc_kpi_card(n_corrected,'Saved corrections')))
    z <- qc_all_results()$candidates
    count <- function(typ) sum(z$Type==typ)
    tags$div(class='medits-kpis',
      qc_kpi_card(format(count('LENGTH'),big.mark=','),'Length flags'),
      qc_kpi_card(format(count('WEIGHT'),big.mark=','),'Weight flags'),
      qc_kpi_card(format(count('AGE'),big.mark=','),'Age flags'),
      qc_kpi_card(n_corrected,'Saved corrections'))
  })
  output$qc_status <- renderText({
    if (is.null(input$qc_scan) || input$qc_scan==0)
      return('No scan yet. Upload the source files and click Scan biological records.')
    x <- qc_all_results()
    nflag <- nrow(x$candidates)
    n_age <- sum(x$catalog$Type=='AGE')
    paste0('Flagged values: ',nflag,
      ' | Visible after filters: ',nrow(qc_visible()),
      ' | Saved manual corrections: ',nrow(corrections()),'.\n',
      'Selected species: ',toupper(trimws(input$species)),
      ' | GSA: ',paste(input$gsa,collapse=','),
      ' | Countries: ',paste(input$country,collapse=','),'.\n',
      'Screened input fields: ',nrow(x$catalog),'. ',
      if (n_age) paste0('Age fields available: ',n_age,'.') else
        'No age observations found. Standard MEDITS TA/TB/TC do not contain individual age; upload an age-readings CSV with an AGE column if available.',
      '\nLimits are review rules only. Zero flags does not mean zero records: check the observed ranges in the field summary. No automatic removal.',
      if(length(attr(qc_selected_rows(),'notices')))
        paste0('\n',paste(attr(qc_selected_rows(),'notices'),collapse='; ')) else '')
  })
  observeEvent(input$qc_scan,{
    x <- qc_all_results()
    showNotification(paste('QC scan:',nrow(x$candidates),'flagged values across',
                           nrow(x$catalog),'fields.'),
                     type=if(nrow(x$candidates)) 'warning' else 'message',duration=6)
  },ignoreInit=TRUE)
  output$qc_fields_table <- renderTable({
    x <- qc_all_results()$fields
    if(!nrow(x)) return(data.frame(Status='No supported biological fields found'))
    x
  },striped=TRUE,spacing='xs',digits=4)
  output$qc_candidates_table <- renderTable({
    x <- qc_visible()
    if(!nrow(x)) return(data.frame(Status='No flagged records match the current filters'))
    head(x,500)
  },striped=TRUE,spacing='xs',digits=5)
  # Server-side searchable selectize: do not send thousands of options to browser.
  observeEvent(qc_visible(), {
    x <- qc_visible()
    if (!nrow(x)) {
      shiny::updateSelectizeInput(session,'qc_candidate_id',choices=character(),selected=character(),server=TRUE)
      return()
    }
    labels <- paste(x$Source,'|',x$Type,'| GSA',x$GSA,'|',x$Year,
                    '| Haul',x$Haul,'|',x$Field,'=',x$Value,
                    '| Source row',x$Row)
    choices <- setNames(as.character(x$ID),labels)
    previous <- isolate(input$qc_candidate_id)
    selected <- if (!is.null(previous) && previous %in% x$ID) previous else as.character(x$ID[[1]])
    shiny::updateSelectizeInput(session,'qc_candidate_id',choices=choices,
                                selected=selected,server=TRUE)
  },ignoreInit=FALSE)
  qc_selected <- reactive({
    x <- qc_candidates(); req(input$qc_candidate_id)
    q <- x[x$ID==suppressWarnings(as.integer(input$qc_candidate_id)),,drop=FALSE]
    req(nrow(q)==1L)
    q
  })
  output$qc_record_detail <- renderTable({
    q <- qc_selected(); src <- as.character(q$Source); row <- as.integer(q$Row)
    original <- raw_input()[[src]]; current <- qc_data()[[src]]
    if(is.null(original) || row>nrow(original)) return(data.frame())
    fields <- c('COUNTRY','AREA','GSA','YEAR','HAUL_NUMBER','GENUS','SPECIES',
                'SEX','LENGTH_CLASS','AGE','AGE_YEARS','AGE_CLASS','AGE_GROUP',
                'TOTAL_WEIGHT_IN_THE_HAUL','WEIGHT_OF_THE_SAMPLE_MEASURED',
                'WEIGHT_OF_THE_FRACTION',as.character(q$Field))
    cols <- which(medits_qc_norm(names(original)) %in% medits_qc_norm(fields))
    cols <- unique(c(cols,match(as.character(q$Field),names(original))))
    cols <- cols[!is.na(cols)]
    data.frame(Field=names(original)[cols],Original=as.character(unlist(original[row,cols,drop=FALSE],use.names=FALSE)),
               Current=as.character(unlist(current[row,cols,drop=FALSE],use.names=FALSE)),
               stringsAsFactors=FALSE)
  },striped=TRUE,spacing='xs')
  observeEvent(input$qc_candidate_id,{
    q <- tryCatch(qc_selected(),error=function(e) NULL)
    if(!is.null(q) && nrow(q))
      updateNumericInput(session,'qc_replacement',value=if(is.finite(q$Value)) q$Value else NA_real_)
  },ignoreInit=TRUE)
  observeEvent(input$qc_apply,{
    q <- tryCatch(qc_selected(),error=function(e) NULL)
    if(is.null(q) || !nrow(q)) {
      showNotification('Select a flagged observation first.',type='warning'); return()
    }
    value <- suppressWarnings(as.numeric(input$qc_replacement))
    if(!length(value) || !is.finite(value)) {
      showNotification('Enter a finite numeric correction.',type='error'); return()
    }
    if ((q$Type=='LENGTH' && value<=0) ||
        (q$Type=='WEIGHT' && value<0) || (q$Type=='AGE' && value<0)) {
      showNotification('Length must be positive; weight and age cannot be negative.',type='error'); return()
    }
    src <- as.character(q$Source); row <- as.integer(q$Row); field <- as.character(q$Field)
    original <- raw_input()[[src]][[field]][row]
    original <- suppressWarnings(as.numeric(as.character(original)))
    prev <- corrections()
    key <- prev$Source==src & prev$Row==row & prev$Field==field
    if(any(key)) prev <- prev[!key,,drop=FALSE]
    # A correction identical to the raw value is not stored.
    if (!is.finite(original) || !isTRUE(all.equal(value,original))) {
      prev <- rbind(prev,data.frame(Source=src,Row=row,Field=field,Type=as.character(q$Type),
        Original=original,Corrected=value,Country=as.character(q$Country),
        GSA=as.character(q$GSA),Year=as.character(q$Year),
        Haul=as.character(q$Haul),Species=as.character(q$Species),
        stringsAsFactors=FALSE))
    }
    corrections(prev)
    showNotification(paste('Correction registered:',src,'row',row,field,
                           '→',value,'. Rerun the analysis to refresh results.'),type='message',duration=7)
  })
  observeEvent(input$qc_undo,{
    corrections(qc_empty_register())
    showNotification('All in-memory QC corrections cleared. Raw uploaded files unchanged.',type='message')
  })
  # Importare correzioni è sicuro solo se ogni riga corrisponde alla sorgente
  # effettivamente caricata (campo, numero di riga, vecchio valore e chiave).
  observeEvent(input$qc_import, {
    req(input$qc_import$datapath)
    incoming <- tryCatch(utils::read.csv(input$qc_import$datapath,
      stringsAsFactors=FALSE,check.names=FALSE),error=function(e)e)
    if(inherits(incoming,"error")) {
      showNotification("Unable to read the correction log CSV",type="error");return()
    }
    needed <- names(qc_empty_register())
    if(!all(needed %in% names(incoming))) {
      showNotification("Incompatible log: required fields are missing",type="error");return()
    }
    incoming <- incoming[,needed,drop=FALSE]
    raw <- raw_input()
    valid <- vapply(seq_len(nrow(incoming)),function(i){
      r <- incoming[i,,drop=FALSE]
      src <- as.character(r$Source);field <- as.character(r$Field)
      ix <- suppressWarnings(as.integer(r$Row))
      if(length(ix)!=1L || !is.finite(ix) || !src %in% names(raw) ||
         !field %in% names(raw[[src]]) || ix<1L || ix>nrow(raw[[src]])) return(FALSE)
      orig <- suppressWarnings(as.numeric(as.character(raw[[src]][[field]][ix])))
      expected <- suppressWarnings(as.numeric(as.character(r$Original)))
      revised <- suppressWarnings(as.numeric(as.character(r$Corrected)))
      if(!is.finite(revised)) return(FALSE)
      if(is.finite(orig) != is.finite(expected)) return(FALSE)
      if(is.finite(orig) && !isTRUE(all.equal(orig,expected,tolerance=0))) return(FALSE)
      for(z in list(c("COUNTRY","Country"),c("AREA","GSA"),
                    c("YEAR","Year"),c("HAUL_NUMBER","Haul"))){
        field_raw <- medits_qc_col(raw[[src]],z[1])
        if(length(field_raw)>=ix) {
          a <- trimws(as.character(field_raw[ix])); b <- trimws(as.character(r[[z[2]]]))
          if(!is.na(a) && nzchar(a) && !is.na(b) && nzchar(b) && a!=b) return(FALSE)
        }
      }
      TRUE
    },logical(1))
    if(!length(valid) || any(!valid)) {
      showNotification(paste0("Import cancelled: ",sum(!valid),
         " corrections do not match the current inputs."),type="error",duration=10)
      return()
    }
    current <- corrections()
    both <- rbind(current,incoming)
    key <- paste(both$Source,both$Row,both$Field,sep="::")
    both <- both[!duplicated(key,fromLast=TRUE),,drop=FALSE]
    corrections(both)
    showNotification(paste(nrow(incoming),"corrections restored and validated."),type="message")
  },ignoreInit=TRUE)
  observeEvent(input$qc_undo_selected, {
    q <- tryCatch(qc_selected(),error=function(e)NULL)
    if(is.null(q) || !nrow(q)) return()
    prev <- corrections()
    keep <- !(prev$Source==q$Source & prev$Row==q$Row & prev$Field==q$Field)
    corrections(prev[keep,,drop=FALSE])
    showNotification("Original value restored for the selected record.",type="message")
  })
  output$qc_changes_table <- renderTable({
    x <- corrections()
    if(!nrow(x)) return(data.frame(Status='No manual corrections saved'))
    x
  },striped=TRUE,spacing='xs')
  output$qc_changes_download <- downloadHandler(
    filename=function() 'MEDITS_QC_CORRECTIONS.csv',
    content=function(file) utils::write.csv(corrections(),file,row.names=FALSE))
  output$qc_corrected_download <- downloadHandler(
    filename=function() paste0('MEDITS_CORRECTED_TA_TB_TC_AGE_',Sys.Date(),'.zip'),
    contentType='application/zip',
    content=function(file) {
      a <- qc_data()
      stage <- tempfile('medits_qc_csv_')
      dir.create(stage,recursive=TRUE)
      on.exit(unlink(stage,recursive=TRUE,force=TRUE),add=TRUE)
      for(src in names(a)) data.table::fwrite(a[[src]],file.path(stage,paste0(src,'_corrected.csv')))
      utils::write.csv(corrections(),file.path(stage,'MEDITS_QC_CORRECTIONS.csv'),row.names=FALSE)
      if(isTRUE(input$qc_scan>0)) {
        utils::write.csv(qc_all_results()$fields,file.path(stage,'MEDITS_QC_FIELD_SUMMARY.csv'),row.names=FALSE)
        utils::write.csv(qc_all_results()$candidates,file.path(stage,'MEDITS_QC_FLAGGED_RECORDS.csv'),row.names=FALSE)
      }
      writeLines(c('MEDITS corrected inputs, record-level audit',
        'Source files TA/TB/TC are copies with approved manual corrections applied.',
        'The optional AGE CSV is included only if it was uploaded.',
        'AGE readings are quality-controlled and exported but are NOT passed to the MEDITS LFD algorithms.',
        'No changes are made to original uploaded source files.'),
        file.path(stage,'README_CORRECTED_INPUTS.txt'))
      if(!requireNamespace('zip',quietly=TRUE)) stop("Install the 'zip' package to export source data: install.packages('zip')")
      zip::zipr(zipfile=file,files=list.files(stage,full.names=FALSE),root=stage,
                include_directories=FALSE)
    })

  output$display_gsa_ui <- renderUI({
    gsas <- input$gsa
    if (!length(gsas)) return(NULL)
    selectInput('display_gsa','Geographic view',
                choices=c(setNames('AGGREGATED',paste0('Aggregated (',paste(gsas,collapse='+'),')')),
                          setNames(paste0('GSA_',gsas),paste('GSA',gsas))),
                selected='AGGREGATED')
  })

  medits_run_state <- reactiveVal("Ready. Upload TA, TB and TC, then click Run MEDITS analysis.")
  output$medits_run_status <- renderText(medits_run_state())

  ## --- full analysis reactive (analisi sempre su tutti gli strati) --
  all_analyses <- eventReactive(input$run_btn, {
    req(input$ta_file, input$tb_file, input$tc_file)
    medits_run_state("Analysis running: computing GSA results and MEDITS indices...")
    tryCatch(shiny::withProgress(message="MEDITS analysis", detail="Preparing input records", value=0, {
    
    inspected <- qc_data()
    TAraw <- inspected$TA
    TBraw <- inspected$TB
    TCraw <- inspected$TC
    
    ## 1) Costruisco un ID di haul comune a TA/TB/TC
    TA_country <- get_col(TAraw, "COUNTRY", "country")
    TA_area    <- get_col(TAraw, "AREA",    "area")
    TA_year    <- get_col(TAraw, "YEAR",    "year")
    TA_haul    <- get_col(TAraw, "HAUL_NUMBER", "haul_number")
    TAraw$id_haul <- paste(TA_country, TA_area, TA_year, TA_haul, sep = "_")
    
    TB_country <- get_col(TBraw, "COUNTRY", "country")
    TB_area    <- get_col(TBraw, "AREA",    "area")
    TB_year    <- get_col(TBraw, "YEAR",    "year")
    TB_haul    <- get_col(TBraw, "HAUL_NUMBER", "haul_number")
    TBraw$id_haul <- paste(TB_country, TB_area, TB_year, TB_haul, sep = "_")
    
    TC_country <- get_col(TCraw, "COUNTRY", "country")
    TC_area    <- get_col(TCraw, "AREA",    "area")
    TC_year    <- get_col(TCraw, "YEAR",    "year")
    TC_haul    <- get_col(TCraw, "HAUL_NUMBER", "haul_number")
    TCraw$id_haul <- paste(TC_country, TC_area, TC_year, TC_haul, sep = "_")
    
    ## 2) Profondità media nel TA
    shooting_depth <- get_col(TAraw, "SHOOTING_DEPTH", "shooting_depth")
    hauling_depth  <- get_col(TAraw, "HAULING_DEPTH",  "hauling_depth")
    
    if (!is.null(shooting_depth) && !is.null(hauling_depth)) {
      TAraw$MEAN_DEPTH <- rowMeans(
        cbind(shooting_depth, hauling_depth),
        na.rm = TRUE
      )
    } else if (!is.null(hauling_depth)) {
      TAraw$MEAN_DEPTH <- hauling_depth
    } else {
      TAraw$MEAN_DEPTH <- NA_real_
    }
    
    ## 3) Macro-strata in base alla profondità media:
    ##    - shelf: 1–200 m
    ##    - slope: 201 m – max
    TAraw$MACRO_STRATA <- NA_character_
    
    TAraw$MACRO_STRATA[
      !is.na(TAraw$MEAN_DEPTH) &
        TAraw$MEAN_DEPTH >= 1 &
        TAraw$MEAN_DEPTH <= 200
    ] <- "shelf"
    
    TAraw$MACRO_STRATA[
      !is.na(TAraw$MEAN_DEPTH) &
        TAraw$MEAN_DEPTH >= 201
    ] <- "slope"
    
    ## Copio HAULING_DEPTH anche su TB e TC se ho la colonna
    if (!is.null(hauling_depth)) {
      TBraw$HAULING_DEPTH <- hauling_depth[match(TBraw$id_haul, TAraw$id_haul)]
      TCraw$HAULING_DEPTH <- hauling_depth[match(TCraw$id_haul, TAraw$id_haul)]
    }
    
    ## Parametri per run_medits_from_raw
    areacode    <- unique(as.character(unlist(input$gsa)))
    countrycode <- unique(as.character(unlist(input$country)))
    sspp        <- toupper(as.character(input$species))  # HKE, ARS, ecc.
    namesurvey  <- as.character(input$survey)
    
    if (!length(areacode) || !length(countrycode))
      stop("Select at least one GSA and one country before Run analysis")
    base_out <- file.path(app_dir, paste0(namesurvey,"_output"))
    dir.create(base_out, recursive=TRUE, showWarnings=FALSE)
    shiny::incProgress(0.08,detail="Computing aggregate GSA results")
    ## 5) Lancio l’analisi usando TUTTI gli haul (tutti gli strati)
    res <- run_medits_from_raw(
      TAraw       = TAraw,
      TBraw       = TBraw,
      TCraw       = TCraw,
      areacode    = areacode,
      countrycode = countrycode,
      sspp        = sspp,
      namesurvey  = namesurvey,
      data_format = "DCF",
      output_base_dir = base_out
    )
    
    
    results <- list(AGGREGATED=res)
    gsa_status <- data.frame(GSA=as.character(areacode),Status="pending",Details="",
                             stringsAsFactors=FALSE)
    if (length(areacode)>1L) {
      for (j in seq_along(areacode)) {
        one_gsa <- areacode[[j]]
        shiny::incProgress(0.80/length(areacode),
          detail=paste0("Individual GSA ",one_gsa," (",j,"/",length(areacode),")"))
        message("Running individual GSA ",one_gsa)
        single <- tryCatch(run_medits_from_raw(
          TAraw=TAraw, TBraw=TBraw, TCraw=TCraw,
          areacode=one_gsa, countrycode=countrycode,
          sspp=sspp, namesurvey=namesurvey, data_format="DCF",
          global_sex_ratio=res$SR_vec,
          output_base_dir=base_out
        ), error=function(e) e)
        if (inherits(single,"error")) {
          gsa_status$Status[j] <- "ERROR"
          gsa_status$Details[j] <- conditionMessage(single)
          warning("Individual GSA ",one_gsa,": ",conditionMessage(single))
        } else {
          gsa_status$Status[j] <- "OK"
          gsa_status$Details[j] <- "Individual run computed"
          results[[paste0("GSA_",one_gsa)]] <- single
        }
      }
    } else {
      # With a single selected GSA, aggregate and individual results are identical.
      results[[paste0("GSA_",areacode[1])]] <- res
      gsa_status$Status[] <- "OK"
      gsa_status$Details[] <- "Single selected GSA; no duplicate calculation"
    }
    if (!is.null(res$output_run_dir) && dir.exists(res$output_run_dir)) {
      write.csv(corrections(),file.path(res$output_run_dir,'QC_MANUAL_CORRECTIONS.csv'),row.names=FALSE)
      write.csv(gsa_status,file.path(res$output_run_dir,'GSA_RUN_STATUS.csv'),row.names=FALSE)
    }
    for (nm in names(results)) {
      r <- results[[nm]]
      r$report_species <- sspp
      r$report_gsas <- if(nm=='AGGREGATED') areacode else sub('GSA_','',nm)
      r$report_years <- sort(unique(as.integer(r$TAn$year)))
      r$report_survey <- namesurvey
      results[[nm]] <- r
    }
    attr(results,"gsa_status") <- gsa_status
    failures <- sum(gsa_status$Status != "OK")
    medits_run_state(paste0("Analysis completed · Aggregated + ",sum(gsa_status$Status=="OK"),
        " GSA elaborate; ",failures," errori.",
        if(failures) " Consulta GSA_RUN_STATUS.csv nello ZIP." else ""))
    shiny::incProgress(0.12,detail="Ready: CSV/JPEG exports")
    results
    }), error=function(e) {
      msg <- paste("Run analysis ERROR:", conditionMessage(e))
      medits_run_state(msg)
      shiny::showNotification(msg, type="error", duration=NULL)
      stop(e)
    })  # end withProgress/tryCatch
  })

  # IMPORTANT: eventReactive is lazy. Without a consumer in the active QC tab,
  # pressing Run would otherwise not execute anything until another tab opened.
  observeEvent(input$run_btn, {
    tryCatch(invisible(all_analyses()), error=function(e) {
      message("MEDITS run failed: ",conditionMessage(e))
    })
  }, priority=-100, ignoreInit=TRUE)

  analysis <- reactive({
    results <- all_analyses()
    key <- input$display_gsa
    if (is.null(key) || !key %in% names(results)) key <- 'AGGREGATED'
    results[[key]]
  })

  # All GSA lines in a single three-panel comparison; the old per-GSA
  # selector and its independent trends are retained unchanged.
  gsa_comparison_data <- reactive({
    medits_collect_gsa_trends(all_analyses(),get_view())
  })
  gsa_comparison_plot <- reactive({
    d <- gsa_comparison_data()
    validate(need(nrow(d)>0,'No GSA trend series available for this domain. Check GSA_RUN_STATUS.'))
    ggplot2::ggplot(d,ggplot2::aes(x=Year,y=Value,colour=GSA,group=GSA))+
      ggplot2::geom_line(linewidth=1,na.rm=TRUE)+
      ggplot2::geom_point(size=1.9,na.rm=TRUE)+
      ggplot2::facet_wrap(~Indicator,scales='free_y',ncol=1)+
      ggplot2::scale_x_continuous(breaks=scales::pretty_breaks(10))+
      ggplot2::labs(title=paste0('MEDITS: aggregate versus individual GSAs (',get_view(),')'),
                    x='Year',y=NULL,colour='Series')+
      ggplot2::theme_minimal(base_size=13)+
      ggplot2::theme(plot.title=ggplot2::element_text(face='bold',size=16),
                     legend.position='bottom',strip.text=ggplot2::element_text(face='bold',size=13))
  })
  output$gsa_run_table <- shiny::renderTable({
    st <- attr(all_analyses(), 'gsa_status')
    if (!is.data.frame(st)) return(data.frame(Status='No individual-GSA run status available'))
    st
  },striped=TRUE,spacing='xs')
  output$gsa_comparison_plot <- shiny::renderPlot(gsa_comparison_plot())
  output$dl_gsa_comparison_csv <- shiny::downloadHandler(
    filename=function() paste0('MEDITS_GSA_trends_',get_view(),'.csv'),
    content=function(file) utils::write.csv(gsa_comparison_data(),file,row.names=FALSE))
  output$dl_gsa_comparison_png <- shiny::downloadHandler(
    filename=function() paste0('MEDITS_GSA_comparison_',get_view(),'.png'),
    content=function(file) ggplot2::ggsave(file,gsa_comparison_plot(),width=13,height=11,dpi=300))

  ## --- year selector for haul map ------------------------------------
  output$year_map_ui <- renderUI({
    a <- analysis()
    years <- sort(unique(a$TAn$year))
    
    sliderInput(
      "map_year",
      label   = "Year to display",
      min     = min(years),
      max     = max(years),
      value   = min(years),
      step    = 1,
      animate = animationOptions(interval = 1200, loop = TRUE)
    )
  })
  
  ## Haul map: same coordinate limits in all facets; coast segments independent.
  ## 'n' means unique VALID TA hauls with usable coordinates in displayed depth domain.
  medits_haul_map_data <- function(a, view = "total", all_years = TRUE, map_year = NULL) {
    ta <- medits_survey_unique_hauls(a$TAn)
    validate(need(nrow(ta) > 0, "No valid TA hauls available"))
    required <- c("year", "TA_lon_fin", "TA_lat_fin", "meandepth")
    validate(need(all(required %in% names(ta)),
                  paste("Missing columns:", paste(setdiff(required, names(ta)), collapse=", "))))
    ta$YEAR <- suppressWarnings(as.integer(as.character(ta$year)))
    ta$TA_lon_fin <- suppressWarnings(as.numeric(ta$TA_lon_fin))
    ta$TA_lat_fin <- suppressWarnings(as.numeric(ta$TA_lat_fin))
    ta$meandepth <- suppressWarnings(as.numeric(ta$meandepth))
    ta <- ta[is.finite(ta$YEAR) & is.finite(ta$TA_lon_fin) & is.finite(ta$TA_lat_fin) &
             abs(ta$TA_lon_fin) <= 180 & abs(ta$TA_lat_fin) <= 90 &
             ta$TA_lon_fin != 0 & ta$TA_lat_fin != 0, ,drop=FALSE]
    validate(need(nrow(ta) > 0, "No valid geographical positions found in TA"))
    all_valid <- ta
    if (view == "shelf") ta <- ta[is.finite(ta$meandepth) & ta$meandepth >= 1 & ta$meandepth <= 200,,drop=FALSE]
    if (view == "slope") ta <- ta[is.finite(ta$meandepth) & ta$meandepth > 200,,drop=FALSE]
    if (!isTRUE(all_years)) {
      validate(need(!is.null(map_year), "Select a map year"))
      ta <- ta[ta$YEAR == as.integer(map_year),,drop=FALSE]
    }
    validate(need(nrow(ta) > 0, "No mapped hauls in the selected display stratum/year"))
    # Global limits always originate from all valid positions, ensuring comparability.
    xlim <- range(pretty(all_valid$TA_lon_fin))
    ylim <- range(pretty(all_valid$TA_lat_fin))
    xpad <- max(.08,diff(xlim)*.035); ypad <- max(.08,diff(ylim)*.035)
    xlim <- xlim+c(-xpad,xpad); ylim <- ylim+c(-ypad,ypad)
    years <- if (isTRUE(all_years)) sort(unique(all_valid$YEAR)) else sort(unique(ta$YEAR))
    n_by_year <- as.integer(tabulate(match(ta$YEAR, years), nbins=length(years)))
    year_levels <- sprintf("%d\nn = %d",years,n_by_year)
    year_labels <- setNames(year_levels,as.character(years))
    ta$Year_panel <- factor(year_labels[as.character(ta$YEAR)], levels=year_levels)
    # A dummy panel data frame keeps years without hauls in the selected stratum visible.
    panel <- data.frame(YEAR=years, x=rep(mean(xlim),length(years)),y=rep(mean(ylim),length(years)),
                        Year_panel=factor(year_levels,levels=year_levels))
    list(data=ta,panel=panel,xlim=xlim,ylim=ylim,years=years,counts=n_by_year,labels=year_labels)
  }

  make_medits_haul_plot <- function(a, view="total", all_years=TRUE, map_year=NULL,
                                   report=FALSE) {
    d <- medits_haul_map_data(a,view,all_years,map_year)
    gsa <- gsa_label(a$gsa)
    title <- paste0("Haul positions - GSA",gsa," ",a$gen,a$spec," (",a$state,") [",view,"]")
    p <- ggplot2::ggplot() +
      ggplot2::geom_segment(
        data=medits_coast_segments(d$xlim,d$ylim),
        ggplot2::aes(x=LONG,y=LAT,xend=LONG_END,yend=LAT_END),
        inherit.aes=FALSE,colour="grey50",linewidth=.24) +
      ggplot2::geom_point(
        data=d$data,
        ggplot2::aes(x=TA_lon_fin,y=TA_lat_fin),
        inherit.aes=FALSE,colour="firebrick",
        size=if(isTRUE(report)) 0.40 else 1.0,
        alpha=if(isTRUE(report)) 0.78 else 0.72) +
      ggplot2::geom_blank(data=d$panel,
                          ggplot2::aes(x=x,y=y),inherit.aes=FALSE) +
      ggplot2::coord_quickmap(xlim=d$xlim,ylim=d$ylim,expand=FALSE) +
      ggplot2::labs(title=title,subtitle="n = valid unique TA hauls mapped in the displayed depth domain",
                    x="Longitude",y="Latitude") +
      ggplot2::theme_bw(base_size=12) +
      ggplot2::theme(plot.title=ggplot2::element_text(size=14,face="bold",hjust=.5),
                     plot.subtitle=ggplot2::element_text(size=9,hjust=.5))
    if (isTRUE(report)) {
      p <- p +
        ggplot2::scale_x_continuous(breaks=pretty(d$xlim,n=2)) +
        ggplot2::scale_y_continuous(breaks=pretty(d$ylim,n=2)) +
        ggplot2::theme(
          strip.text=ggplot2::element_text(size=8,face="bold"),
          axis.text=ggplot2::element_text(size=7),
          panel.grid.minor=ggplot2::element_blank(),
          panel.spacing=grid::unit(1.1,"mm"))
    }
    if (isTRUE(all_years)) {
      p <- p + ggplot2::facet_wrap(~Year_panel, ncol=6,drop=FALSE) +
        ggplot2::theme(strip.text=ggplot2::element_text(size=if(isTRUE(report)) 8 else 9,face="bold"),
                       axis.text=ggplot2::element_text(size=if(isTRUE(report)) 7 else 8))
    } else {
      p <- p + ggplot2::labs(subtitle=paste0("Year ",d$years[1],"; n = ",d$counts[1],
                                             " georeferenced valid hauls in the selected depth domain"))
    }
    p
  }

  haul_map_plot <- reactive({
    make_medits_haul_plot(analysis(),get_view(),
                          all_years=isTRUE(input$haul_show_all_years),
                          map_year=isolate(input$map_year))
  })
  # The PDF always contains all survey years, independently of the selected UI year.
  haul_map_report_plot <- reactive({
    make_medits_haul_plot(analysis(),get_view(),all_years=TRUE,report=TRUE)
  })

  survey_period_plot <- reactive({
    a <- analysis()
    medits_make_survey_period_plot(a$TAn,
                                   gsa=gsa_label(a$gsa),state=a$state,
                                   ncol=4)
  })
  # Layout app: 4 colonne e altezza dinamica (INVARIATO).
  # Il PDF usa invece un calendario su 31 righe: un anno per riga.
  survey_period_report_plot <- reactive({
    a <- analysis()
    medits_make_survey_period_report_plot(a$TAn,
                                          gsa=gsa_label(a$gsa),state=a$state)
  })

  output$survey_plot_ui <- renderUI({
    a <- analysis()
    n_years <- length(unique(a$TAn$year[!is.na(a$TAn$year)]))
    n_rows <- max(1L, ceiling(n_years/4L))
    # Four columns, ~180 screen px per annual panel; no compression.
    plotOutput("survey_period",
               height=paste0(max(660L, as.integer(n_rows*190L+155L)),"px"),
               width="100%")
  })


  # Shared map helper: no inherited ggplot aesthetics and no `group` aesthetic.
  # Coast segments, not polygons, prevent the group-length error in ggplot2.
  make_medits_spatial_plot <- function(a, metric, view = "total", report=FALSE) {
    if (!metric %in% c("BIOMASS", "DENSITY")) stop("Unknown map metric")
    df <- medits_spatial_prepare(a,metric,view)
    xlims <- range(df$LONG) + c(-0.3, 0.3)
    ylims <- range(df$LAT) + c(-0.3, 0.3)
    # Guard against zero-width limits when the dataset has a single position.
    if (diff(xlims) <= 0) xlims <- mean(xlims) + c(-0.3, 0.3)
    if (diff(ylims) <= 0) ylims <- mean(ylims) + c(-0.3, 0.3)
    coast <- medits_coast_segments(xlims, ylims)
    df$log_value <- log1p(df[[metric]])
    view_label <- if (view == "total") "all strata" else view
    title <- paste0(if (metric == "BIOMASS") "Biomass" else "Density",
                    " (log1p) - ", a$report_species, " - GSA ",
                    paste(a$report_gsas, collapse = ", "), " [", view_label, "]")
    # The map shows per-haul catch/area, not the area-stratified annual index.
    # The map uses independent datasets for coastlines and hauls.
    # Absolutely no polygon/`group` aesthetics are mapped.
    p <- ggplot2::ggplot() +
      ggplot2::geom_segment(
        data = coast,
        mapping = ggplot2::aes(x = LONG, y = LAT, xend = LONG_END, yend = LAT_END),
        inherit.aes = FALSE, colour = "grey50", linewidth = 0.25
      ) +
      ggplot2::geom_point(
        data = df,
        mapping = ggplot2::aes(x = LONG, y = LAT, fill = log_value),
        inherit.aes = FALSE, shape = 21, colour = "#173749", stroke = 0.24,
        size = if(isTRUE(report)) 1.2 else 2.9,
        alpha = if(isTRUE(report)) 0.96 else 0.98
      ) +
      ggplot2::facet_wrap(~YEAR, ncol = 6) +
      ggplot2::coord_quickmap(xlim = xlims, ylim = ylims, expand = FALSE) +
      ggplot2::scale_fill_gradientn(
        colours = c("#2671ac", "#40c1bd", "#f1d85a", "#f28a43", "#b82b46"),
        name = paste0("log1p(", metric, ")")) +
      ggplot2::labs(title = title, x = "Longitude", y = "Latitude") +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        legend.position = "bottom",
        strip.text = ggplot2::element_text(size=11, face="bold"),
        plot.title = ggplot2::element_text(size=14, face="bold", hjust=0.5)
      )
    if(isTRUE(report)) {
      # Few ticks: avoid overlapped labels in the 31 compact annual panels.
      p <- p +
        ggplot2::scale_x_continuous(breaks=pretty(xlims,n=2)) +
        ggplot2::scale_y_continuous(breaks=pretty(ylims,n=2)) +
        ggplot2::theme(
          strip.text=ggplot2::element_text(size=8,face="bold"),
          axis.text=ggplot2::element_text(size=7),
          panel.spacing=grid::unit(1,"mm"),
          panel.grid.minor=ggplot2::element_blank(),
          legend.key.width=grid::unit(12,"mm"))
    }
    p
  }

  map_export_data <- function(a, metric, view="total") {
    medits_spatial_prepare(a,metric,view)
  }

  biomass_map_plot <- reactive({
    make_medits_spatial_plot(analysis(), "BIOMASS", get_view())
  })

  density_map_plot <- reactive({
    make_medits_spatial_plot(analysis(), "DENSITY", get_view())
  })

  # Dedicated compact maps for the PDF: smaller markers, fewer axis ticks.
  biomass_map_report_plot <- reactive({
    make_medits_spatial_plot(analysis(), "BIOMASS", get_view(), report=TRUE)
  })
  density_map_report_plot <- reactive({
    make_medits_spatial_plot(analysis(), "DENSITY", get_view(), report=TRUE)
  })

  output$lfd_allocation_status <- renderUI({
    a <- analysis()
    z <- medits_lfd_allocate(a,get_view())
    msg <- z$message
    if(nrow(z$audit)) {
      bad <- sum(z$audit$Negative_residual > 1e-5 * pmax(1,z$audit$Total),na.rm=TRUE)
      if(bad>0) msg <- paste0(msg,". WARNING: ",bad,
         " length-year cells have observed F+M exceeding total; inspect allocation QC CSV.")
    }
    tags$p(class="medits-hint",msg)
  })

  # Age readings, if provided, are for a descriptive ALK only. The assessment
  # algorithms and sex-specific raised LFDs do not consume the AGE table.
  alk_data <- reactive({
    req(input$age_file)
    a <- qc_data()
    req("AGE" %in% names(a))
    medits_alk_prepare(a$AGE,
       species=toupper(trimws(as.character(input$species))),
       gsa=if(identical(input$display_gsa,'AGGREGATED')) as.character(input$gsa)
           else sub('^GSA_','',input$display_gsa),
       country=as.character(input$country))
  })
  output$alk_data_status <- renderUI({
    if(is.null(input$age_file)) return(tags$p(class="medits-hint",
      "No age readings uploaded. Add the age-readings CSV to enable ALKs."))
    tryCatch({
      d <- alk_data()
      tags$p(class="medits-hint", paste(nrow(d),
        "valid aged observations available for the selected species, GSA and country. ",
        "The display is descriptive and not raised by survey area."))
    }, error=function(e) tags$p(class="medits-hint", conditionMessage(e)))
  })
  output$alk_year_ui <- renderUI({
    req(input$age_file)
    z <- tryCatch(alk_data(),error=function(e) NULL)
    if(is.null(z) || !nrow(z)) return(NULL)
    yy <- sort(unique(z$Year[is.finite(z$Year)]))
    selectInput("alk_year", "Year",choices=c("All years (pooled)"="ALL",
         setNames(as.character(yy),as.character(yy))),selected="ALL")
  })
  alk_table <- reactive({
    d <- alk_data()
    medits_alk_probabilities(d,sex=input$alk_sex,
        year=if(is.null(input$alk_year)) 'ALL' else input$alk_year)
  })
  alk_plot <- reactive({ medits_alk_plot(alk_table(),input$alk_sex) })
  output$alk_plot <- renderPlot({
    validate(need(!is.null(input$age_file),"Upload age readings CSV to generate the ALK."))
    alk_plot()
  })
  output$dl_alk_png <- downloadHandler(
    filename=function() paste0("MEDITS_ALK_",input$alk_sex,".png"),
    content=function(file) ggplot2::ggsave(file,alk_plot(),width=12,height=8,dpi=300))
  output$dl_alk_csv <- downloadHandler(
    filename=function() paste0("MEDITS_ALK_",input$alk_sex,".csv"),
    content=function(file) utils::write.csv(alk_table(),file,row.names=FALSE))

  ## ================= 1) HAUL MAP =====================================
  output$map_hauls <- renderPlot({
    haul_map_plot()
  })
  
  ## Export respects the all-years / one-year choice and displayed depth stratum.
  output$dl_haulmap_png <- downloadHandler(
    filename = function() if(isTRUE(input$haul_show_all_years)) "hauls_all_years.png" else paste0("haul_map_",input$map_year,".png"),
    content = function(file) {
      ggplot2::ggsave(file,plot=haul_map_plot(),
                      width=if(isTRUE(input$haul_show_all_years)) 14 else 12,
                      height=if(isTRUE(input$haul_show_all_years)) 12 else 9,
                      units="in",dpi=300)
    }
  )
  output$dl_haulmap_csv <- downloadHandler(
    filename = function() if(isTRUE(input$haul_show_all_years)) "haul_positions_all_years.csv" else paste0("haul_positions_",input$map_year,".csv"),
    content = function(file) {
      d <- medits_haul_map_data(analysis(),get_view(),
                                 all_years=isTRUE(input$haul_show_all_years),
                                 map_year=isolate(input$map_year))
      utils::write.csv(d$data,file,row.names=FALSE,na="")
    }
  )

  ## ================ 2) SURVEY PERIOD =================================
  # Un solo plot per interfaccia Shiny, PNG e report PDF.
  output$survey_period <- renderPlot({ survey_period_plot() })
  output$dl_surveyp_png <- downloadHandler(
    filename=function() "survey_period_hauls_by_year.png",
    content=function(file) {
      ggplot2::ggsave(file,plot=survey_period_plot(),width=12,height=9,dpi=300,units="in")
    }
  )
  output$dl_surveyp_csv <- downloadHandler(
    filename=function() "survey_period_hauls_by_year.csv",
    content=function(file) {
      write.csv(medits_survey_period_summary(analysis()$TAn),file,row.names=FALSE,na="")
    }
  )

  ## ================ 3) BIOMASS MAP (LOG) ============================
  output$map_biomass <- renderPlot({
    biomass_map_plot()
  })
  output$dl_biomap_png <- downloadHandler(
    filename = function() "biomass_map.png",
    content = function(file){
      ggsave(file, biomass_map_plot(), width = 14, height = 10, dpi = 300)
    }
  )
  output$dl_biomap_csv <- downloadHandler(
    filename = function() "biomass_map_data.csv",
    content = function(file) {
      write.csv(map_export_data(analysis(), "BIOMASS", get_view()), file, row.names=FALSE)
    }
  )

  ## ================ 4) DENSITY MAP (LOG) ============================
  output$map_density <- renderPlot({
    density_map_plot()
  })
  output$dl_denmap_png <- downloadHandler(
    filename = function() "density_map.png",
    content = function(file){
      ggsave(file, density_map_plot(), width = 14, height = 10, dpi = 300)
    }
  )
  output$dl_denmap_csv <- downloadHandler(
    filename = function() "density_map_data.csv",
    content = function(file) {
      write.csv(map_export_data(analysis(), "DENSITY", get_view()), file, row.names=FALSE)
    }
  )

  ## ================ 5) HAULS BY STRATUM =============================
  hauls_stratum_plot <- reactive({
    a   <- analysis()
    TAn <- a$TAn
    gsa_raw <- a$gsa
    gsa   <- gsa_label(gsa_raw)
    state <- a$state
    gen   <- a$gen
    spec  <- a$spec
    
    hauls_table_std <- aggregate(TAn$haul_number,
                                 list(TAn$year, TAn$strata),
                                 length)
    names(hauls_table_std) <- c("Year", "Stratum", "Hauls")
    
    hauls_table_std$Stratum <- as.character(hauls_table_std$Stratum)
    
    shelf_codes <- c("a","b","c")
    slope_codes <- c("d","e")
    
    hauls_shelf <- hauls_table_std %>%
      dplyr::filter(tolower(Stratum) %in% shelf_codes) %>%
      dplyr::group_by(Year) %>%
      dplyr::summarise(Hauls = sum(Hauls), .groups = "drop") %>%
      dplyr::mutate(Stratum = "Shelf (a+b+c)")
    
    hauls_slope <- hauls_table_std %>%
      dplyr::filter(tolower(Stratum) %in% slope_codes) %>%
      dplyr::group_by(Year) %>%
      dplyr::summarise(Hauls = sum(Hauls), .groups = "drop") %>%
      dplyr::mutate(Stratum = "Slope (d+e)")
    
    hauls_all <- dplyr::bind_rows(hauls_table_std, hauls_shelf, hauls_slope)
    
    ggplot(hauls_all,
           aes(x = Year, y = Hauls, colour = Stratum, group = Stratum)) +
      geom_line(size = 1.2) +
      geom_point(size = 3) +
      facet_wrap(~Stratum, scales = "free_y", ncol = 2) +
      theme_bw() +
      ggtitle(paste0(gen, spec,
                     " - Hauls by year & stratum (GSA", gsa, " ", state, ")")) +
      scale_x_continuous(breaks = seq(min(TAn$year), max(TAn$year), by = 2)) +
      theme(
        legend.position  = "none",
        strip.text       = element_text(size = 18, face = "bold"), 
        axis.title       = element_text(size = 18),
        axis.text        = element_text(size = 16),
        axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 16),
        plot.title       = element_text(size = 22, face = "bold", hjust = 0.5)
      )
  })
  
  output$hauls_by_stratum <- renderPlot({
    hauls_stratum_plot()
  })
  output$dl_haulstrat_png <- downloadHandler(
    filename = function() "hauls_by_stratum.png",
    content = function(file){
      ggsave(file, hauls_stratum_plot(), width = 14, height = 9, dpi = 300)
    }
  )
  output$dl_haulstrat_csv <- downloadHandler(
    filename = function() "hauls_by_stratum_data.csv",
    content = function(file){
      a   <- analysis()
      TAn <- a$TAn
      hauls_table_std <- aggregate(TAn$haul_number,
                                   list(TAn$year, TAn$strata),
                                   length)
      names(hauls_table_std) <- c("Year", "Stratum", "Hauls")
      write.csv(hauls_table_std, file, row.names = FALSE)
    }
  )
  
  ## ================ 6) BIOMASS / DENSITY / MEAN W ===================
  biom_trend_plot <- reactive({
    a <- analysis()
    B <- a$BIOMASS
    
    gsa_raw <- a$gsa
    gsa   <- gsa_label(gsa_raw)
    state <- a$state
    gen   <- a$gen
    spec  <- a$spec
    
    view <- get_view()
    
    # colonne originali
    ycol  <- "total_biomass"
    sdcol <- "stdev_tot"
    titolo_strato <- "Total (shelf + slope)"
    
    
    if (view == "shelf") {
      ycol  <- "shelf_biomass"
      sdcol <- "stdev_shelf"
      titolo_strato <- "Shelf (1–200 m)"
    } else if (view == "slope") {
      ycol  <- "slope_biomass"
      sdcol <- "stdev_slope"
      titolo_strato <- "Slope (201 m – max)"
    }
    
    Bplot <- B
    Bplot$plot_biomass <- Bplot[[ycol]]
    Bplot$plot_stdev   <- Bplot[[sdcol]]
    
    p <- ggplot(Bplot, aes(x = year, y = plot_biomass)) +
      geom_line(colour = "blue", size = 1.1) +
      geom_point(colour = "red", size = 3) +
      geom_errorbar(aes(ymin = plot_biomass - plot_stdev,
                        ymax = plot_biomass + plot_stdev),
                    width = 0.2,
                    linetype = "dashed",
                    colour = "red",
                    size = 0.8) +
      ylab("Biomass (kg/km²)") +
      xlab("Year") +
      ggtitle(paste0(gen, spec,
                     " - ", titolo_strato,
                     " biomass index (GSA", gsa, " ", state, ")")) +
      theme_bw() +
      theme(
        legend.position  = "bottom",
        strip.text       = element_text(size = 12, face = "bold"),
        axis.title       = element_text(size = 13),
        axis.text        = element_text(size = 11),
        plot.title       = element_text(size = 16, face = "bold", hjust = 0.5)
      )
    
    p
  })
  
  output$biom_trend <- renderPlot({ biom_trend_plot() })
  output$dl_biomtrend_png <- downloadHandler(
    filename = function() "biomass_trend.png",
    content = function(file){
      ggsave(file, biom_trend_plot(), width = 8, height = 4, dpi = 300)
    }
  )
  output$dl_biomtrend_csv <- downloadHandler(
    filename = function() "biomass_trend_data.csv",
    content = function(file){
      a <- analysis()
      write.csv(a$BIOMASS, file, row.names = FALSE)
    }
  )
  
  dens_trend_plot <- reactive({
    a <- analysis()
    D <- a$DENSITY
    
    gsa_raw <- a$gsa
    gsa   <- gsa_label(gsa_raw)
    state <- a$state
    gen   <- a$gen
    spec  <- a$spec
    
    view <- get_view()
    
    # colonne originali
    ycol  <- "total_density"
    sdcol <- "stdev_tot"
    titolo_strato <- "Total (shelf + slope)"
    
    
    if (view == "shelf") {
      ycol  <- "shelf_density"
      sdcol <- "stdev_shelf"
      titolo_strato <- "Shelf (1–200 m)"
    } else if (view == "slope") {
      ycol  <- "slope_density"
      sdcol <- "stdev_slope"
      titolo_strato <- "Slope (201 m – max)"
    }
    
    Dplot <- D
    Dplot$plot_density <- Dplot[[ycol]]
    Dplot$plot_stdev   <- Dplot[[sdcol]]
    
    p <- ggplot(Dplot, aes(x = year, y = plot_density)) +
      geom_line(colour = "blue", size = 1.1) +
      geom_point(colour = "red", size = 3) +
      geom_errorbar(aes(ymin = plot_density - plot_stdev,
                        ymax = plot_density + plot_stdev),
                    width = 0.2,
                    linetype = "dashed",
                    colour = "red",
                    size = 0.8) +
      ylab("Density (n/km²)") +
      xlab("Year") +
      ggtitle(paste0(gen, spec,
                     " - ", titolo_strato,
                     " density index (GSA", gsa, " ", state, ")")) +
      theme_bw() +
      theme(
        legend.position  = "bottom",
        strip.text       = element_text(size = 12, face = "bold"),
        axis.title       = element_text(size = 13),
        axis.text        = element_text(size = 11),
        plot.title       = element_text(size = 16, face = "bold", hjust = 0.5)
      )
    
    p
  })
  
  
  output$dens_trend <- renderPlot({ dens_trend_plot() })
  output$dl_denstrend_png <- downloadHandler(
    filename = function() "density_trend.png",
    content = function(file){
      ggsave(file, dens_trend_plot(), width = 14, height = 10, dpi = 300)
    }
  )
  output$dl_denstrend_csv <- downloadHandler(
    filename = function() "density_trend_data.csv",
    content = function(file){
      a <- analysis()
      write.csv(a$DENSITY, file, row.names = FALSE)
    }
  )
  
  mean_weight_plot <- reactive({
    a <- analysis()
    M <- medits_mean_weight_view(a,get_view())
    gsa_raw <- a$gsa
    gsa   <- gsa_label(gsa_raw)
    state <- a$state
    gen   <- a$gen
    spec  <- a$spec
    
    ggplot(M, aes(x = Year, y = `Mean_weight(g)`)) +
      geom_line(size = 1.1) +
      geom_point(size = 3, colour = "red") +
      theme_bw() +
      ggtitle(paste0(gen, spec,
                     " - Mean weight (GSA", gsa, " ", state, "; ", get_view(), ")")) +
      xlab("Year") + ylab("Mean weight (g)")+
      theme(
        legend.position  = "bottom",
        strip.text       = element_text(size = 12, face = "bold"),
        axis.title       = element_text(size = 13),
        axis.text        = element_text(size = 11),
        plot.title       = element_text(size = 16, face = "bold", hjust = 0.5))
  })
  
  output$mean_weight <- renderPlot({ mean_weight_plot() })
  output$dl_meanw_png <- downloadHandler(
    filename = function() "mean_weight.png",
    content = function(file){
      ggsave(file, mean_weight_plot(), width = 14, height = 10, dpi = 300)
    }
  )
  output$dl_meanw_csv <- downloadHandler(
    filename = function() "mean_weight_data.csv",
    content = function(file){
      a <- analysis()
      write.csv(medits_mean_weight_view(a,get_view()), file, row.names = FALSE)
    }
  )
  
  ## ================================================================
  ##   DA QUI IN POI: SOLO OGGETTI DI OUTPUT (ML, LFD, LFDF, LFDM, LFDSR, TCn)
  ## ================================================================
  
  ## ----------------------------------------------------------
  ## helpers per FLQuant da LFD
  ## ----------------------------------------------------------
  
  break_mode_ties <- function(m) {
    m2  <- as.matrix(m)
    eps <- 1e-9
    off <- eps * seq_len(nrow(m2))
    for (j in seq_len(ncol(m2))) {
      m2[, j] <- m2[, j] - off
    }
    m2
  }
  
  compute_length_stats_yearly <- medits_length_stats_yearly

  matrix_to_FLQuant <- function(mat, units = "mm") {
    mat <- as.matrix(mat)
    nlen <- nrow(mat)
    ny   <- ncol(mat)
    arr  <- array(mat, dim = c(nlen, ny, 1, 1, 1, 1))
    dimnames(arr) <- list(
      len    = rownames(mat),
      year   = colnames(mat),
      unit   = "all",
      season = "all",
      area   = "all",
      iter   = "1"
    )
    FLQuant(arr, units = units)
  }
  
  build_flq_total_from_LFD <- function(a, view) {
    # CRITICAL: a$LFD already contains all sex categories. Do NOT bind F+M again.
    LFD <- select_lfd(a,view,'C')
    if (!nrow(LFD)) stop('Empty total LFD for selected macro-stratum.')
    LFD1 <- stats::aggregate(Frequency ~ Year + Length,data=LFD,FUN=sum)
    names(LFD1) <- c('year','len','data')
    LFD1$unit <- 'C'
    LFD1$area <- paste0('GSA',gsa_label(a$gsa))
    LFD1$iter <- 1
    FLCore::as.FLQuant(LFD1,unit='mm')
  }

  build_flq_sex_from_LFD <- function(a, view, sex = c("F","M")) {
    sex <- match.arg(sex)
    LFD <- select_lfd(a,view,sex)
    unit_lab <- sex
    LFD$Len  <- LFD$Length
    LFD$Numb <- LFD$Frequency
    
    if (nrow(LFD) == 0) {
      stop(paste("Empty", unit_lab, "LFD for selected macro-stratum."))
    }
    
    LFD1 <- aggregate(Numb ~ Year * Len, data = LFD, FUN = sum)
    flip <- reshape2::dcast(LFD1, Year ~ Len, value.var = "Numb")
    
    lfd <- reshape2::melt(
      flip,
      id.vars       = "Year",
      variable.name = "Len",
      value.name    = "Numb"
    )
    
    names(lfd) <- c("year","len","data")
    lfd$unit  <- unit_lab
    lfd$area  <- paste0("GSA", gsa_label(a$gsa))
    lfd$iter  <- 1
    
    flq <- as.FLQuant(lfd, unit = "mm")
    flq[flq == 0] <- NA
    nmin <- apply(flq, c(2,5), function(x) {
      if (all(is.na(x))) return(NA_real_)
      min(x, na.rm = TRUE)
    })
    nmin[is.na(nmin)] <- 1
    flqs <- round(flq %/% nmin)
    flqs[is.na(flqs)] <- 0
    flq[is.na(flq)]   <- 0
    
    flq
  }
  
  ## ================ 7) MEAN LENGTH (usa ML dell’output) ==============
  
  mean_length_plot <- reactive({
    a    <- analysis()
    view <- get_view()
    d <- select_lfd(a,view,'C')
    d <- d[is.finite(d$Frequency) & d$Frequency >= 0,,drop=FALSE]
    ML <- medits_lfd_summary(d,'total')
    if (!nrow(ML)) stop('No valid mean length for selected area.')
    ML$Mean_Length_mm <- ML$Weighted_mean_length_mm
    view_label <- view
    
    if ("Mean_Length_mm" %in% names(ML)) {
      ML$Mean_Length <- ML$Mean_Length_mm
    } else if ("Mean_Length(mm)" %in% names(ML)) {
      ML$Mean_Length <- ML[["Mean_Length(mm)"]]
    } else if ("Mean_length" %in% names(ML)) {
      ML$Mean_Length <- ML$Mean_length
    } else {
      stop("Colonna di mean length non trovata in ML.")
    }
    
    gsa_raw <- a$gsa
    gsa   <- gsa_label(gsa_raw)
    state <- a$state
    gen   <- a$gen
    spec  <- a$spec
    
    ggplot(ML, aes(x = Year, y = Mean_Length)) +
      geom_line(size = 1.1) +
      geom_point(size = 3, colour = "red") +
      theme_bw() +
      ggtitle(paste0(gen, spec,
                     " - Mean length (", view_label,
                     ", GSA", gsa, " ", state, ")")) +
      xlab("Year") + ylab("Mean length (mm)")
  })
  
  output$mean_length <- renderPlot({
    mean_length_plot()
  })
  
  output$dl_meanlen_png <- downloadHandler(
    filename = function() "mean_length.png",
    content = function(file){
      ggsave(file, mean_length_plot(), width = 7, height = 5, dpi = 300)
    }
  )
  
  output$dl_meanlen_csv <- downloadHandler(
    filename = function() "mean_length_data.csv",
    content = function(file){
      a    <- analysis()
      ML <- medits_lfd_summary(select_lfd(a,get_view(),'C'),'total')
      write.csv(ML, file, row.names = FALSE)
    }
  )
  
  ## ================ 8) SAMPLE SIZE  TOTALE ==========================
  
  nsample_total_report_plot <- reactive({
    a   <- analysis()
    TCn <- a$TCn   # TC già filtrato per specie, GSA e country
    
    ns_tot <- aggregate(nblon ~ year, data = TCn, sum, na.rm = TRUE)
    
    ggplot(ns_tot, aes(x = year, y = nblon)) +
      geom_line(size = 1.1) +
      geom_point(size = 3, colour = "red") +
      theme_bw() +
      labs(
        title = "Total sample size per year (all sexes, raw TCn counts)",
        x = "Year",
        y = "Number of individuals"
      )
  })
  output$nsample_total <- renderPlot({ nsample_total_report_plot() })
  
  output$dl_nsample_tot_csv <- downloadHandler(
    filename = function() "nsample_total_data.csv",
    content = function(file){
      a   <- analysis()
      TCn <- a$TCn
      ns_tot <- aggregate(nblon ~ year, data = TCn, sum, na.rm = TRUE)
      write.csv(ns_tot, file, row.names = FALSE)
    }
  )
  
  output$dl_nsample_tot_png <- downloadHandler(
    filename = function() "nsample_total.png",
    content = function(file){
      a   <- analysis()
      TCn <- a$TCn
      ns_tot <- aggregate(nblon ~ year, data = TCn, sum, na.rm = TRUE)
      
      p <- ggplot(ns_tot, aes(x = year, y = nblon)) +
        geom_line(size = 1.1) +
        geom_point(size = 3, colour = "red") +
        theme_bw() +
        labs(
          title = "Total sample size per year (all sexes, raw TCn counts)",
          x = "Year",
          y = "Number of individuals"
        )
      
      ggsave(file, p, width = 10, height = 6, dpi = 300)
    }
  )
  
  ## ================ 8.1) LFD BLOCKS – TOTALE ========================
  
  lfd_blocks_total_report_plot <- reactive({
    medits_lfd_legacy_plot(analysis(),get_view(),'C')
  })
  output$lfd_blocks_total <- renderPlot({ lfd_blocks_total_report_plot() })

  output$dl_lfd_tot_csv <- downloadHandler(
    filename = function() paste0('LFD_total_',get_view(),'.csv'),
    content = function(file) write.csv(select_lfd(analysis(),get_view(),'C'),
                                      file,row.names=FALSE)
  )

  output$dl_lfd_tot_png <- downloadHandler(
    filename=function() paste0('LFD_total_',get_view(),'.png'),
    content=function(file) ggplot2::ggsave(file,lfd_blocks_total_report_plot(),
                                            width=14,height=16,dpi=300)
  )

  lfd_undet_reactive <- reactive({
    medits_lfd_legacy_plot(analysis(),get_view(),'U')
  })
  output$lfd_undet_plot <- renderPlot(lfd_undet_reactive())
  output$dl_lfd_u_csv <- downloadHandler(
    filename=function() paste0('LFD_undetermined_',get_view(),'.csv'),
    content=function(file) write.csv(select_lfd(analysis(),get_view(),'U'),file,row.names=FALSE)
  )
  output$dl_lfd_u_png <- downloadHandler(
    filename=function() paste0('LFD_undetermined_',get_view(),'.png'),
    content=function(file) ggplot2::ggsave(file,lfd_undet_reactive(),width=11,height=8,dpi=300)
  )

  ## ================ 8.2) FEMALE – sample size (numeri grezzi TCn) ====
  
  nsample_female_report_plot <- reactive({
    a   <- analysis()
    TCn <- a$TCn
    
    TCn_F <- subset(TCn, sex == "F")
    if (nrow(TCn_F) == 0) {
      validate(need(FALSE, "No females (sex == 'F') found in the data: plot not available."))
    }
    
    ns_F <- aggregate(
      nblon ~ year,
      data = TCn_F,
      sum, na.rm = TRUE
    )
    
    names(ns_F) <- c("Year", "N")
    
    gsa   <- gsa_label(a$gsa)
    state <- a$state
    gen   <- a$gen
    spec  <- a$spec
    
    ggplot(ns_F, aes(x = Year, y = N)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 3, colour = "red") +
      labs(
        title = paste0(gen, spec,
                       " - Female sample size (GSA", gsa, " ", state, ")"),
        x = "Year",
        y = "N. specimens"
      ) +
      theme_bw(base_size = 13) +
      theme(
        plot.title  = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
      )
  })
  output$nsample_female <- renderPlot({ nsample_female_report_plot() })
  
  output$dl_nsample_f_png <- downloadHandler(
    filename = function() "nsample_female.png",
    content = function(file){
      a   <- analysis()
      TCn <- a$TCn
      ns_F <- aggregate(
        nblon ~ year,
        data = subset(TCn, sex == "F"),
        sum, na.rm = TRUE
      )
      names(ns_F) <- c("Year", "N")
      
      p <- ggplot(ns_F, aes(x = Year, y = N)) +
        geom_line(linewidth = 1.2) +
        geom_point(size = 3, colour = "red") +
        labs(
          title = "Female Sample size",
          x = "Year", y = "N. specimens"
        ) +
        theme_bw(base_size = 13)
      
      ggsave(file, p, width = 12, height = 8, dpi = 300)
    }
  )
  
  output$dl_nsample_f_csv <- downloadHandler(
    filename = function() "nsample_female_data.csv",
    content = function(file){
      a   <- analysis()
      TCn <- a$TCn
      ns_F <- aggregate(
        nblon ~ year,
        data = subset(TCn, sex == "F"),
        sum, na.rm = TRUE
      )
      names(ns_F) <- c("Year", "N")
      write.csv(ns_F, file, row.names = FALSE)
    }
  )
  
  lfd_blocks_female_report_plot <- reactive({
    medits_lfd_legacy_plot(analysis(),get_view(),'F',
      corrected=identical(input$lfd_sex_mode,'corrected'))
  })
  output$lfd_blocks_female <- renderPlot({ lfd_blocks_female_report_plot() })

  output$dl_lfd_f_csv <- downloadHandler(
    filename=function() paste0("LFD_female_",input$lfd_sex_mode,"_",get_view(),".csv"),
    content=function(file) {
      z <- medits_lfd_allocate(analysis(),get_view())
      d <- if(identical(input$lfd_sex_mode,'corrected')) z$F_corrected else z$F_observed
      utils::write.csv(d,file,row.names=FALSE)
    })

  output$dl_lfd_f_png <- downloadHandler(
    filename=function() paste0('LFD_female_',input$lfd_sex_mode,'_',get_view(),'.png'),
    content=function(file) ggplot2::ggsave(file,lfd_blocks_female_report_plot(),
                                            width=14,height=16,dpi=300)
  )

  ## ================ 8.3) MALE – nsample + blocks ====================
  
  nsample_male_report_plot <- reactive({
    a   <- analysis()
    TCn <- a$TCn
    
    TCn_M <- subset(TCn, sex == "M")
    if (nrow(TCn_M) == 0) {
      validate(need(FALSE, "No males (sex == 'M') found in the data: plot not available."))
    }
    
    ns_M<- aggregate(
      nblon ~ year,
      data = TCn_M,
      sum, na.rm = TRUE
    )
    names(ns_M) <- c("Year", "N")
    
    gsa   <- gsa_label(a$gsa)
    state <- a$state
    gen   <- a$gen
    spec  <- a$spec
    
    ggplot(ns_M, aes(x = Year, y = N)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 3, colour = "red") +
      labs(
        title = paste0(gen, spec,
                       " - Male sample size (GSA", gsa, " ", state, ")"),
        x = "Year",
        y = "N. specimens"
      ) +
      theme_bw(base_size = 13) +
      theme(
        plot.title  = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
      )
  })
  output$nsample_male <- renderPlot({ nsample_male_report_plot() })
  
  output$dl_nsample_m_png <- downloadHandler(
    filename = function() "nsample_male.png",
    content = function(file){
      a   <- analysis()
      TCn <- a$TCn
      ns_M <- aggregate(
        nblon ~ year,
        data = subset(TCn, sex == "M"),
        sum, na.rm = TRUE
      )
      names(ns_M) <- c("Year", "N")
      
      p <- ggplot(ns_M, aes(x = Year, y = N)) +
        geom_line(linewidth = 1.2) +
        geom_point(size = 3, colour = "red") +
        labs(
          title = "Male sample size",
          x = "Year", y = "N. specimens"
        ) +
        theme_bw(base_size = 13)
      
      ggsave(file, p, width = 12, height = 8, dpi = 300)
    }
  )
  
  output$dl_nsample_m_csv <- downloadHandler(
    filename = function() "nsample_male_data.csv",
    content = function(file){
      a   <- analysis()
      TCn <- a$TCn
      ns_M <- aggregate(
        nblon ~ year,
        data = subset(TCn, sex == "M"),
        sum, na.rm = TRUE
      )
      names(ns_M) <- c("Year", "N")
      write.csv(ns_M, file, row.names = FALSE)
    }
  )
  
  lfd_blocks_male_report_plot <- reactive({
    medits_lfd_legacy_plot(analysis(),get_view(),'M',
      corrected=identical(input$lfd_sex_mode,'corrected'))
  })
  output$lfd_blocks_male <- renderPlot({ lfd_blocks_male_report_plot() })

  output$dl_lfd_m_csv <- downloadHandler(
    filename=function() paste0("LFD_male_",input$lfd_sex_mode,"_",get_view(),".csv"),
    content=function(file) {
      z <- medits_lfd_allocate(analysis(),get_view())
      d <- if(identical(input$lfd_sex_mode,'corrected')) z$M_corrected else z$M_observed
      utils::write.csv(d,file,row.names=FALSE)
    })

  output$dl_lfd_m_png <- downloadHandler(
    filename=function() paste0('LFD_male_',input$lfd_sex_mode,'_',get_view(),'.png'),
    content=function(file) ggplot2::ggsave(file,lfd_blocks_male_report_plot(),
                                            width=14,height=16,dpi=300)
  )

  ## ================ 9) SEX RATIO – usa LFDSR ========================
  
  sex_ratio_plot <- reactive({
    a    <- analysis()
    view <- get_view()
    
    SRdf <- select_sr(a,view)
    view_label <- if ("macro" %in% names(SRdf)) view else "all strata (source)"
    
    if ("macro" %in% names(SRdf)) {
      SRdf <- subset(SRdf, macro == view)
    }
    
    if (!("Year" %in% names(SRdf))) {
      stop("LFDSR deve contenere la colonna 'Year'.")
    }
    
    len_col <- if ("Length" %in% names(SRdf)) "Length" else "Len"
    
    # Length classes are the same as those used in LFD (mm), not cm.
    SRdf$Length_mm <- suppressWarnings(as.numeric(SRdf[[len_col]]))

    gsa   <- gsa_label(a$gsa)
    state <- a$state
    gen   <- a$gen
    spec  <- a$spec
    
    ggplot(SRdf, aes(x = Length_mm, y = SR, colour = factor(Year))) +
      geom_line(size = 0.8, na.rm = TRUE) +
      facet_wrap(~Year) +
      geom_hline(yintercept = 0.5, linetype = "dashed", colour = "black") +
      theme_bw() +
      ggtitle(paste0(gen, spec,
                     " - Sex ratio (F/(F+M)) (",
                     view_label, ", GSA", gsa, " ", state, ")")) +
      xlab("Length (mm)") + ylab("SR (F/(F+M))") +
      theme(legend.position = "none")
  })
  
  output$sex_ratio <- renderPlot({
    sex_ratio_plot()
  })
  
  output$dl_sr_png <- downloadHandler(
    filename = function() "sex_ratio.png",
    content = function(file){
      ggsave(file, sex_ratio_plot(), width = 10, height = 8, dpi = 300)
    }
  )
  
  output$dl_sr_csv <- downloadHandler(
    filename = function() "sex_ratio_data.csv",
    content = function(file){
      a    <- analysis()
      view <- get_view()
      SRdf <- select_sr(a,view)
      if ("macro" %in% names(SRdf)) {
        SRdf <- subset(SRdf, macro == view)
      }
      len_col <- if ("Length" %in% names(SRdf)) "Length" else "Len"
      SRdf$Length_mm <- suppressWarnings(as.numeric(SRdf[[len_col]]))
      write.csv(SRdf, file, row.names = FALSE)
    }
  )

  # Download every CSV and JPEG produced by the CURRENT run (not previous runs).
  # It includes all nine corrected/combined LFD, raw F/M/U, LOESS and QC,
  # together with the other MEDITS indices, tables and graphics.
  output$download_all_csv_jpeg <- downloadHandler(
    filename = function() {
      a <- analysis()
      yr <- paste(range(a$report_years, na.rm=TRUE), collapse="-")
      gs <- paste(a$report_gsas, collapse="-")
      paste0("MEDITS_ALL_CSV_JPEG_", a$report_species,
             "_GSA", gs, "_", yr, ".zip")
    },
    contentType = "application/zip",
    content = function(file) {
      a <- analysis()
      run_dir <- a$output_run_dir
      all_runs <- all_analyses()
      paths <- unlist(lapply(all_runs,function(r) {
        if (!dir.exists(r$output_run_dir)) return(character())
        list.files(r$output_run_dir, pattern="\\.(csv|jpe?g)$",
                   full.names=TRUE, recursive=TRUE, ignore.case=TRUE)
      }),use.names=FALSE)
      paths <- unique(paths[file.exists(paths) & !dir.exists(paths)])
      if(!length(paths)) stop('No CSV or JPEG files found')
      # Organize archive into CSV/ and JPEG/, retaining any subfolder names.
      stage <- tempfile("medits_zip_")
      dir.create(stage, recursive=TRUE)
      on.exit(unlink(stage,recursive=TRUE,force=TRUE),add=TRUE)
      for (p in paths) {
        origin <- names(all_runs)[vapply(all_runs,function(r)
          startsWith(normalizePath(p,winslash='/'),
                     paste0(normalizePath(r$output_run_dir,winslash='/'),'/')),logical(1))][1]
        if(is.na(origin)) next
        root <- all_runs[[origin]]$output_run_dir
        rel <- substring(normalizePath(p,winslash='/'),nchar(normalizePath(root,winslash='/'))+2L)
        ext <- tolower(tools::file_ext(p))
        bucket <- if(ext=='csv') 'CSV' else 'JPEG'
        dst <- file.path(stage,origin,bucket,rel)
        dir.create(dirname(dst),recursive=TRUE,showWarnings=FALSE)
        if (!file.copy(p,dst,overwrite=TRUE))
          stop("Cannot stage file for ZIP: ",p)
      }
      # Also export all available dashboard figures (current selected domain)
      # in JPEG format, without changing the other download buttons.
      dashboard_plots <- list(
        Haul_positions=haul_map_report_plot,
        Survey_period=survey_period_report_plot,
        Biomass_spatial_map=biomass_map_report_plot,
        Density_spatial_map=density_map_report_plot,
        Hauls_by_stratum=hauls_stratum_plot,
        Biomass_trend=biom_trend_plot,
        Density_trend=dens_trend_plot,
        Mean_weight=mean_weight_plot,
        Mean_length=mean_length_plot,
        Total_sample_size=nsample_total_report_plot,
        Female_sample_size=nsample_female_report_plot,
        Male_sample_size=nsample_male_report_plot,
        LFD_combined=lfd_blocks_total_report_plot,
        LFD_undetermined=lfd_undet_reactive,
        LFD_female_corrected=lfd_blocks_female_report_plot,
        LFD_male_corrected=lfd_blocks_male_report_plot,
        Sex_ratio=sex_ratio_plot
      )
      plot_issues <- character()
      view_name <- get_view()
      for (nm in names(dashboard_plots)) {
        tryCatch({
          plot_obj <- dashboard_plots[[nm]]()
          if (!inherits(plot_obj,"ggplot")) stop("Not a ggplot")
          outfile <- file.path(stage,"JPEG", "DASHBOARD",
                               paste0(nm,"_",view_name,".jpeg"))
          dir.create(dirname(outfile),recursive=TRUE,showWarnings=FALSE)
          ggplot2::ggsave(filename=outfile,plot=plot_obj,
                          width=12,height=9,dpi=300,units="in",limitsize=FALSE)
        }, error=function(e) {
          plot_issues <<- c(plot_issues,paste0(nm,": ",conditionMessage(e)))
        })
      }
      # Export the exact figure series used by the PDF report, at print quality.
      # These reactives are the same ones passed to medits_report_write_pdf().
      # The report is for the user-selected bathymetric domain.
      report_plots <- list(
        "Haul positions" = haul_map_report_plot,
        "Survey period" = survey_period_report_plot,
        "Biomass spatial map" = biomass_map_report_plot,
        "Density spatial map" = density_map_report_plot,
        "Hauls by stratum" = hauls_stratum_plot,
        "Biomass trend" = biom_trend_plot,
        "Density trend" = dens_trend_plot,
        "Mean weight" = mean_weight_plot,
        "Mean length" = mean_length_plot,
        "Total sample size" = nsample_total_report_plot,
        "Female sample size" = nsample_female_report_plot,
        "Male sample size" = nsample_male_report_plot,
        "LFD total" = lfd_blocks_total_report_plot,
        "LFD undetermined" = lfd_undet_reactive,
        "LFD females" = lfd_blocks_female_report_plot,
        "LFD males" = lfd_blocks_male_report_plot,
        "Sex ratio" = sex_ratio_plot
      )
      report_issues <- character()
      for (i in seq_along(report_plots)) {
        nm <- names(report_plots)[i]
        tryCatch({
          plot_obj <- report_plots[[i]]()
          if (!inherits(plot_obj, "ggplot")) stop("Not a ggplot")
          safe_name <- gsub("[^A-Za-z0-9]+", "_", nm)
          outfile <- file.path(stage, "JPEG", "REPORT_FIGURES",
                               sprintf("Figure_%02d_%s_%s.jpeg", i, safe_name, view_name))
          dir.create(dirname(outfile), recursive=TRUE, showWarnings=FALSE)
          ggplot2::ggsave(filename=outfile, plot=plot_obj,
                          width=12, height=9, dpi=300, units="in",
                          limitsize=FALSE)
        }, error=function(e) {
          report_issues <<- c(report_issues, paste0("Figure ",i," (",nm,"): ",conditionMessage(e)))
        })
      }
      write.csv(data.frame(Figure=seq_along(report_plots),
                           Title=names(report_plots),
                           Bathymetric_view=rep(view_name,length(report_plots))),
                file.path(stage,"REPORT_FIGURES_INDEX.csv"),row.names=FALSE)
      write.csv(corrections(),file.path(stage,'QC_MANUAL_CORRECTIONS.csv'),row.names=FALSE)
      if(isTRUE(input$qc_scan>0)) {
        try(utils::write.csv(qc_all_results()$candidates,
          file.path(stage,'QC_FLAGGED_RECORDS.csv'),row.names=FALSE),silent=TRUE)
        try(utils::write.csv(qc_all_results()$fields,
          file.path(stage,'QC_FIELD_SUMMARY.csv'),row.names=FALSE),silent=TRUE)
      }
      # GSA-specific LFDs and JPEGs, plus one combined long-format CSV.
      medits_export_lfd_by_gsa(all_runs,stage,plots=TRUE)
      medits_export_metrics_by_gsa(all_runs,stage)
      included <- sort(list.files(stage,recursive=TRUE,full.names=FALSE))
      writeLines(c(
        "MEDITS - files from ONE analysis run",
        paste("Species:",a$report_species),
        paste("GSA:",paste(a$report_gsas,collapse=",")),
        paste("Years:",paste(range(a$report_years),collapse="-")),
        paste("Source run:",basename(run_dir)),
        "Sex-ratio allocation: one global LOESS (10-800 m, all observed years)",
        "Original model LFDs are preserved; separate LOESS-corrected F/M CSVs and allocation QC are included.",
        paste0("PDF report JPEG figure set: JPEG/REPORT_FIGURES (selected domain: ",view_name,")."),
        "Observed F/M/U and allocation QC files are also included.",
        "See CSV/LFD_BY_GSA: one GSA per file and LFD_ALL_SELECTED_GSAS.csv.",
        "Shelf/slope LFD are exported only if actual standardized objects exist.",
        paste("Files:",length(included)),
        if(length(plot_issues)) c("Dashboard plot warnings:",plot_issues) else "All dashboard plots generated successfully.",
        if(length(report_issues)) c("Report figure warnings:",report_issues) else "All PDF report figures generated successfully.",
        "",included), file.path(stage,"MANIFEST.txt"))
      # zip::zipr works on Windows and does not require a system zip binary.
      if (!requireNamespace("zip",quietly=TRUE))
        stop("Install the zip package once: install.packages('zip')")
      zip::zipr(zipfile=file,
                files=c("MANIFEST.txt",included),root=stage,
                include_directories=FALSE)
      if (!file.exists(file) || file.info(file)$size<=0)
        stop("ZIP creation failed")
    }
  )

  # Explicit LFD data for EACH selected GSA, never inferred by splitting pooled LFD.
  # Area-stratified estimates cannot be obtained by dividing aggregate estimates.
  medits_export_lfd_by_gsa <- function(all_runs, stage, plots=TRUE) {
    output_path <- file.path(stage,"CSV","LFD_BY_GSA")
    dir.create(output_path,recursive=TRUE,showWarnings=FALSE)
    sex_objects <- c(TOTAL="LFD",FEMALE="LFDF",MALE="LFDM",UNDETERMINED="LFDU")
    all_data <- list()
    catalog <- list()
    selected <- names(all_runs)[grepl("^GSA_",names(all_runs))]
    for (nm in selected) {
      a <- all_runs[[nm]]
      gsa <- sub("^GSA_","",nm)
      for (view in c("total","shelf","slope")) {
        # Separate, explicitly named derived output; leave original model LFDs untouched.
        alloc <- medits_lfd_allocate(a,view)
        if(nrow(alloc$audit)) {
          utils::write.csv(alloc$audit,
            file.path(output_path,paste0("LFD_GSA_",gsa,"_",view,"_LOESS_ALLOCATION_QC.csv")),
            row.names=FALSE)
          if(isTRUE(alloc$allocated)) {
            utils::write.csv(alloc$F_corrected,
              file.path(output_path,paste0("LFD_GSA_",gsa,"_",view,"_FEMALE_CORRECTED.csv")),row.names=FALSE)
            utils::write.csv(alloc$M_corrected,
              file.path(output_path,paste0("LFD_GSA_",gsa,"_",view,"_MALE_CORRECTED.csv")),row.names=FALSE)
          }
        }
        for (sex in names(sex_objects)) {
          key <- sex_objects[[sex]]
          if (view!="total") key <- paste0(key,"_",view)
          df <- a[[key]]
          valid <- is.data.frame(df) && all(c("Year","Length","Frequency") %in% names(df)) && nrow(df)>0
          catalog[[length(catalog)+1L]] <- data.frame(GSA=gsa,Domain=view,Sex=sex,
            SourceObject=key,Status=if(valid) "EXPORTED" else "NO_DATA_OR_NOT_IMPLEMENTED")
          if(!valid) next
          df <- as.data.frame(df)
          df$GSA <- as.character(gsa)
          df$Domain <- view
          df$Sex <- sex
          first_cols <- c("GSA","Domain","Sex","Year","Length","Frequency")
          df <- df[,c(first_cols,setdiff(names(df),first_cols)),drop=FALSE]
          outfile <- file.path(output_path,paste0("LFD_GSA_",gsa,"_",view,"_",sex,".csv"))
          utils::write.csv(df,outfile,row.names=FALSE)
          all_data[[length(all_data)+1L]] <- df[,first_cols,drop=FALSE]
          if (plots) {
            tryCatch({
              years <- length(unique(df$Year))
              fig <- ggplot2::ggplot(df,ggplot2::aes(x=Length,y=Frequency))+
                ggplot2::geom_col(fill="#2A7182",width=0.9)+
                ggplot2::facet_wrap(~Year,scales="free_y",ncol=4)+
                ggplot2::labs(title=paste("MEDITS LFD – GSA",gsa,"–",sex,"–",view),
                              x="Length (source units)",y="Frequency (standardized)")+
                ggplot2::theme_minimal(base_size=12)+
                ggplot2::theme(strip.text=ggplot2::element_text(size=10,face="bold"))
              dest <- file.path(stage,"JPEG","LFD_BY_GSA",
                            paste0("LFD_GSA_",gsa,"_",view,"_",sex,".jpeg"))
              dir.create(dirname(dest),recursive=TRUE,showWarnings=FALSE)
              ggplot2::ggsave(dest,plot=fig,width=15,
                              height=min(34,max(8,ceiling(years/4)*2.0)),dpi=220,units="in",limitsize=FALSE)
            },error=function(e) warning("GSA plot ",nm," ",sex,": ",conditionMessage(e)))
          }
        }
      }
    }
    if(length(all_data)) {
      common <- do.call(rbind,all_data)
      common <- common[order(common$GSA,common$Domain,common$Sex,common$Year,common$Length),,drop=FALSE]
      utils::write.csv(common,file.path(output_path,"LFD_ALL_SELECTED_GSAS.csv"),row.names=FALSE)
    }
    status <- attr(all_runs,"gsa_status")
    if(is.data.frame(status)) utils::write.csv(status,file.path(output_path,"GSA_RUN_STATUS.csv"),row.names=FALSE)
    if(length(catalog)) utils::write.csv(do.call(rbind,catalog),
                                         file.path(output_path,"LFD_AVAILABILITY_BY_GSA.csv"),row.names=FALSE)
    invisible(length(all_data))
  }

  # Save original independently calculated model output tables for each GSA
  # (and the genuine aggregate run) without reconstructing from pooled results.
  medits_export_metrics_by_gsa <- function(all_runs,stage) {
    base <- file.path(stage,'CSV','INDICES_BY_GSA')
    dir.create(base,recursive=TRUE,showWarnings=FALSE)
    table_keys <- c('BIOMASS','DENSITY',
                    'Mweight','Mweight_domain','ML','LFDSR','TAn')
    for(nm in names(all_runs)) {
      a <- all_runs[[nm]]
      for(key in table_keys) {
        d <- a[[key]]
        if(!is.data.frame(d)||!nrow(d)) next
        d <- as.data.frame(d); d$GSA_Analysis <- nm
        name <- paste0('MEDITS_',nm,'_',key,'.csv')
        utils::write.csv(d,file.path(base,name),row.names=FALSE)
      }
    }
    for(domain in c('total','shelf','slope')) {
      d <- medits_collect_gsa_trends(all_runs,domain)
      if(nrow(d)) utils::write.csv(d,file.path(base,paste0('MEDITS_TRENDS_',domain,'_ALL_GSAS.csv')),row.names=FALSE)
    }
    st <- attr(all_runs,'gsa_status')
    if(is.data.frame(st)) utils::write.csv(st,file.path(base,'GSA_RUN_STATUS.csv'),row.names=FALSE)
    invisible(TRUE)
  }

  output$download_lfd_gsa <- downloadHandler(
    filename=function() paste0("MEDITS_LFD_BY_GSA_",analysis()$report_species,".zip"),
    contentType="application/zip",
    content=function(file) {
      all_runs <- all_analyses()
      stage <- tempfile("medits_gsa_lfd_")
      dir.create(stage,recursive=TRUE)
      on.exit(unlink(stage,recursive=TRUE,force=TRUE),add=TRUE)
      medits_export_lfd_by_gsa(all_runs,stage,plots=TRUE)
      medits_export_metrics_by_gsa(all_runs,stage)
      files <- list.files(stage,recursive=TRUE,full.names=FALSE)
      if(!length(files)) stop("No GSA LFD export available. Check GSA_RUN_STATUS.csv")
      if(!requireNamespace("zip",quietly=TRUE)) stop("Install package 'zip'")
      zip::zipr(zipfile=file,files=files,root=stage,include_directories=FALSE)
    })

  ## PDF includes cover, figure/table index, all figures and scientific captions.
  ## Automatic interpretations use the numeric analysis results, NOT image recognition.
  ## Detailed length-frequency and sex-ratio tables can be disabled for a shorter report.
  output$download_full_report <- downloadHandler(
    filename = function() {
      a <- analysis()
      yr <- paste(range(a$report_years, na.rm = TRUE), collapse = "-")
      gs <- paste(a$report_gsas, collapse = "-")
      paste0("MEDITS_Report_", a$report_species, "_GSA", gs, "_", yr, ".pdf")
    },
    contentType = "application/pdf",
    content = function(file) {
      a <- analysis()
      # Use exactly the same plot reactives as the application.
      plots <- list(
        "Haul positions" = haul_map_report_plot,
        "Survey period" = survey_period_report_plot,
        "Biomass spatial map" = biomass_map_report_plot,
        "Density spatial map" = density_map_report_plot,
        "Hauls by stratum" = hauls_stratum_plot,
        "Biomass trend" = biom_trend_plot,
        "Density trend" = dens_trend_plot,
        "Mean weight" = mean_weight_plot,
        "Mean length" = mean_length_plot,
        "Total sample size" = nsample_total_report_plot,
        "Female sample size" = nsample_female_report_plot,
        "Male sample size" = nsample_male_report_plot,
        "LFD total" = lfd_blocks_total_report_plot,
        "LFD undetermined" = lfd_undet_reactive,
        "LFD females" = lfd_blocks_female_report_plot,
        "LFD males" = lfd_blocks_male_report_plot,
        "Sex ratio" = sex_ratio_plot
      )
      medits_report_write_pdf(
        file = file,
        a = a,
        view = get_view(),
        plot_reactives = plots,
        map_year = isolate(input$map_year),
        include_detail = isTRUE(input$report_include_detailed_tables),
        include_interpretation = isTRUE(input$report_include_interpretation)
      )
    }
  )

}

shinyApp(ui = ui, server = server)