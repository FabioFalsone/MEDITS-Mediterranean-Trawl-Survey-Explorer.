# Run in RStudio in the same folder as the app. Static syntax checks and
# source contract assertions; does not run the full assessment.
for (file in c('Medits_app_GSA.R','Medits_function_GSAMedits_gsa.R',
               'medits_report_helpers.R')) {
  invisible(parse(file=file))
  cat('Parsed:', file, '\n')
}
app <- paste(readLines('Medits_app_GSA.R',warn=FALSE),collapse='\n')
engine <- paste(readLines('Medits_function_GSAMedits_gsa.R',warn=FALSE),collapse='\n')
stopifnot(grepl('Mweight_domain = Mweight_domain',engine,fixed=TRUE))
stopifnot(grepl('medits_mean_weight_view(a,get_view())',app,fixed=TRUE))
stopifnot(grepl('medits_mean_weight_view(a,domain)',app,fixed=TRUE))
stopifnot(grepl("'Mweight_domain'",app,fixed=TRUE))
cat('Mean-weight domain source checks: OK\n')