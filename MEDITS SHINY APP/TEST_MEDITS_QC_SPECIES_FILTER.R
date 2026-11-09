# Run from the MEDITS application folder after installing data.table.
source('medits_qc_helpers.R')
# Direct filtering test uses a tiny local lookup with same structure as MEDITS.
wd <- tempfile('medits_species_test_');dir.create(wd)
on.exit(unlink(wd,recursive=TRUE),add=TRUE)
sp <- data.frame(A=1:2, Scientific_name=c('Species alpha','Species beta'),
                 C=1:2,D=1:2,E=1:2,F=1:2,code=c('0001001','0002002'))
as <- data.frame(A=1:2,B=1:2,Scientific_name=sp$Scientific_name,
                 `3A_CODE`=c('HKE','DPS'),check.names=FALSE)
utils::write.csv(sp,file.path(wd,'Sp_Medits.csv'),row.names=FALSE)
utils::write.csv(as,file.path(wd,'ASFIS_2020.csv'),row.names=FALSE)
x <- data.frame(GENUS=c('0001','0002','0001'),SPECIES=c('001','002','001'),
                AREA=c(16,16,15),COUNTRY=c('ITA','ITA','MLT'), YEAR=2025,
                LENGTH_CLASS=c(20,200,10), WEIGHT_G=c(5,5000,2),stringsAsFactors=FALSE)
rows <- medits_qc_species_rows(list(TC=x), 'HKE',gsa=16,country='ITA',
  species_file=file.path(wd,'Sp_Medits.csv'),asfis_file=file.path(wd,'ASFIS_2020.csv'))
stopifnot(identical(rows$TC,1L))
scan <- medits_qc_scan_data(list(TC=x),max_length=10,max_weight=1,max_age=20,
                            selected_rows=rows)
stopifnot(nrow(scan$candidates)>0,all(scan$candidates$Row==1L),
          !any(scan$candidates$Row %in% c(2L,3L)))
message('PASS: QC returns selected species/GSA/country only, preserving source row indices.')