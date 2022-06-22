
download_and_untar <- function(filepath, outdir = FALSE) {

  file <- paste0("REF/", basename(filepath))
  if (!fs::dir_exists("REF")) fs::dir_create("REF")
  download.file(filepath, destfile = file)
  if (!outdir) untar(file, exdir = "REF")
  else untar(file, exdir = paste0("REF/", stringr::str_extract(basename(filepath), ".*(?=.tgz)")))
  fs::file_delete(file)
  NULL

}


#' Download reference files
#'
#' Downloads needed reference files from \url{https://alkesgroup.broadinstitute.org/LDSCORE/} and \url{https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/eur_w_ld_chr.tar.bz2}
#'
#' @return NULL
#' @export
#'
#' @examples NA
download_reference_data <- function() {

  download_and_untar("https://alkesgroup.broadinstitute.org/LDSCORE/1000G_Phase3_weights_hm3_no_MHC.tgz")
  download_and_untar("https://alkesgroup.broadinstitute.org/LDSCORE/1000G_Phase3_frq.tgz")
  download_and_untar("https://alkesgroup.broadinstitute.org/LDSCORE/1000G_Phase3_baselineLD_v2.2_ldscores.tgz", outdir = TRUE)
  download_and_untar("https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/eur_w_ld_chr.tar.bz2")

}
