

#' Genetic covariance between traits using LDScore regression
#'
#' @param paths Paths to munged summary statistics
#' @param ld_score_annotation LDScores for annotation
#' @param sample_prev Sample prevalences for GWASs. NA for continuous
#' @param population_prev Population prevalences for traits. NA for continuous
#' @param trait_names Specify names for traits
#'
#' @return Stratified genetic covariance matrix estimated by GenomicSEM::s_ldsc
#' @export
#'
#' @examples NA
est_stratified_gen_cov <- function(paths, ld_score_annotation = "REF/1000G_Phase3_baselineLD_v2.2_ldscores/", sample_prev, population_prev, trait_names) {

  wld <- "REF/1000G_Phase3_weights_hm3_no_MHC/"
  frq <- "REF/1000G_Phase3_frq/"

  if(!fs::dir_exists(ld_score_annotation))  stop("LDscores not present, please run download_reference_data()")
  if(!fs::dir_exists(wld)) stop("Weights not present, please run download_reference_data()")
  if(!fs::dir_exists(frq)) stop("Frequencies not present, please run download_reference_data()")


  GenomicSEM::s_ldsc(paths,
         sample.prev = sample_prev,
         population.prev = population_prev,
         ld = ld_score_annotation,
         wld = wld,
         frq = frq,
         trait.names = trait_names)


}


#' Stratified genetic covariance between traits using LDScore regression
#'
#' @param paths Paths to munged summary statistics
#' @param sample_prev Sample prevalences for GWASs. NA for continuous
#' @param population_prev Population prevalences for traits. NA for continuous
#' @param trait_names Specify names for traits
#'
#' @return Genetic covariance object estimated from GenomicSEM::ldsc
#' @export
#'
#' @examples NA
est_gen_cov <- function(paths, sample_prev, population_prev, trait_names) {


  ld <- "REF/eur_w_ld_chr/"
  if(!fs::dir_exists(ld)) stop("LD not present, please run download_reference_data()")
  wld <- ld
  GenomicSEM::ldsc(paths,
       sample.prev = sample_prev,
       population.prev = population_prev,
       ld = ld,
       wld = wld,
       trait.names = trait_names)


}
