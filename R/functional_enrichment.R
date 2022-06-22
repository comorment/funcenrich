

#' Enrichment in genomic regions of a variance of a common latent factor
#'
#' @param traits Traits making up the common latent factor, should beat least 3
#' @param gen_cov_strat Stratified genetic covariance matrices
#'
#' @return Tibble of enrichment in regulatory genomic regions
#' @export
#'
#' @examples NA
est_enrichment <- function(traits, gen_cov_strat) {

  if( length(traits) < 3 ) stop("Need at least three traits to model common factor")

  traits_sum <- stringr::str_c(traits, collapse = " + ")
  traits_sum <- glue::glue("F =~ {traits_sum}")
  traits_var <- purrr::map2(traits, letters[1:length(traits)], ~ glue::glue("{.x} ~~ {.y} * {.x}"))
  traits_const <- glue::glue("{letters[1:length(traits)]} > 1 / 1e4")
  spec_list <- c(traits_sum, traits_var, traits_const)
  spec <- rlang::inject(glue::glue(!!!spec_list, .sep = "\n"))
  params <- "F ~~ F"
  res <- GenomicSEM::enrich(s_covstruc = gen_cov_strat, model = spec, params = params, std.lv = TRUE)[[1]]
  dplyr::arrange(tibble::as_tibble(res), desc(Enrichment - 1.96 * Enrichment_SE))

}
