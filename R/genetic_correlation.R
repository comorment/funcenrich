

get_cor <- function(m) {

  res <- m[grepl("T1~~T2", paste0(m$lhs, m$op, m$rhs)), c("STD_Genotype" , "STD_Genotype_SE")]
  tibble::tibble(rg = res[[1]], rg_se = as.numeric(res[[2]]))

}



#' Genetic correlation
#'
#' @param trait1 Name of trait 1
#' @param trait2 Name of trait 2
#' @param gen_cov Genetic covariance matrix estimated by ldsc
#' @param adjust Name of trait to adjust for
#'
#' @return Tibble with genetic correlation between trait1 and trait2 adjusted for genetic components of trait in adjust
#' @export
#'
#' @examples NA
est_genetic_correlation <- function(trait1, trait2, gen_cov, adjust = NULL) {


  if (is.null(adjust)) {

    res_traits <- tibble::tibble(Trait1 = trait1, Trait2 = trait2, Adjust_for = "Nothing")
    spec <- glue::glue("T1 =~ {trait1}", "T2 =~ {trait2}", .sep = "\n")

  } else {

    res_traits <- tibble::tibble(Trait1 = trait1, Trait2 = trait2, Adjust_for = adjust)
    spec <- glue::glue("T1 =~ {trait1}", "T2 =~ {trait2}", "T1 ~ {adjust}", "T2 ~ {adjust}", .sep = "\n")

  }

  m <- GenomicSEM::usermodel(gen_cov,
                             estimation = "DWLS",
                             model = spec,
                             CFIcalc = TRUE,
                             std.lv = TRUE,
                             imp_cov = FALSE)

  res <- get_cor(m$results)
  dplyr::bind_cols(res_traits, res)


}
