
<!-- README.md is generated from README.Rmd. Please edit that file -->

# funcenrich

<!-- badges: start -->
<!-- badges: end -->

funcenrich implements some genetic correlation and functional enrichment
analysis performed for work package 4 of the CoMorMent H2020 EU grant.

## Installation

The package can be installed by

``` r
devtools::install_github("comorment/funcenrich")
```

## Example

The package comes with three prepared munged summary stats from publicly
available summary statistics of GWASs of coronary artery disease (CAD),
anxiety (ANX) and BMI, which have been munged to Hapmap3 SNPs. We have
done some internal processing of the data so it wont be exactly like
what you would get if you download it yourself. The summary statistics
are supplied as tsv files and can be access from the paths.

``` r
library(funcenrich)
path_anx <- system.file("extdata", "anx", package = "funcenrich")
path_bmi <- system.file("extdata", "bmi", package = "funcenrich")
path_cad <- system.file("extdata", "cad", package = "funcenrich")
paths <- c(path_anx, path_bmi, path_cad)
```

The CAD data comes from

Nelson, C., Goel, A., Butterworth, A. et al. Association analyses based
on false discovery rate implicate new loci for coronary artery disease.
Nat Genet 49, 1385–1391 (2017). <https://doi.org/10.1038/ng.3913>

and can be downloaded from
<http://www.cardiogramplusc4d.org/media/cardiogramplusc4d-consortium/data-downloads/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt.gz>.
The ANX data comes from

Purves, K.L., Coleman, J.R.I., Meier, S.M. et al. A major role for
common genetic variation in anxiety disorders. Mol Psychiatry 25,
3292–3303 (2020). <https://doi.org/10.1038/s41380-019-0559-1> and can be
downloaded from <https://www.kcl.ac.uk/people/kirstin-purves>. Finally,
the BMI data comes from

Loic Yengo, Julia Sidorenko, Kathryn E Kemper, Zhili Zheng, Andrew R
Wood, Michael N Weedon, Timothy M Frayling, Joel Hirschhorn, Jian Yang,
Peter M Visscher, the GIANT Consortium, Meta-analysis of genome-wide
association studies for height and body mass index in ∼700000
individuals of European ancestry, Human Molecular Genetics, Volume 27,
Issue 20, 15 October 2018, Pages 3641–3649,
<https://doi.org/10.1093/hmg/ddy271>

and can be downloaded from
<http://cnsgenomics.com/data/yengo_et_al_2018_hmg/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt.gz>

### Reference files

First we need a bunch of reference files. The function
`download_reference_data()` will download the necessary files from
<https://alkesgroup.broadinstitute.org/LDSCORE> and
<https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/> and
will put them in a `REF` directory in the working directory.

``` r
download_reference_data()
#> NULL
```

Let’s define trait names and include sample- and population prevalences
taken from the original publications and the global burden of disease
project.

``` r
pop_prev <- c(0.0582, NA, 0.215)
sample_prev <- c(0.3, NA, 0.11)
trait_names <- c("ANX", "BMI", "CAD")
```

We can now compute the genetic covariance matrix between the traits and
compute the genetic correlation between CAD and ANX as well as the
genetic correlation between those two traits adjusted for the genetic
component of BMI.

``` r
  gen_cov <- est_gen_cov(paths, sample_prev = sample_prev, population_prev = pop_prev, trait_names = trait_names)
#> Warning: replacing previous import 'gdata::nobs' by 'lavaan::nobs' when loading
#> 'GenomicSEM'
#> Warning: replacing previous import 'gdata::last' by 'data.table::last' when
#> loading 'GenomicSEM'
#> Warning: replacing previous import 'gdata::first' by 'data.table::first' when
#> loading 'GenomicSEM'
#> Warning: replacing previous import 'gdata::env' by 'R.utils::env' when loading
#> 'GenomicSEM'
#> Warning: replacing previous import 'gdata::resample' by 'R.utils::resample' when
#> loading 'GenomicSEM'
#> Warning: replacing previous import 'data.table::last' by 'dplyr::last' when
#> loading 'GenomicSEM'
#> Warning: replacing previous import 'plyr::summarize' by 'dplyr::summarize' when
#> loading 'GenomicSEM'
#> Warning: replacing previous import 'plyr::mutate' by 'dplyr::mutate' when
#> loading 'GenomicSEM'
#> Warning: replacing previous import 'plyr::id' by 'dplyr::id' when loading
#> 'GenomicSEM'
#> Warning: replacing previous import 'plyr::arrange' by 'dplyr::arrange' when
#> loading 'GenomicSEM'
#> Warning: replacing previous import 'plyr::summarise' by 'dplyr::summarise' when
#> loading 'GenomicSEM'
#> Warning: replacing previous import 'data.table::first' by 'dplyr::first' when
#> loading 'GenomicSEM'
#> Warning: replacing previous import 'plyr::rename' by 'dplyr::rename' when
#> loading 'GenomicSEM'
#> Warning: replacing previous import 'plyr::desc' by 'dplyr::desc' when loading
#> 'GenomicSEM'
#> Warning: replacing previous import 'plyr::count' by 'dplyr::count' when loading
#> 'GenomicSEM'
#> Warning: replacing previous import 'gdata::combine' by 'dplyr::combine' when
#> loading 'GenomicSEM'
#> Warning: replacing previous import 'data.table::between' by 'dplyr::between'
#> when loading 'GenomicSEM'
#> Warning: replacing previous import 'plyr::failwith' by 'dplyr::failwith' when
#> loading 'GenomicSEM'
#> Warning: replacing previous import 'Matrix::tail' by 'utils::tail' when loading
#> 'GenomicSEM'
#> Warning: replacing previous import 'Rcpp::.DollarNames' by 'utils::.DollarNames'
#> when loading 'GenomicSEM'
#> Warning: replacing previous import 'Rcpp::prompt' by 'utils::prompt' when
#> loading 'GenomicSEM'
#> Warning: replacing previous import 'Matrix::head' by 'utils::head' when loading
#> 'GenomicSEM'
#> Warning: replacing previous import 'R.utils::timestamp' by 'utils::timestamp'
#> when loading 'GenomicSEM'
#> Warning: replacing previous import 'gdata::object.size' by 'utils::object.size'
#> when loading 'GenomicSEM'
#> Multivariate ld-score regression of 3 traits (C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/anx C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/bmi C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/cad) began at: 2022-06-28 10:13:13
#> Reading in LD scores
#> Read in summary statistics [1/3] from: C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/anx
#> Out of 1100641 SNPs, 1097932 remain after merging with LD-score files
#> Removing 0 SNPs with Chi^2 > 83.566; 1097932 remain
#> Read in summary statistics [2/3] from: C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/bmi
#> Out of 1019865 SNPs, 1014995 remain after merging with LD-score files
#> Removing 6 SNPs with Chi^2 > 795.64; 1014989 remain
#> Read in summary statistics [3/3] from: C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/cad
#> Out of 1186097 SNPs, 1173849 remain after merging with LD-score files
#> Removing 23 SNPs with Chi^2 > 154.654; 1173826 remain
#> Estimating heritability [1/6] for: C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/anx
#> Please note that the results initially printed to the screen and log file reflect the NON-liability h2 and cov_g. However, a liability conversion is being used for trait C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/anx when creating the genetic covariance matrix used as input for Genomic SEM and liability scale results are printed at the end of the log file.
#> Heritability Results for trait: C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/anx
#> Mean Chi^2 across remaining SNPs: 1.1835
#> Lambda GC: 1.1714
#> Intercept: 1.008 (0.0089)
#> Ratio: 0.0436 (0.0484)
#> Total Observed Scale h2: 0.1048 (0.0076)
#> h2 Z: 13.9
#> Calculating genetic covariance [2/6] for traits: C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/anx and C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/bmi
#> 954338 SNPs remain after merging C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/anx and C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/bmi summary statistics
#> Results for genetic covariance between: C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/anx and C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/bmi
#> Mean Z*Z: 0.0732
#> Cross trait Intercept: 0.0168 (0.0097)
#> Total Observed Scale Genetic Covariance (g_cov): 0.0121 (0.0043)
#> g_cov Z: 2.85
#> g_cov P-value: 0.0044091
#> Calculating genetic covariance [3/6] for traits: C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/anx and C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/cad
#> 1091859 SNPs remain after merging C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/anx and C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/cad summary statistics
#> Results for genetic covariance between: C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/anx and C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/cad
#> Mean Z*Z: 0.0446
#> Cross trait Intercept: -0.0098 (0.0053)
#> Total Observed Scale Genetic Covariance (g_cov): 0.0245 (0.0038)
#> g_cov Z: 6.53
#> g_cov P-value: 6.7903e-11
#> Estimating heritability [4/6] for: C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/bmi
#> Heritability Results for trait: C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/bmi
#> Mean Chi^2 across remaining SNPs: 3.9344
#> Lambda GC: 2.7869
#> Intercept: 1.0199 (0.0277)
#> Ratio: 0.0068 (0.0094)
#> Total Observed Scale h2: 0.2091 (0.0063)
#> h2 Z: 33.3
#> Calculating genetic covariance [5/6] for traits: C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/bmi and C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/cad
#> 1011369 SNPs remain after merging C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/bmi and C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/cad summary statistics
#> Results for genetic covariance between: C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/bmi and C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/cad
#> Mean Z*Z: 0.2977
#> Cross trait Intercept: 0.0175 (0.0088)
#> Total Observed Scale Genetic Covariance (g_cov): 0.0423 (0.0032)
#> g_cov Z: 13.1
#> g_cov P-value: 5.2908e-39
#> Estimating heritability [6/6] for: C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/cad
#> Please note that the results initially printed to the screen and log file reflect the NON-liability h2 and cov_g. However, a liability conversion is being used for trait C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/cad when creating the genetic covariance matrix used as input for Genomic SEM and liability scale results are printed at the end of the log file.
#> Heritability Results for trait: C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/cad
#> Mean Chi^2 across remaining SNPs: 1.1379
#> Lambda GC: 1.0465
#> Intercept: 0.8515 (0.0092)
#> Ratio: -1.0768 (0.0665)
#> Total Observed Scale h2: 0.0945 (0.0059)
#> h2 Z: 16
#> Liability Scale Results
#> Liability scale results for: ANX
#> Total Liability Scale h2: 0.1108 (0.008)
#> Total Liability Scale Genetic Covariance between ANX and BMI: 0.0125 (0.0044)
#> Total Liability Scale Genetic Covariance between ANX and CAD: 0.0465 (0.0071)
#> Liability scale results for: BMI
#> Total Liability Scale h2: 0.2091 (0.0063)
#> Total Liability Scale Genetic Covariance between BMI and CAD: 0.0781 (0.006)
#> Liability scale results for: CAD
#> Total Liability Scale h2: 0.322 (0.0202)
#> Genetic Correlation Results
#> Genetic Correlation between ANX and BMI: 0.082 (0.0288)
#> Genetic Correlation between ANX and CAD: 0.2461 (0.0377)
#> Genetic Correlation between BMI and CAD: 0.3011 (0.023)
#> LDSC finished running at 2022-06-28 10:14:17
#> Running LDSC for all files took 1 minutes and 3 seconds
  est_genetic_correlation(trait1 = "CAD", trait2 = "ANX", gen_cov = gen_cov)
#> [1] "Running primary model"
#> [1] "Calculating CFI"
#> [1] "Calculating Standardized Results"
#> [1] "Calculating SRMR"
#> elapsed 
#>    0.47 
#> [1] "Model fit statistics are all printed as NA as you have specified a fully saturated model (i.e., df = 0)"
#> # A tibble: 1 × 5
#>   Trait1 Trait2 Adjust_for    rg  rg_se
#>   <chr>  <chr>  <chr>      <dbl>  <dbl>
#> 1 CAD    ANX    Nothing    0.246 0.0382
  est_genetic_correlation(trait1 = "CAD", trait2 = "ANX", gen_cov = gen_cov, adjust = "BMI")
#> [1] "Running primary model"
#> [1] "Calculating CFI"
#> [1] "Calculating Standardized Results"
#> [1] "Calculating SRMR"
#> elapsed 
#>    0.61 
#> [1] "Model fit statistics are all printed as NA as you have specified a fully saturated model (i.e., df = 0)"
#> # A tibble: 1 × 5
#>   Trait1 Trait2 Adjust_for    rg  rg_se
#>   <chr>  <chr>  <chr>      <dbl>  <dbl>
#> 1 CAD    ANX    BMI        0.233 0.0394
```

Finally, we can compute genomic regions enriched for the shared genetic
component of CAD, ANX, and BMI

``` r
  gen_cov_stratified <- est_stratified_gen_cov(paths, sample_prev = sample_prev, population_prev = pop_prev, trait_names = trait_names)
#> Analysis started at2022-06-28 10:14:18
#> The following traits are being analyzed analyzed:ANX
#> 
#>  The following traits are being analyzed analyzed:BMI
#> 
#>  The following traits are being analyzed analyzed:CAD
#> 
#> The following annotations were added to the model: 
#> REF/1000G_Phase3_baselineLD_v2.2_ldscores/
#> 
#> Reading in LD scores from REF/1000G_Phase3_baselineLD_v2.2_ldscores/.[1-22]
#> 
#> LD scores contain1190321SNPs and97annotations
#> 
#> Reading in weighted LD scores fromREF/1000G_Phase3_weights_hm3_no_MHC/.[1-22]
#> 
#> Weighted LD scores contain 1187349 SNPs
#>  
#> Reading in annotation files. This step may take a few minutes.
#> Read in summary statistics [1/3] from: C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/anx
#> Out of 1100641 SNPs, 1087360 remain after merging with LD-score files
#> Removing 0 SNPs with Chi^2 > 83.566; 1087360 remain
#> Read in summary statistics [2/3] from: C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/bmi
#> Out of 1019865 SNPs, 1009250 remain after merging with LD-score files
#> Removing 6 SNPs with Chi^2 > 795.64; 1009244 remain
#> Read in summary statistics [3/3] from: C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/cad
#> Out of 1186097 SNPs, 1172697 remain after merging with LD-score files
#> Removing 23 SNPs with Chi^2 > 154.654; 1172674 remain
#> Liability h2:0.1229(0.0121)
#> 
#> Lambda GC:1.1715
#> Mean Chi^2:1.1835
#> Intercept: 1.0018(0.0113)
#> Ratio: 0.0096(0.0617)
#> Partitioning the heritability over the annotations
#> 946781 SNPs remain after merging C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/anx and C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/bmi summary statistics
#> Results for covariance between:C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/anxandC:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/bmi
#> Mean Z*Z:0.0726
#> Cross trait Intercept: 0.0328(0.0123)
#> cov_g:0.005(0.0057)
#> Partitioning the genetic covariance over the annotations
#> 1083181 SNPs remain after merging C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/anx and C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/cad summary statistics
#> Results for covariance between:C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/anxandC:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/cad
#> Mean Z*Z:0.0448
#> Cross trait Intercept: -0.0039(0.0065)
#> cov_g:0.0212(0.0059)
#> Partitioning the genetic covariance over the annotations
#> h2:0.2168(0.0067)
#> 
#> Lambda GC:2.7872
#> Mean Chi^2:3.9302
#> Intercept: 1.0607(0.0309)
#> Ratio: 0.0207(0.0105)
#> Partitioning the heritability over the annotations
#> 1006198 SNPs remain after merging C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/bmi and C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/cad summary statistics
#> Results for covariance between:C:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/bmiandC:/Users/jacber/AppData/Local/Temp/RtmpMdyWEx/temp_libpath14805ff96614/funcenrich/extdata/cad
#> Mean Z*Z:0.2972
#> Cross trait Intercept: 0.0169(0.0094)
#> cov_g:0.0456(0.0035)
#> Partitioning the genetic covariance over the annotations
#> Liability h2:0.3495(0.022)
#> 
#> Lambda GC:1.0466
#> Mean Chi^2:1.1359
#> Intercept: 0.8518(0.0099)
#> Ratio: -1.0906(0.0731)
#> Partitioning the heritability over the annotations
#> Analysis ended at 2022-06-28 10:23:28
#> Analysis took 9 minutes and 10.4758009910583 seconds
  est_enrichment(trait_names, gen_cov_stratified)
#> [1] "baseL2 is assumed to be the baseline annotation that includes all SNPs."
#> [1] "Running model for baseline annotation"
#> [1] "Confirming fixed model reproduces estimate from freely estimated model for baseline annotation."
#> [1] "Beginning estimation of enrichment for 97 functional annotations."
#> [1] "46 annotations were removed from the output because they were either continuous or flanking window annotatoins."
#> # A tibble: 51 × 11
#>    Annotation     lhs   op    rhs   Cov_Smooth Z_smooth Enrichment Enrichment_SE
#>    <chr>          <chr> <chr> <chr>      <dbl>    <dbl>      <dbl>         <dbl>
#>  1 Conserved_Mam… F     ~~    F              0        0      9.92         3.11  
#>  2 Conserved_Pri… F     ~~    F              0        0      9.45         3.36  
#>  3 H3K9ac_peaks_… F     ~~    F              0        0      4.97         1.91  
#>  4 MAFbin9L2      F     ~~    F              0        0      1.54         0.334 
#>  5 baseL2         F     ~~    F              0        0      1.00         0.0740
#>  6 Intron_UCSCL2  F     ~~    F              0        0      1.01         0.120 
#>  7 Transcr_Hoffm… F     ~~    F              0        0      1.22         0.261 
#>  8 SuperEnhancer… F     ~~    F              0        0      1.10         0.203 
#>  9 H3K4me1_Trynk… F     ~~    F              0        0      1.18         0.255 
#> 10 H3K27ac_Hnisz… F     ~~    F              0        0      0.920        0.159 
#> # … with 41 more rows, and 3 more variables: Enrichment_p_value <dbl>,
#> #   Error <dbl>, Warning <chr>
```
