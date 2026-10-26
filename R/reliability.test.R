#' Reliability Significance Test
#'
#' Performs the reliability significance test to determine whether the inferred phenotype-to-cell
#' associations are reliable (statistical p-value less than 0.05) or are false positives.
#'
#' This statistical test can determine whether the inferred phenotype-to-cell associations are
#' reliable (statistical p-value less than 0.05) or are false positives.
#'
#' The motivation for the reliability significance test is: if the chosen single-cell and bulk data are not suitable for
#' the phenotype-to-cell associations, the correlations would be less informative and not well associated with the phenotype
#' labels. Thus, the corresponding prediction performance would be poor and not be significantly distinguishable from the
#' randomly permutated labels.
#'
#' The evaluation measurements used in the reliability significance test are the mean squared error (MSE) for linear regression,
#' the area under the ROC curve (AUC) for classification, and the concordance index (c-index) for Cox regression.
#'
#' @param X Scissor calculated correlation matrix for each pair of cells and bulk samples.
#' @param Y Scissor calculated response variable in the regression model.
#' @param network Scissor calculated cell-cell similarity network.
#' @param alpha Same parameter alpha used in Scissor.
#' @param family Same parameter family used in Scissor.
#' @param cell_num The number of the Scissor selected cells.
#' @param n Permutation times.
#' @param nfold The fold number in cross-validation.
#' @param bkg_stat Optional precomputed background statistics, one mean measurement per permutation, e.g. the
#' \code{background} element of a previous result. If supplied, the permutations are skipped.
#' @param workers Number of parallel workers for the permutations. Defaults to 8, or fewer if fewer
#' cores are available (respects Slurm allocations).
#'
#' @return A list containing the following components:
#'   \item{statistic}{The test statistic.}
#'   \item{p}{The test p-value.}
#'   \item{Measurement_test_real}{The evaluation measurement on each fold calculated using the real label. Named \code{MSE_test_real},
#'   \code{AUC_test_real} or \code{c_index_test_real} for \code{family = gaussian}, \code{binomial} or \code{cox}.}
#'   \item{Measurement_test_back}{A list with each component contains the evaluation measurements on each fold calculated using the
#'   permutated labels (\code{MSE_test_back}, \code{AUC_test_back} or \code{c_index_test_back}). \code{NA} if \code{bkg_stat} was supplied.}
#'   \item{background}{The mean measurement of each permutation. Can be passed as \code{bkg_stat} to rerun the test without permutations.}
#'
#' @import progress Matrix
#' @export
reliability.test <- function(X, Y, network, alpha, family = c("gaussian","binomial","cox"), cell_num, n = 100, nfold = 10,
                             bkg_stat = NULL, workers = min(8, future::availableCores())){ ### MOD: Pass precomputed background ###

    family <- match.arg(family) ### MOD: Defaults to "gaussian" if not given ###
    library(progress)
    library(Matrix)

    if (family == 'gaussian'){
        result <- test_lm(X, Y, network, alpha, cell_num, n, nfold, bkg_stat, workers)
    }
    if (family == 'binomial'){
        library(pROC)
        result <- test_logit(X, Y, network, alpha, cell_num, n, nfold, bkg_stat, workers)
    }
    if (family == 'cox'){
        library(survival)
        result <- test_cox(X, Y, network, alpha, cell_num, n, nfold, bkg_stat, workers)
    }
    return(result)
}
