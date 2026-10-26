#' Scissor: Single-Cell Identification of Subpopulations with bulk Sample phenOtype coRrelation
#'
#' \code{Scissor} is a novel approach that utilizes the phenotypes, such as disease stage, tumor metastasis, treatment response, and survival
#' outcomes, collected from bulk assays to identify the most highly phenotype-associated cell subpopulations from single-cell data.
#'
#' Scissor is a novel algorithm to identify cell subpopulations from single-cell data that are most highly associated with the given phenotypes.
#' The three data sources for Scissor inputs are a single-cell expression matrix, a bulk expression matrix, and a phenotype of interest.
#' The phenotype annotation of each bulk sample can be a continuous dependent variable, binary group indicator vector, or clinical survival data.
#' The key step of Scissor is to quantify the similarity between the single-cell and bulk samples by Pearson correlation for each pair of cells and bulk samples.
#' After this, Scissor optimizes a regression model on the correlation matrix with the sample phenotype.
#' The selection of the regression model depends on the type of the input phenotype, i.e., linear regression for continuous variables,
#' logistic regression for dichotomous variables, and Cox regression for clinical survival data.
#' Based on the signs of the estimated regression coefficients, the cells with non-zero coefficients can be indicated as
#' Scissor positive (Scissor+) cells and Scissor negative (Scissor-) cells, which are positively and negatively associated
#' with the phenotype of interest, respectively.
#'
#' @param bulk_dataset Bulk expression matrix of related disease. Each row represents a gene and each column represents a sample.
#' Counts, TPM or CPM all work, since the data are quantile-normalized together with the single-cell data. Genes expressed in fewer
#' than 10 samples and samples with fewer than 10 expressed genes are removed.
#' @param sc_dataset A Seurat (v5) object with raw counts. It is preprocessed with \code{\link{sc_preprocess}} (HVG selection,
#' optional Harmony batch correction, cell-cell graph construction) on every run without \code{Load_file}.
#' @param phenotype Phenotype annotation of each bulk sample, named by bulk sample. It can be a continuous dependent variable,
#' binary group indicator vector, or clinical survival data:
#'   \itemize{
#'   \item Continuous dependent variable. Should be a named quantitative vector for \code{family = gaussian}.
#'   \item Binary group indicator vector. Should be a named 0-1 encoded vector or a factor with two levels for \code{family = binomial}.
#'   For a factor, the first level is coded 0 and the second 1.
#'   \item Clinical survival data. Should be a two-column data frame or matrix with bulk sample names as rownames and columns named
#'   'time' and 'status'. The latter is a binary variable, with '1' indicating event (e.g.recurrence of cancer or death), and '0'
#'   indicating right censored.
#'   }
#' @param tag Names for each phenotypic group. Required for logistic regression. For linear regression, give one tag per distinct
#' phenotype value (e.g. disease stages), or leave \code{NULL} for a continuous phenotype.
#' @param alpha Parameter used to balance the effect of the l1 norm and the network-based penalties. It can be a number or a searching vector.
#' If \code{alpha = NULL}, a default searching vector is used. The range of alpha is in \code{[0,1]}. A larger alpha lays more emphasis on the l1 norm.
#' @param cutoff Cutoff for the percentage of the Scissor selected cells in total cells. This parameter is used to restrict the number of the
#' Scissor selected cells. A cutoff less than \code{50\%} (default \code{20\%}) is recommended depending on the input data.
#' @param family Response type for the regression model. It depends on the type of the given phenotype and
#' can be \code{family = gaussian} for linear regression, \code{family = binomial} for classification, or \code{family = cox} for Cox regression.
#' @param sc_correct_batch Logical. Run Harmony batch correction on \code{sc_dataset}. Only used if \code{sc_batch} is found in its metadata.
#' @param sc_batch Name of the \code{sc_dataset} metadata column giving the batch. If \code{NULL}, no batch correction is done.
#' @param only_hvg Logical. Use only the highly variable genes of \code{sc_dataset} (shared with the bulk data) to correlate cells and bulk samples.
#' If \code{FALSE}, all shared genes are used.
#' @param n_hvg Number of highly variable genes to select in \code{sc_dataset}.
#' @param Save_file File name for saving the preprocessed regression inputs into an RDS file.
#' @param Load_file File name of an RDS file saved by a previous run (\code{Save_file}). Skips preprocessing; it can help to tune the
#' model parameter \code{alpha}. \code{.RData} files saved by the original Scissor cannot be loaded.
#' Please see Scissor Tutorial for more details.
#'
#' @return This function returns a list with the following components:
#'   \item{para}{A list contains the final model parameters.}
#'   \item{Coefs}{The regression coefficient for each cell.}
#'   \item{Scissor_pos}{The cell IDs of Scissor+ cells.}
#'   \item{Scissor_neg}{The cell IDs of Scissor- cells.}
#'
#' @references Duanchen Sun and Zheng Xia (2021): Phenotype-guided subpopulation identification from single-cell sequencing data. Nature Biotechnology.
#' @import Seurat Matrix preprocessCore
#' @export
Scissor <- function(bulk_dataset, sc_dataset, phenotype, tag = NULL,
                    alpha = NULL, cutoff = 0.2, family = c("gaussian","binomial","cox"),
                    sc_correct_batch = TRUE, sc_batch = NULL, only_hvg = TRUE, n_hvg = 4000, ### MOD: add mew parameters for sc_preprocess ### 
                    Save_file = "Scissor_inputs.rds", Load_file = NULL){                     ### MOD: save RDS file rather than Rda ### 
    
    family = match.arg(family) ### MOD: defaults to "gaussian" if not given ###
    library(Seurat)
    library(Matrix)
    library(preprocessCore)


    if (is.null(Load_file)){
        bulk_dataset = as.matrix(bulk_dataset)

        ### MOD: Preprocess the sc_dataset with batch correction ###
        if(sc_correct_batch && !is.null(sc_batch) && sc_batch %in% colnames(sc_dataset@meta.data)){
            sc_dataset = sc_preprocess(sc_dataset, assay=DefaultAssay(sc_dataset), correct_batch=sc_correct_batch, group.by.vars=sc_batch, nfeatures_batch=n_hvg+1000, nfeatures=n_hvg)
        } else{
            message("No batch correction. Possible reasons: (1) sc_correct_batch set to FALSE; (2) cannot find ",sc_batch," in sc_dataset metadata.")
            sc_dataset = sc_preprocess(sc_dataset, assay=DefaultAssay(sc_dataset), correct_batch=F, nfeatures_batch=n_hvg+1000, nfeatures=n_hvg)
        }

        ### MOD: Keep the cell-cell graph sparse ###
        network = as(sc_dataset@graphs[[paste0(DefaultAssay(sc_dataset), "_snn")]],"dgCMatrix")
        diag(network) = 0
        network = drop0(network)

        ### MOD: Basic QC on bulk expression matrix ###
        bulk_dataset = .qc_on_mtx(bulk_dataset,row_min_cnt = 10,col_min_cnt = 10)

        ### MOD: Use all features or only the highly variable features for integrating single-cell and bulk data ###
        if(only_hvg){
            common = intersect(rownames(bulk_dataset), VariableFeatures(sc_dataset))
        } else{
            common <- intersect(rownames(bulk_dataset), rownames(sc_dataset))
        }
        
        if (length(common) == 0) {
            stop("There is no common genes between the given single-cell and bulk samples.")
        }

        # Align the gene and sample lists
        bulk_dataset = bulk_dataset[common,]
        sc_exprs = as.matrix(sc_dataset[[DefaultAssay(sc_dataset)]]$counts[common,])
        pheno_samples = if(family == "cox") rownames(phenotype) else names(phenotype)
        if(is.null(pheno_samples)){
            stop("phenotype must be named by bulk sample: rownames for family = 'cox', names for 'gaussian'/'binomial'.")
        }
        samples = intersect(pheno_samples, colnames(bulk_dataset))
        bulk_dataset = bulk_dataset[,samples]
        phenotype = if(family == "cox") phenotype[samples,] else phenotype[samples]
        
        # if (class(sc_dataset) == "Seurat"){
        #     sc_exprs <- as.matrix(sc_dataset@assays$RNA@data)
        #     network  <- as.matrix(sc_dataset@graphs$RNA_snn)
        # }else{
        #     sc_exprs <- as.matrix(sc_dataset)
        #     Seurat_tmp <- CreateSeuratObject(sc_dataset)
        #     Seurat_tmp <- FindVariableFeatures(Seurat_tmp, selection.method = "vst", verbose = F)
        #     Seurat_tmp <- ScaleData(Seurat_tmp, verbose = F)
        #     Seurat_tmp <- RunPCA(Seurat_tmp, features = VariableFeatures(Seurat_tmp), verbose = F)
        #     Seurat_tmp <- FindNeighbors(Seurat_tmp, dims = 1:10, verbose = F)
        #     network  <- as.matrix(Seurat_tmp@graphs$RNA_snn)
        # }
        # diag(network) <- 0
        # network[which(network != 0)] <- 1
        #############################################################################################################

        dataset0 <- cbind(bulk_dataset[common,], sc_exprs[common,])         # Dataset before quantile normalization.
        dataset1 <- normalize.quantiles(dataset0)                           # Dataset after  quantile normalization.
        rownames(dataset1) <- rownames(dataset0)
        colnames(dataset1) <- colnames(dataset0)

        Expression_bulk <- dataset1[,1:ncol(bulk_dataset)]
        Expression_cell <- dataset1[,(ncol(bulk_dataset) + 1):ncol(dataset1)]
        X <- cor(Expression_bulk, Expression_cell)

        quality_check <- quantile(X)
        print("|**************************************************|")
        print("Performing quality-check for the correlations")
        print("The five-number summary of correlations:")
        print(quality_check)
        print("|**************************************************|")
        if (quality_check[3] < 0.01){
            warning("The median correlation between the single-cell and bulk samples is relatively low.")
        }
        if (family == "binomial"){
            ##### MOD: code factor/character phenotypes as 0/1 #####
            # Y <- as.numeric(phenotype)
            if (is.factor(phenotype) || is.character(phenotype)){
                phenotype <- factor(phenotype)
                Y <- as.numeric(phenotype) - 1   # first level -> 0, second level -> 1
            }else{
                Y <- as.numeric(phenotype)
            }
            if (!all(Y %in% c(0, 1))){
                stop("For family = 'binomial', phenotype must be 0/1 or a factor with two levels.")
            }
            ########################################################
            z <- table(Y)
            if (length(z) != length(tag)){
                stop("The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type.")
            }else{
                print(sprintf("Current phenotype contains %d %s and %d %s samples.", z[1], tag[1], z[2], tag[2]))
                print("Perform logistic regression on the given phenotypes:")
            }
        }
        if (family == "gaussian"){
            Y <- as.numeric(phenotype)
            ### MOD: allow continuous phenotypes; tags are only needed for grouped (e.g. stage) phenotypes ###
            # z <- table(Y)
            # if (length(z) != length(tag)){
            #     stop("The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type.")
            # }else{
            #     tmp <- paste(z, tag)
            #     print(paste0("Current phenotype contains ", paste(tmp[1:(length(z)-1)], collapse = ", "), ", and ", tmp[length(z)], " samples."))
            #     print("Perform linear regression on the given phenotypes:")
            # }
            if (any(is.na(Y))){
                stop("For family = 'gaussian', phenotype must be numeric without missing values.")
            }
            if (is.null(tag)){
                # Continuous phenotype: no groups to name
                print(sprintf("Current phenotype is continuous across %d samples (min = %s, median = %s, max = %s).",
                              length(Y), signif(min(Y), 3), signif(median(Y), 3), signif(max(Y), 3)))
            }else{
                # Grouped phenotype: one tag per distinct value, as in the original Scissor
                z <- table(Y)
                if (length(z) != length(tag)){
                    stop("The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type.")
                }
                tmp <- paste(z, tag)
                print(paste0("Current phenotype contains ", paste(tmp[1:(length(z)-1)], collapse = ", "), ", and ", tmp[length(z)], " samples."))
            }
            print("Perform linear regression on the given phenotypes:")
            ####################################################################################################
        }
        if (family == "cox"){
            Y <- as.matrix(phenotype)
            if (ncol(Y) != 2){
                stop("The size of survival data is wrong. Please check Scissor inputs and selected regression type.")
            }else{
                print("Perform cox regression on the given clinical outcomes:")
            }
        }

        ### MOD: Save the intermediate results for the first run ####
        # save(X, Y, network, Expression_bulk, Expression_cell, file = Save_file)
        res = list(X = X, Y = Y, network = network, Expression_bulk = Expression_bulk, Expression_cell = Expression_cell)
        saveRDS(res,file = Save_file) 
        #############################################################

    }else{
        ###### MOD: Load the preprocessed data ######
        # load(Load_file)
        res = readRDS(Load_file)
        X = res$X
        Y = res$Y
        network = res$network
        Expression_bulk = res$Expression_bulk
        Expression_cell = res$Expression_cell
        rm(res)
        print('Loading finished!')
        ##############################################
    }

    if (is.null(alpha)){
        alpha <- c(0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
    }
    for (i in 1:length(alpha)){
        set.seed(123)
        fit0 <- APML1(X, Y, family = family, penalty = "Net", alpha = alpha[i], Omega = network, nlambda = 100, nfolds = min(10,nrow(X)))
        fit1 <- APML1(X, Y, family = family, penalty = "Net", alpha = alpha[i], Omega = network, lambda = fit0$lambda.min)
        if (family == "binomial"){
            Coefs <- as.numeric(fit1$Beta[2:(ncol(X)+1)])
        }else{
            Coefs <- as.numeric(fit1$Beta)
        }
        Cell1 <- colnames(X)[which(Coefs > 0)]
        Cell2 <- colnames(X)[which(Coefs < 0)]
        percentage <- (length(Cell1) + length(Cell2)) / ncol(X)
        print(sprintf("alpha = %s", alpha[i]))
        print(sprintf("Scissor identified %d Scissor+ cells and %d Scissor- cells.", length(Cell1), length(Cell2)))
        print(sprintf("The percentage of selected cell is: %s%%", formatC(percentage*100, format = 'f', digits = 3)))

        if (percentage < cutoff){
            break
        }
        cat("\n")
    }
    print("|**************************************************|")

    return(list(para = list(alpha = alpha[i], lambda = fit0$lambda.min, family = family),
                Coefs = Coefs,
                Scissor_pos = Cell1,
                Scissor_neg = Cell2))
}

