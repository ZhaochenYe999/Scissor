test_logit <- function(X, Y, network, alpha, cell_num, n = 100, nfold = 10,
                       bkg_stat=NULL, workers=min(8, future::availableCores())){ ### MOD: parallelizing the permutations ###

    set.seed(2)
    m1 <- sum(Y == 1)
    m2 <- sum(Y == 0)
    index1 <- sample(cut(seq(m1), breaks = nfold, labels = F))
    index2 <- sample(cut(seq(m2), breaks = nfold, labels = F))

    print("|**************************************************|")
    print("Perform cross-validation on X with true label")
    AUC_test_real <- NULL
    pb1 <- progress_bar$new(total = nfold)
    for (j in 1:nfold){
        c_index <- c(which(Y == 1)[which(index1 == j)], which(Y == 0)[which(index2 == j)])
        X_train <- X[-c_index,]
        Y_train <- Y[-c_index]

        ####### MOD: single fit, no randomness when nfolds = 1, so while could only loop forever ######
        # fit <- NULL
        # while (is.null(fit$fit)){
        #     set.seed(123)
        #     fit <- APML1(X_train, Y_train, family = "binomial", penalty = "Net", alpha = alpha, Omega = network, nlambda = 100)
        # }
        fit <- APML1(X_train, Y_train, family = "binomial", penalty = "Net", alpha = alpha, Omega = network, nlambda = 100)
        if (is.null(fit$fit)) stop("APML1 failed to fit fold ", j, " with the true labels.")
        ################################################################################################

        index <- which.min(abs(fit$fit$nzero - cell_num))
        Coefs <- as.numeric(fit$Beta[2:(ncol(X_train)+1), index])
        Cell1 <- Coefs[which(Coefs > 0)]
        Cell2 <- Coefs[which(Coefs < 0)]

        X_test <- X[c_index,]
        Y_test <- Y[c_index]
        score_test <- 1/(1+exp(-X_test%*%Coefs-fit$Beta[1,index]))[,1]
        AUC_test_real[j] <- roc(Y_test, score_test, direction = "<", quiet = T)$auc

        #pb1$tick()
        Sys.sleep(1 / 100)
        if (j == nfold) cat("Finished!\n")
    }

    ###################### MODIFIED: parallelized permutations ##################
    # print("|**************************************************|")
    # print("Perform cross-validation on X with permutated label")
    # AUC_test_back <- list()
    # pb2 <- progress_bar$new(total = n)
    # for (i in 1:n){
    #     set.seed(i+100)
    #     AUC_test_back[[i]] <- matrix(0, nfold, 1, dimnames = list(paste0("Testing_",  1:nfold), "AUC"))
    #     Y2 <- sample(Y)
    #     names(Y2) <- rownames(X)
    #     for (j in 1:nfold){
    #         c_index <- c(which(Y2 == 1)[which(index1 == j)], which(Y2 == 0)[which(index2 == j)])
    #         X_train <- X[-c_index,]
    #         Y_train <- Y2[-c_index]
    #         fit <- NULL
    #         while (is.null(fit$fit)){
    #             set.seed(123)
    #             fit <- APML1(X_train, Y_train, family = "binomial", penalty = "Net", alpha = alpha, Omega = network, nlambda = 100)
    #         }
    #         index <- which.min(abs(fit$fit$nzero - cell_num))
    #         Coefs <- as.numeric(fit$Beta[2:(ncol(X_train)+1), index])
    #         Cell1 <- Coefs[which(Coefs > 0)]
    #         Cell2 <- Coefs[which(Coefs < 0)]
    #
    #         X_test <- X[c_index,]
    #         Y_test <- Y2[c_index]
    #         score_test <- 1/(1+exp(-X_test%*%Coefs-fit$Beta[1,index]))[,1]
    #         AUC_test_back[[i]][j] <- roc(Y_test, score_test, direction = "<", quiet = T)$auc
    #     }
    #     #pb2$tick()
    #     Sys.sleep(1 / 100)
    #     if (i == n) cat("Finished!\n")
    # }
    # statistic  <- mean(AUC_test_real)
    # background <- NULL
    # for (i in 1:n){
    #     background[i] <- mean(AUC_test_back[[i]][,1])
    # }
    # p <- sum(background > statistic)/n
    #
    # print(sprintf("Test statistic = %s", formatC(statistic, format = "f", digits = 3)))
    # print(sprintf("Reliability significance test p = %s", formatC(p, format = "f", digits = 3)))

    if(is.null(bkg_stat)){
        print("|**************************************************|")
        print("Perform cross-validation on X with permutated label")

        old_plan <- future::plan(future::multisession, workers = workers)
        on.exit(future::plan(old_plan), add = TRUE) # restore the caller's plan
        # future_lapply replaces the 'for (i in 1:n)' loop
        AUC_test_back <- future.apply::future_lapply(1:n, function(i) {
            set.seed(i + 100)
            Y2 <- sample(Y) # Permute labels
            names(Y2) <- rownames(X)

            # Local nfold loop inside each permutation
            fold_results <- numeric(nfold)
            for (j in 1:nfold) {
                c_index <- c(which(Y2 == 1)[which(index1 == j)], which(Y2 == 0)[which(index2 == j)])
                X_train <- X[-c_index, ]
                Y_train <- Y2[-c_index]

                fit <- APML1(X_train, Y_train, family = "binomial", penalty = "Net",
                             alpha = alpha, Omega = network, nlambda = 100)
                if (is.null(fit$fit)) {
                    warning("APML1 failed on permutation ", i, ", fold ", j, "; recorded as NA")
                    fold_results[j] <- NA
                    next
                }
                index <- which.min(abs(fit$fit$nzero - cell_num))
                Coefs <- as.numeric(fit$Beta[2:(ncol(X_train)+1), index])

                X_test <- X[c_index, ]
                Y_test <- Y2[c_index]
                score_test <- 1/(1+exp(-X_test%*%Coefs-fit$Beta[1,index]))[,1]
                fold_results[j] <- pROC::roc(Y_test, score_test, direction = "<", quiet = T)$auc
            }
            return(fold_results)
        }, future.seed = TRUE) # CRITICAL for correct random numbers in parallel

        background <- sapply(AUC_test_back, mean, na.rm = TRUE) # skip failed folds
        background <- background[!is.na(background)] # drop permutations where every fold failed
    } else{
        # Reuse a previously computed background (one mean AUC per permutation)
        background = bkg_stat
        AUC_test_back = NA
    }
    statistic  <- mean(AUC_test_real)
    p <- sum(background > statistic)/length(background) # bkg_stat may not have length n

    print(sprintf("Test statistic = %s", formatC(statistic, format = "f", digits = 8)))          ### MOD: Increase digits
    print(sprintf("Reliability significance test p = %s", formatC(p, format = "f", digits = 8))) ### MOD: Increase digits
    ############################################################################

    return(list(statistic = statistic,
                p = p,
                AUC_test_real = AUC_test_real,
                AUC_test_back = AUC_test_back,
                background = background)) ### MOD: one mean per permutation; pass back as bkg_stat ###
}


