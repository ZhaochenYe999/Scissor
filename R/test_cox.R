test_cox <- function(X, Y, network, alpha, cell_num, n = 100, nfold = 10,
                     bkg_stat=NULL){ # MOD: parallelizing the permutations

    set.seed(1)
    m1 <- sum(Y[,2] == 1)
    m2 <- sum(Y[,2] == 0)
    index0 <- sample(cut(seq(m1+m2), breaks = nfold, labels = F))
    #index1 <- sample(cut(seq(m1), breaks = nfold, labels = F))
    #index2 <- sample(cut(seq(m2), breaks = nfold, labels = F))

    print("|**************************************************|")
    print("Perform cross-validation on X with true label")
    c_index_test_real <- NULL
    pb1 <- progress_bar$new(total = nfold)
    for (j in 1:nfold){
        c_index <- which(index0 == j)
        #c_index <- c(which(Y[,2] == 1)[which(index1 == j)], which(Y[,2] == 0)[which(index2 == j)])
        X_train <- X[-c_index,]
        Y_train <- Y[-c_index,]
        fit <- NULL
        while (is.null(fit$fit)){
            set.seed(123)
            fit <- APML1(X_train, Y_train, family = "cox", penalty = "Net", alpha = alpha, Omega = network, nlambda = 100)
        }
        index <- which.min(abs(fit$fit$nzero - cell_num))
        Coefs <- as.numeric(fit$Beta[,index])
        Cell1 <- Coefs[which(Coefs > 0)]
        Cell2 <- Coefs[which(Coefs < 0)]

        X_test <- X[c_index,]
        Y_test <- Y[c_index,]
        test_data <- data.frame(cbind(Y_test, X_test%*%Coefs))
        colnames(test_data) <- c("OS_time", "Status", "Prediction")
        res.cox <- coxph(Surv(OS_time, Status) ~ Prediction, data = test_data)
        c_index_test_real[j] <- concordance(res.cox)$concordance

        #pb1$tick()
        Sys.sleep(1 / 100)
        if (j == nfold) cat("Finished!\n")
    }

    ###################### MODIFIED ############################
    if(is.null(bkg_stat)){
        print("|**************************************************|")
        print("Perform cross-validation on X with permutated label")
        # c_index_test_back <- list()
        # pb2 <- progress_bar$new(total = n)
        # for (i in 1:n){
        #     set.seed(i+100)
        #     c_index_test_back[[i]] <- matrix(0, nfold, 1, dimnames = list(paste0("Testing_", 1:nfold),  "Concordance"))
        #     Y2 <- Y[sample(nrow(Y)),]
        #     for (j in 1:nfold){
        #         c_index <- which(index0 == j)
        #         #c_index <- c(which(Y2[,2] == 1)[which(index1 == j)], which(Y2[,2] == 0)[which(index2 == j)])
        #         X_train <- X[-c_index,]
        #         Y_train <- Y2[-c_index,]
        #         fit <- NULL
        #         while (is.null(fit$fit)){
        #             set.seed(123)
        #             fit <- APML1(X_train, Y_train, family = "cox", penalty = "Net", alpha = alpha, Omega = network, nlambda = 100)
        #         }
        #         index <- which.min(abs(fit$fit$nzero - cell_num))
        #         Coefs <- as.numeric(fit$Beta[,index])
        #         Cell1 <- Coefs[which(Coefs > 0)]
        #         Cell2 <- Coefs[which(Coefs < 0)]

        #         X_test <- X[c_index,]
        #         Y_test <- Y2[c_index,]
        #         test_data <- data.frame(cbind(Y_test, X_test%*%Coefs))
        #         colnames(test_data) <- c("OS_time", "Status", "Prediction")
        #         res.cox <- coxph(Surv(OS_time, Status) ~ Prediction, data = test_data)
        #         c_index_test_back[[i]][j] <- concordance(res.cox)$concordance
        #     }
        #     #pb2$tick()
        #     Sys.sleep(1 / 100)
        #     if (i == n) cat("Finished!\n")

        plan(multisession, workers = 8)
        # future_lapply replaces the 'for (i in 1:n)' loop
        c_index_test_back <- future_lapply(1:n, function(i) {
            set.seed(i + 100)
            Y2 <- Y[sample(nrow(Y)), ] # Permute labels
        
            # Local nfold loop inside each permutation
            fold_results <- numeric(nfold)
            for (j in 1:nfold) {
                c_index <- which(index0 == j)
                X_train <- X[-c_index, ]; Y_train <- Y2[-c_index, ]
            
                fit <- NULL
                while (is.null(fit$fit)) {
                    set.seed(123)
                    fit <- APML1_2025(X_train, Y_train, family = "cox", penalty = "Net", 
                                      alpha = alpha, Omega = network, nlambda = 100)
                }
                index <- which.min(abs(fit$fit$nzero - cell_num))
                Coefs <- as.numeric(fit$Beta[, index])
                #Cell1 <- Coefs[which(Coefs > 0)]
                #Cell2 <- Coefs[which(Coefs < 0)]
            
                test_data <- data.frame(cbind(Y2[c_index,], X[c_index,] %*% Coefs))
                colnames(test_data) <- c("OS_time", "Status", "Prediction")
                res.cox <- survival::coxph(survival::Surv(OS_time, Status) ~ Prediction, data = test_data)
                fold_results[j] <- survival::concordance(res.cox)$concordance
            }
            return(fold_results)
        }, future.seed = TRUE) # CRITICAL for correct random numbers in parallel

        background <- sapply(c_index_test_back, mean) # MOD
    } else{
        background = bkg_stats
        c_index_test_back = NA
    }
    statistic  <- mean(c_index_test_real)
    # background <- NULL
    # for (i in 1:n){
    #     background[i] <- mean(c_index_test_back[[i]][,1])
    # }
    ############################################################
    
    p <- sum(background > statistic)/n

    print(sprintf("Test statistic = %s", formatC(statistic, format = "f", digits = 8))) # MOD: Increase digits
    print(sprintf("Reliability significance test p = %s", formatC(p, format = "f", digits = 8))) # MOD: Increase digits

    return(list(statistic = statistic,
                p = p,
                c_index_test_real = c_index_test_real,
                c_index_test_back = c_index_test_back))
}




