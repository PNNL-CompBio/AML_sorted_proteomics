mDEG <- function(data.list, factor.info,
                 feature.names = rep("Gene", length(data.list)), p = 0.05, 
                 FDR.features = 0.05, n.dot.features = 10) {
  #### Step 1. Check if formats are correct ####
  # check that there are as many types as data.list inputs
  if (length(feature.names) != length(data.list)) {
    stop("Length of feature.names vector must match that of data.list")
  }
  types <- names(data.list)

  #### Step 2. Differential expression analysis ####
  if (length(levels(factor.info[,1])) != 2) {
    stop("Exactly 2 levels are required for the first column of factor.info")
  } else {
    # use parallel computing if possible
    if (requireNamespace("parallel") &
        requireNamespace("snow") &
        requireNamespace("doSNOW")) {
      cores <- parallel::detectCores() # number of cores available
      if (cores[1] > 1) {
        cl <- snow::makeCluster(cores[1] - 1) # cluster using all but 1 core
        doSNOW::registerDoSNOW(cl) # register cluster
      }
    }
    
    deg <- list()
    for (i in 1:length(types)) {
      # make sure there are no duplicated feature names
      if (nrow(data.list[[i]]) > length(unique(rownames(data.list[[i]])))) {
        stop(paste(
          "There must be only one expression value per sample",
          "for each feature"
        ))
      } else {
        # select annotated samples and set feature names as rownames
        all.data.list <- data.list[[i]][,rownames(factor.info)]
        
        rownames(all.data.list) <- data.list[[i]][ , feature.names[i]]

        # create expression sets
        eset <- Biobase::ExpressionSet(as.matrix(all.data.list))

        # log2 transformation if distribution isn't normal
        ex <- Biobase::exprs(eset)
        qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0),
          na.rm = TRUE
        ))
        LogC <- (qx[5] > 100) ||
          (qx[6] - qx[1] > 50 && qx[2] > 0)
        if (LogC) {
          ex[which(ex <= 0)] <- NaN
          Biobase::exprs(eset) <- log2(ex)
        }

        # identify sample phenotypes and set up design matrix
        eset$group <- factor.info[,1]
        #levels(eset$group) <- levels(factor.info[,1])
        if (ncol(factor.info) > 1) {
          eset$batch <- factor.info[,2]
          levels(eset$batch) <- levels(factor.info[,2])
          design <- stats::model.matrix(~ group + batch + 0, eset)
          colnames(design)[1:2] <- substr(colnames(design)[1:2], 6, nchar(colnames(design)[1:2]))
        } else {
          design <- stats::model.matrix(~ group + 0, eset)
          colnames(design) <- substr(colnames(design), 6, nchar(colnames(design)))
        }
        
        # fit linear model
        fit <- limma::lmFit(eset, design)

        # set up contrasts of interest and redo fit
        cts <- paste(levels(factor.info[,1]), collapse = "-")
        cont.matrix <- limma::makeContrasts(contrasts = cts, levels = design)
        fit <- limma::contrasts.fit(fit, cont.matrix)

        # calculate differential expression and statistics
        fit <- limma::eBayes(fit, 0.01)
        deg[[types[i]]] <- limma::topTable(fit,
          adjust = "fdr",
          number = nrow(fit)
        )
        
        deg[[types[i]]][, feature.names[i]] <- rownames(deg[[types[i]]])
        colnames(deg[[types[i]]])[1] <- "Log2FC"
        deg[[types[i]]] <-
          deg[[types[i]]][, c(
            feature.names[i],
            colnames(deg[[types[i]]])
            [1:(ncol(deg[[types[i]]]) - 1)]
          )]
        
        # add log2transform info as last column
        deg[[types[i]]]$Log2Transformed <- LogC
      }
    }
    
    # compile DEG results across omics types if less than all of them are phospho
    if (any(grepl("phospho", names(deg), ignore.case = TRUE))) {
      compileDEGs <- all(grepl("phospho", names(deg), ignore.case = TRUE))
    } else {
      compileDEGs <- TRUE
    }
    if (length(types) > 1 & compileDEGs) {
      compiled.DEGs <- panSEA::compile_mDEG(deg, p, FDR.features, 
                                            n.dot.features)
    } else {
      compiled.DEGs <- NA
    }
    
    # stop cluster if relevant
    if (requireNamespace("parallel") &
        requireNamespace("snow") &
        requireNamespace("doSNOW")) {
      if (cores[1] > 1) {
        snow::stopCluster(cl) # stop cluster
        rm(cl)
      }
    }

  return(list(compiled.results = compiled.DEGs, all.results = deg
  ))
  }
}


