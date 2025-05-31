mDEG <- function(data.list, types, group.names = c("Diseased", "Healthy"),
                 group.samples = list(
                   2:(0.5 * (ncol(data.list[[1]]) + 1)),
                   (0.5 * (ncol(data.list[[1]]) + 1) + 1):ncol(data.list[[1]])
                 ), feature.names = rep("Gene", length(types)), p = 0.05, 
                 FDR.features = 0.05, n.dot.features = 10) {
  #### Step 1. Check if formats are correct ####
  # check that there are as many types as data.list inputs
  if (length(types) != length(data.list)) {
    stop("Length of types vector must match that of data.list")
  }

  #### Step 2. Differential expression analysis ####
  if (length(group.names) != 2) {
    stop("Only 2 groups are allowed for differential expression analysis")
  } else {
    deg <- list()
    for (i in 1:length(types)) {
      # make sure there are no duplicated feature names
      if (nrow(data.list[[i]]) > length(unique(rownames(data.list[[i]])))) {
        stop(paste(
          "There must be only one expression value per sample",
          "for each feature"
        ))
      } else {
        # separate data.list based on group (i.e., categorical phenotype)
        data.list1 <- as.data.frame(data.list[[i]][, group.samples[[1]]])
        data.list2 <- as.data.frame(data.list[[i]][, group.samples[[2]]])
        all.data.list <- cbind(data.list1, data.list2)
        
        # set feature names as rownames
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

        # identify sample phenotypes
        gs <- factor(c(rep(0, ncol(data.list1)), rep(1, ncol(data.list2))))
        levels(gs) <- make.names(group.names)
        eset$group <- gs

        # set up design matrix
        design <- stats::model.matrix(~ group + 0, eset)
        colnames(design) <- levels(gs)

        # fit linear model
        fit <- limma::lmFit(eset, design)

        # set up contrasts of interest and redo fit
        cts <- paste(group.names[1], group.names[2], sep = "-")
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
    
    # compile DEG results across omics types
    if (length(types) > 1) {
      compiled.DEGs <- panSEA::compile_mDEG(deg, p, FDR.features, 
                                            n.dot.features)
    } else {
      compiled.DEGs <- NA
    }

  return(list(compiled.results = compiled.DEGs, all.results = deg
  ))
  }
}


