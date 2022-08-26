
## Adapted from https://github.com/LieberInstitute/spatialDLPFC/blob/main/code/analysis/12_spatial_registration_sc/01_spatial_registration_sc.R

computeEnrichment <- function(spe, var_oi, covars) {
    mat_formula <- eval(str2expression(paste("~", "0", "+", var_oi, "+", paste(covars, collapse = " + "))))

    ## Pseudo-bulk for our current BayesSpace cluster results
    message("Make psuedobulk object - ", Sys.time())
    ## I think this needs counts assay
    spe_pseudo <- scuttle::aggregateAcrossCells(
        spe,
        DataFrame(
            BayesSpace = colData(spe)[[var_oi]],
            sample_id = spe$sample_id
        )
    )

    ###############################
    message("Filter lowly expressed genes - ", Sys.time())
    rowData(spe_pseudo)$filter_expr <- edgeR::filterByExpr(spe_pseudo)
    summary(rowData(spe_pseudo)$filter_expr)
    spe_pseudo <- spe_pseudo[which(rowData(spe_pseudo)$filter_expr), ]

    message("Normalize expression - ", Sys.time())
    x <- edgeR::cpm(edgeR::calcNormFactors(spe_pseudo), log = TRUE, prior.count = 1)
    ## Verify that the gene order hasn't changed
    stopifnot(identical(rownames(x), rownames(spe_pseudo)))
    ## Fix the column names. DGEList will have samples names as Sample1 Sample2 etc
    dimnames(x) <- dimnames(spe_pseudo)
    ## Store the log normalized counts on the SingleCellExperiment object
    logcounts(spe_pseudo) <- x
    ## We don't need this 'x' object anymore
    rm(x)

    ##### get mean expression  ####
    mat_filter <- assays(spe_pseudo)$logcounts # make matrix of filtered just the log normalized counts

    #####################
    ## Build a group model

    ### access different elements of formula and check to see if they're in colData(spe_pseudo)
    terms <- attributes(terms(mat_formula))$term.labels
    terms <- terms[!grepl(":", terms)]
    for (i in seq_along(terms)) {
        if (!terms[i] %in% colnames(colData(spe_pseudo))) {
            stop("Error: formula term ", terms[i], " is not contained in colData()")
        }
    }

    # create matrix where the rownames are the sample:clusters and the columns are the other variables (spatial.cluster + region + age + sex)
    message("Create model matrix - ", Sys.time())
    mod <- model.matrix(mat_formula,
        data = colData(spe_pseudo)
    ) # binarizes factors

    ## get duplicate correlation #http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/dupcor.html
    message("Run duplicateCorrelation() - ", Sys.time())
    corfit <- limma::duplicateCorrelation(mat_filter, mod,
        block = spe_pseudo$sample_id
    )

    ## Next for each layer test that layer vs the rest
    cluster_idx <- rafalib::splitit(colData(spe_pseudo)[, var_oi])

    message("Run enrichment statistics - ", Sys.time())
    eb0_list_cluster <- lapply(cluster_idx, function(x) {
        res <- rep(0, ncol(spe_pseudo))
        res[x] <- 1
        res_formula <- paste("~", "res", "+", paste(covars, collapse = " + "))
        m <- with(
            colData(spe_pseudo),
            model.matrix(eval(str2expression(res_formula)))
        )

        # josh suggested use top table as a wrapper because it makes the output of eBayes nicer

        limma::eBayes(
            limma::lmFit(
                mat_filter,
                design = m,
                block = spe_pseudo$sample_id,
                correlation = corfit$consensus.correlation
            )
        )
    })


    message("extract and reformat enrichment results - ", Sys.time())
    ##########
    ## Extract the p-values
    pvals0_contrasts_cluster <- sapply(eb0_list_cluster, function(x) {
        x$p.value[, 2, drop = FALSE]
    })
    rownames(pvals0_contrasts_cluster) <- rownames(mat_filter)

    t0_contrasts_cluster <- sapply(eb0_list_cluster, function(x) {
        x$t[, 2, drop = FALSE]
    })
    rownames(t0_contrasts_cluster) <- rownames(mat_filter)

    fdrs0_contrasts_cluster <- apply(pvals0_contrasts_cluster, 2, p.adjust, "fdr")
    rownames(fdrs0_contrasts_cluster) <- rownames(mat_filter)

    data.frame(
        "FDRsig" = colSums(fdrs0_contrasts_cluster < 0.05 &
            t0_contrasts_cluster > 0),
        "Pval10-6sig" = colSums(pvals0_contrasts_cluster < 1e-6 &
            t0_contrasts_cluster > 0),
        "Pval10-8sig" = colSums(pvals0_contrasts_cluster < 1e-8 &
            t0_contrasts_cluster > 0)
    )

    results_specificity <-
        f_merge(spe_pseudo, p = pvals0_contrasts_cluster, fdr = fdrs0_contrasts_cluster, t = t0_contrasts_cluster)
    head(results_specificity)

    # object to return
    results_specificity <- as.data.frame(results_specificity@listData)
    rownames(results_specificity) <- results_specificity$ensembl
    message("return enrichment results")
    return(results_specificity)
}


f_merge <- function(spe_pseudo, p, fdr, t) {
    colnames(p) <- paste0("p_value_", colnames(p))
    colnames(fdr) <- paste0("fdr_", colnames(fdr))
    colnames(t) <- paste0("t_stat_", colnames(t))
    res <- as.data.frame(cbind(t, p, fdr))
    res$ensembl <- rownames(res)
    ## Check it's all in order
    res <- merge(
        x = res, y = rowData(spe_pseudo)[, c("gene_id", "gene_name")], by.x = "ensembl",
        by.y = "gene_id"
    )
    colnames(res)[ncol(res)] <- "gene"
    rownames(res) <- res$ensembl
    return(res)
}
