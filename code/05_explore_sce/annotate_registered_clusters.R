annotate_registered_clusters <- function(cor_matrix, confidence_threshold = 0.25, cutoff_merge_ratio = 0.25) {
    annotated <- apply(cor_matrix, 1, annotate_registered_cluster, cutoff_merge_ratio = cutoff_merge_ratio)

    annotated <- gsub("ayer", "", annotated)
    annotated <- gsub("\\/L", "\\/", annotated)
    annotated <- gsub("^WM\\/", "WM\\/L", annotated)

    data.frame(
        cluster = names(annotated),
        layer_confidence = ifelse(apply(cor_matrix, 1, max) > confidence_threshold, "good", "poor"),
        layer_label = annotated,
        row.names = NULL
    )
}

annotate_registered_cluster <- function(remaining, label = "", current = NULL, cutoff_merge_ratio = 0.25) {
    ## Filter negative correlations
    remaining <- remaining[remaining > 0]

    ## There's nothing else to continue with
    if(length(remaining) == 0) {
        return(label)
    }

    ## Find the next highest correlation
    next_i <- which.max(remaining)
    next_cor <- remaining[next_i]

    if(label == "") {
        ## Initial case when we didn't have a label
        annotate_registered_cluster(
            remaining = remaining[-next_i],
            label = names(next_cor),
            current = next_cor,
            cutoff_merge_ratio = cutoff_merge_ratio
        )
    } else {
        ## Find the difference, then divide by the next correlation
        next_diff_ratio <- (current - next_cor) / next_cor

        if (next_diff_ratio > cutoff_merge_ratio) {
            ## It's above the cutoff, so we don't decide to merge
            ## and are done =)
            return(label)
        } else {
            ## It's below the cutoff, so we need to look at the next one
            annotate_registered_cluster(
                remaining = remaining[-next_i],
                label = paste0(label, "/", names(next_cor)),
                current = next_cor,
                cutoff_merge_ratio = cutoff_merge_ratio
            )
        }
    }
}
