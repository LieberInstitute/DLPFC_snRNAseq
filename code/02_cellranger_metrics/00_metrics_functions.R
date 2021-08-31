find_metrics_csv <- function(x) {
    d <- list.dirs(x, recursive = FALSE)
    m_path <- file.path(d, "outs", "metrics_summary.csv")
    m_exist <- file.exists(m_path)
    res <- m_path[m_exist]
    names(res) <- basename(d[m_exist])
    return(res)
}

metrics_to_numbers <- function(x) {
    res <- sapply(x, function(col) {
        if (any(grepl("%", col))) {
            as.numeric(gsub("%", "", col))
        } else {
            as.numeric(gsub(",", "", col))
        }
    })
    as.data.frame(res)
}
