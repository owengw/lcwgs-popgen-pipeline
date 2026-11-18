ind_dir  <- "/mnt/parscratch/users/bi4og/physalia/angsd_results/roh/individual"
pop_dir  <- "/mnt/parscratch/users/bi4og/physalia/angsd_results/roh/population"
meta_file <- "/mnt/parscratch/users/bi4og/physalia/angsd_results/metadata.tsv"

dir.create(pop_dir, showWarnings = FALSE, recursive = TRUE)

metadata <- read.table(meta_file, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)
metadata$sample_prefix <- sub("_.*", "", metadata$sample_id)

summary_files <- list.files(ind_dir, pattern = "[.]summary[.]txt$",
                             full.names = TRUE)
cat("Summary files found:", length(summary_files), "\n")

parse_rohan_summary <- function(f) {
    sample <- sub("[.]summary[.]txt$", "", basename(f))
    lines  <- tryCatch(readLines(f), error = function(e) NULL)
    if (is.null(lines) || length(lines) == 0) return(NULL)

    get_val <- function(pattern) {
        m <- grep(pattern, lines, value = TRUE, ignore.case = TRUE)
        if (length(m) == 0) return(NA_real_)
        tok <- sub(".*:[[:space:]]*([0-9eE.+-]+).*", "\\1", m[1])
        suppressWarnings(as.numeric(tok))
    }

    data.frame(
        sample_id    = sample,
        theta        = get_val("Genome-wide theta outside ROH"),
        n_roh        = get_val("Segments in ROH[[:space:]]*:"),
        pct_roh      = get_val("Segments in ROH[(]"),
        avg_roh_bp   = get_val("Avg[.] length of ROH"),
        stringsAsFactors = FALSE
    )
}

ind_list <- lapply(summary_files, parse_rohan_summary)
ind_list <- ind_list[!sapply(ind_list, is.null)]
ind_df   <- do.call(rbind, ind_list)

cat("\nExample parsed values (first 5):\n")
print(head(ind_df, 5))

# Match to population
ind_df$sample_prefix <- sub("_.*", "", ind_df$sample_id)
ind_df <- merge(ind_df,
                metadata[, c("sample_prefix", "site")],
                by = "sample_prefix", all.x = TRUE)

ind_df$avg_roh_kb <- ind_df$avg_roh_bp / 1000

pop_order <- c("CAI", "NS", "SB", "SI")

pop_summary <- do.call(rbind, lapply(pop_order, function(pop) {
    d <- ind_df[!is.na(ind_df$site) & ind_df$site == pop, ]
    if (nrow(d) == 0) return(NULL)
    data.frame(
        population      = pop,
        n_individuals   = nrow(d),
        mean_theta      = round(mean(d$theta,      na.rm = TRUE), 6),
        mean_n_roh      = round(mean(d$n_roh,      na.rm = TRUE), 2),
        mean_pct_roh    = round(mean(d$pct_roh,    na.rm = TRUE), 4),
        mean_avg_roh_kb = round(mean(d$avg_roh_kb, na.rm = TRUE), 3),
        stringsAsFactors = FALSE
    )
}))

cat("\n=== Population ROH summary ===\n")
print(pop_summary)

write.table(ind_df, file.path(pop_dir, "individual_ROH_all.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(pop_summary, file.path(pop_dir, "population_ROH_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nDone. Results in:", pop_dir, "\n")
