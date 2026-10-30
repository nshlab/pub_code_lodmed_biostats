if (!requireNamespace("here", quietly = TRUE)) stop("here is required")
actual_path <- here::here("results", "qc", "figure2_measurement_counts.csv")
if (!file.exists(actual_path)) stop("Run Rscript analysis/run.R figure2 first")
a <- read.csv(actual_path, check.names = FALSE)
b <- read.csv(here::here("tests", "application", "expected_figure2_measurement_counts.csv"), check.names = FALSE)
key <- function(x) paste(x$cohort, x$day, x$status, sep = "|")
a <- a[order(key(a)), names(b)]; b <- b[order(key(b)), ]
rownames(a) <- NULL; rownames(b) <- NULL
stopifnot(isTRUE(all.equal(a, b, check.attributes = FALSE, tolerance = 1e-10)))
tab <- read.csv(here::here("results", "tables", "tableS1.csv"), check.names = FALSE)
stopifnot(tab$Overall[tab$Characteristics == "N"] == "749",
          tab$Placebo[tab$Characteristics == "N"] == "371",
          tab[["A/R"]][tab$Characteristics == "N"] == "378")
message("Descriptive aggregate output checks passed; this does not verify full manuscript agreement.")
