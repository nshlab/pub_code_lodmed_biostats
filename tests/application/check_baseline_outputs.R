if (!requireNamespace("here", quietly = TRUE)) stop("here is required")
actual <- read.csv(here::here("results", "qc", "tableS1_baseline_rna.csv"), check.names = FALSE)
expected <- read.csv(here::here("tests", "application", "expected_baseline_rna.csv"), check.names = FALSE)
actual <- actual[match(expected$group, actual$group), names(expected)]
rownames(actual) <- NULL; rownames(expected) <- NULL
stopifnot(isTRUE(all.equal(actual, expected, check.attributes = FALSE, tolerance = 1e-10)))
tab <- read.csv(here::here("results", "tables", "tableS1.csv"), check.names = FALSE)
ix <- which(tab$Characteristics == "AN RNA day 0 (log10 copies/mL)")
stopifnot(length(ix) == 1L, tab$Overall[ix] == "3.84 (2.41)",
          tab$Placebo[ix] == "3.85 (2.39)", tab[["A/R"]][ix] == "3.83 (2.43)")
message("Baseline Table S1 row and group-specific availability QC match the independent benchmark.")
