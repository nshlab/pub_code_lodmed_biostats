if (!requireNamespace("here", quietly = TRUE)) stop("Package 'here' is required.")

tests <- sort(list.files(
  here::here("tests", "simulation"),
  pattern = "^test_.*\\.R$",
  full.names = TRUE
))

status <- data.frame(
  test = basename(tests),
  passed = FALSE,
  message = NA_character_,
  stringsAsFactors = FALSE
)

for (i in seq_along(tests)) {
  cat("\n--- ", basename(tests[[i]]), " ---\n", sep = "")
  ans <- tryCatch(
    {
      sys.source(tests[[i]], envir = new.env(parent = globalenv()))
      list(ok = TRUE, msg = "PASS")
    },
    error = function(e) list(ok = FALSE, msg = conditionMessage(e))
  )
  status$passed[[i]] <- ans$ok
  status$message[[i]] <- ans$msg
}

print(status, row.names = FALSE)

out <- here::here("tests", "simulation", "regression_test_status.csv")
utils::write.csv(status, out, row.names = FALSE)

if (!all(status$passed)) {
  stop("One or more simulation regression tests failed. See: ", out)
}
message("All simulation regression tests passed.")
