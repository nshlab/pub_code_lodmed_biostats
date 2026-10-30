# Historical machine label "<LLoQ" means at OR below 2.0 in this implementation.
classify_rna_numeric <- function(x, lloq = 2.0) {
  ifelse(is.na(x), NA_character_, ifelse(x > lloq, "Observed", "<LLoQ"))
}
rna_status_at <- function(d, day, lloq = 2.0) {
  value <- d[[paste0("day", day, "_value")]]
  if (is.null(value)) stop("Missing RNA value column for Day ", day)
  status <- classify_rna_numeric(value, lloq)
  evidence <- rep(FALSE, nrow(d)); below <- rep(FALSE, nrow(d))
  for (suffix in c("_text", "_cat1")) {
    field <- paste0("day", day, suffix)
    if (!field %in% names(d)) next
    text <- toupper(trimws(as.character(d[[field]])))
    present <- !is.na(text) & nzchar(text)
    present <- present & !text %in% c("NA", "MISSING", "NOT DONE", "NOT COLLECTED")
    evidence <- evidence | present
    below <- below | (present & grepl("<[ =]*LLOQ|BELOW.*(QUANT|LLOQ)|NOT DETECTED|UNDETECT", text))
  }
  absent <- is.na(value)
  status[absent & below] <- "<LLoQ"
  status[absent & !below & evidence] <- "Unclassified"
  status[absent & !evidence] <- "Missing"
  status
}
