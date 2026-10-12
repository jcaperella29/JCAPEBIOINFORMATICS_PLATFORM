args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) {
    stop("Missing argument: ", flag)
  }
  args[[i + 1L]]
}

classification_dir <- get_arg("--classification-dir")
targets_string     <- get_arg("--targets")
available_file     <- get_arg("--available")
skipped_file       <- get_arg("--skipped")
selection_file     <- get_arg("--selection")

requested <- trimws(strsplit(targets_string, ",", fixed = TRUE)[[1]])
requested <- requested[nzchar(requested)]

rds_files <- list.files(
  classification_dir,
  pattern = "mixscape_classified\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

if (!length(rds_files)) {
  stop(
    "Could not find mixscape_classified.rds under: ",
    classification_dir
  )
}

rds_path <- rds_files[[1]]
message("Reading: ", rds_path)

so <- readRDS(rds_path)

md <- so@meta.data

class_col <- if ("mixscape_class_active" %in% colnames(md)) {
  "mixscape_class_active"
} else if ("mixscape_class" %in% colnames(md)) {
  "mixscape_class"
} else {
  stop(
    "No Mixscape class column found. Columns: ",
    paste(colnames(md), collapse = ", ")
  )
}

classes <- as.character(md[[class_col]])
classes <- classes[!is.na(classes)]

class_counts <- table(classes)

available_ko <- names(class_counts)[
  grepl("\\s+KO$", names(class_counts), ignore.case = TRUE)
]

available_targets <- toupper(
  sub("\\s+KO$", "", available_ko, ignore.case = TRUE)
)

requested_upper <- toupper(requested)

is_available <- requested_upper %in% available_targets

selected <- requested[is_available]
skipped  <- requested[!is_available]

writeLines(selected, available_file)

get_count <- function(label) {
  hit <- names(class_counts)[toupper(names(class_counts)) == toupper(label)]

  if (!length(hit)) {
    return(0L)
  }

  as.integer(class_counts[[hit[[1]]]])
}

selection <- data.frame(
  target = requested,
  status = ifelse(is_available, "available", "skipped"),
  ko_cells = vapply(
    requested,
    function(x) get_count(paste0(x, " KO")),
    integer(1)
  ),
  np_cells = vapply(
    requested,
    function(x) get_count(paste0(x, " NP")),
    integer(1)
  ),
  reason = ifelse(
    is_available,
    "Mixscape KO class present",
    "No Mixscape KO class present"
  ),
  stringsAsFactors = FALSE
)

write.table(
  selection,
  selection_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  selection[selection$status == "skipped", , drop = FALSE],
  skipped_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("\nAvailable targets:\n")
print(selected)

cat("\nSkipped targets:\n")
print(skipped)

cat("\nSelection table:\n")
print(selection)
