#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
opts <- list()

idx <- 1
while (idx <= length(args)) {
  key <- args[[idx]]
  if (substr(key, 1, 2) == "--") {
    name <- substring(key, 3)
    if (idx == length(args)) {
      stop(paste("Missing value for", key))
    }
    opts[[name]] <- args[[idx + 1]]
    idx <- idx + 2
    next
  }
  idx <- idx + 1
}

required <- c("rmd", "output")
missing <- required[sapply(required, function(name) {
  is.null(opts[[name]]) || opts[[name]] == ""
})]
if (length(missing) > 0) {
  stop(paste("Missing required arguments:", paste(missing, collapse = ", ")))
}

rmd_path <- opts[["rmd"]]
output_path <- opts[["output"]]

if (!file.exists(rmd_path)) {
  stop(paste("Rmd template not found:", rmd_path))
}

out_dir <- dirname(output_path)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

render_params <- list(
  dataframe_csv = opts[["dataframe_csv"]],
  report_html = opts[["report_html"]],
  data_dir = opts[["data_dir"]],
  output_dir = opts[["output_dir"]],
  pegas_version = opts[["pegas_version"]],
  pegas_install_dir = opts[["pegas_install_dir"]]
)

log_line <- function(text) {
  cat(paste0("[pegas] ", text, "\n"))
}

base_dirs <- c(getwd())
if (!is.null(render_params$output_dir) && render_params$output_dir != "") {
  base_dirs <- c(base_dirs, render_params$output_dir)
}
if (!is.null(render_params$report_html) && render_params$report_html != "") {
  base_dirs <- c(base_dirs, dirname(dirname(render_params$report_html)))
}
if (!is.null(render_params$dataframe_csv) && render_params$dataframe_csv != "") {
  base_dirs <- c(base_dirs, dirname(dirname(render_params$dataframe_csv)))
}
if (!is.null(render_params$pegas_install_dir) && render_params$pegas_install_dir != "") {
  pegas_parent <- dirname(render_params$pegas_install_dir)
  base_dirs <- c(base_dirs, pegas_parent, dirname(pegas_parent))
}
base_dirs <- unique(base_dirs)
log_line(paste("FastQC base dirs:", paste(base_dirs, collapse = " | ")))

fastqc_dirs <- c("fastqc")
for (base in base_dirs) {
  fastqc_dirs <- c(
    fastqc_dirs,
    file.path(base, "fastqc"),
    file.path(base, "out", "fastqc")
  )
}
fastqc_dirs <- unique(fastqc_dirs[dir.exists(fastqc_dirs)])
if (length(fastqc_dirs) == 0) {
  log_line("FastQC search dirs (existing): <none>")
} else {
  log_line(paste("FastQC search dirs (existing):", paste(fastqc_dirs, collapse = " | ")))
}

html_files <- unlist(lapply(fastqc_dirs, function(d) {
  list.files(d, pattern = "\\.html$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
}))
if (length(html_files) == 0) {
  fallback_files <- unlist(lapply(base_dirs, function(d) {
    if (!dir.exists(d)) {
      return(character(0))
    }
    list.files(
      d,
      pattern = "_fastqc\\.html$",
      full.names = TRUE,
      recursive = TRUE,
      ignore.case = TRUE
    )
  }))
  html_files <- unique(fallback_files)
}
html_files <- normalizePath(html_files, winslash = "/", mustWork = FALSE)
html_files <- unique(html_files[file.exists(html_files)])
log_line(paste("FastQC HTML files found:", length(html_files)))
if (length(html_files) > 0) {
  preview <- paste(utils::head(html_files, 3), collapse = " | ")
  log_line(paste("FastQC sample files:", preview))
}

render_params$fastqc_files <- html_files

suppressPackageStartupMessages(library(rmarkdown))

render_env <- new.env(parent = globalenv())

rmarkdown::render(
  input = rmd_path,
  output_file = basename(output_path),
  output_dir = out_dir,
  params = render_params,
  envir = render_env,
  quiet = TRUE
)
