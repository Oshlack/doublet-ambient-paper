get_package_citation <- function(pkg_name, output_dir = NULL) {
  # if output dir is NULL, then make a folder in working dir called package_refs
  if (is.null(output_dir)) {
    output_dir <- paste0(getwd(), "/package_refs/"); dir.create(output_dir, showWarnings = F)
  }

  output_file <- paste0(output_dir, pkg_name, ".bib")

  # Check the package is installed
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    stop(paste("Package", pkg_name, "is not installed."))
  }

  # Get citation object
  cit <- citation(pkg_name)

  # Write to bibtex file
  toBibtex(cit) |> writeLines(output_file)

  message(paste("Citation written to", output_file))
  invisible(output_file)
}


#example
get_package_citation("Seurat", output_dir = "output/package_refs/")


packages <- c("Biobase", "scuttle", "scater", "DropletUtils", "yardstick", "coop", "renv", "SingleCellExperiment",
              "ggplot2", "patchwork", "ggpubr")

lapply(packages, function(pkg) get_package_citation(pkg, output_dir = "output/package_refs/"))

get_package_citation("caret", output_dir = "output/package_refs/")




