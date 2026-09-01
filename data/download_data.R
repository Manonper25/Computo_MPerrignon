# Downloads large industrial-case datasets from Zenodo
# DOI: 10.5281/zenodo.22086876

dir.create("data/industrial_case", showWarnings = FALSE, recursive = TRUE)

zenodo_files <- list(
  "results_industrial_case.rds" = "https://zenodo.org/records/22086876/files/results_industrial_case.rds",
  "bornes_max_optim.rds" = "https://zenodo.org/records/22086876/files/bornes_max_optim.rds",
  "bornes_min_optim.rds" = "https://zenodo.org/records/22086876/files/bornes_min_optim.rds",
  "modY1.rds" = "https://zenodo.org/records/22086876/files/modY1.rds",
  "modY2.rds" = "https://zenodo.org/records/22086876/files/modY2.rds",
  "modY3.rds" = "https://zenodo.org/records/22086876/files/modY3.rds",
  "modY4.rds" = "https://zenodo.org/records/22086876/files/modY4.rds",
  "results_HV.rds" = "https://zenodo.org/records/22086876/files/results_HV.rds",
  "Var_fixe_optim.rds" = "https://zenodo.org/records/22086876/files/Var_fixe_optim.rds"
)

for (fname in names(zenodo_files)) {
  dest <- file.path("data/industrial_case", fname)
  if (!file.exists(dest)) {
    message("Downloading ", fname, "...")
    download.file(zenodo_files[[fname]], destfile = dest, mode = "wb")
  } else {
    message(fname, " already exists, skipping.")
  }
}