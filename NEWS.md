# sfi 0.99.0

* Initial BioC submission.

# sfi 0.99.1

* Fix warning for C++ code
* Remove sfi.Rproj from git
* Change demo data with smaller size

# sfi 0.99.2

* Recover demo data with new algorithm for `getsff`
* Add find_peaks_low_res() function for low-resolution (unit mass) data analysis
* Add run_app() function to launch Shiny app for SFI data analysis
* Add interactive Shiny app with GUI for feature extraction and alignment
* Update NAMESPACE to export new functions and import required dependencies
* Add documentation for new functions and updated existing functions
* Add shiny and DT packages to Imports for interactive functionality
  
# sfi 0.99.3

* Renamed `getqc()` to `get_sfi_params()` to better reflect its purpose. It now returns a named vector `c(window, idelta)` instead of an unnamed vector.
* Renamed `getqcdf()` to `get_qc_features()` for clarity.
* `getmzml()` now returns a `data.frame` with named columns (`mz`, `intensity`, `rt`) instead of a matrix.
* `get_qc_features()` now sets informative row names corresponding to the sample injection index.
* Detailed documentation for return values and column definitions for `getsff`, and `get_qc_features`.
* Updated README installation instructions for Bioconductor.
* Added a comprehensive `testthat` test suite.
* Updated `DESCRIPTION` to set `LazyData: false` per Bioconductor guidelines.
* Updated `find_2d_peaks` to return an object with class c("sfi_peaks", "data.frame"). All the related functions have been updated.

# sfi 0.99.4

* `getsfm()` now returns a `SummarizedExperiment` object instead of a `data.frame`, with intensities in the `assay` slot and feature metadata in `rowData`.
* Expanded vignette with "Using the Results" and "Downstream Analysis" sections, demonstrating the `SummarizedExperiment` output.
* Added `SummarizedExperiment`, `S4Vectors`, and `methods` to `Imports`.
* Updated Shiny app and documentation to handle `SummarizedExperiment` from `getsfm()`.

