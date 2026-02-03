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
* `getsfm()` (Sample Feature Matrix) and `get_qc_features()` now set informative row names corresponding to the sample injection index.
* Expanded vignette with "Using the Results" and "Downstream Analysis" sections, demonstrating integration with `SummarizedExperiment`.
* Detailed documentation for return values and column definitions for `getsfm`, `getsff`, and `get_qc_features`.
* Updated README installation instructions for Bioconductor.
* Added a comprehensive `testthat` test suite.
* Updated `DESCRIPTION` to set `LazyData: false` per Bioconductor guidelines.
* Added `SummarizedExperiment` to `Suggests`.
* Updated `find_2d_peaks` to return an object with class c("sfi_peaks", "data.frame"). All the related functions have been updated.

