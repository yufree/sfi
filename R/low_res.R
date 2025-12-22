#' Find peaks in low-resolution data using the 2D peak finding algorithm
#'
#' This function adapts the fast `find_2d_peaks` function for use with
#' low-resolution (unit mass) data. It does this by first aggregating the
#' signal at each integer mass and then calling the 2D peak finder.
#'
#' @param mz A numeric vector of mass-to-charge ratios.
#' @param rt A numeric vector of retention times.
#' @param intensity A numeric vector of intensities.
#' @param deltart Numeric. Tolerance for retention time matching. Default is 5.
#' @param snr Numeric. Signal-to-noise ratio to find peaks. Default is 3.0.
#'
#' @return A data frame with columns `mz`, `rt`, and `intensity`,
#'         representing the detected peaks. The `mz` values are integer masses.
#' @examples
#' data(sfi)
#' peaks <- find_peaks_low_res(
#'     mz = sfi$mz, rt = sfi$rt,
#'     intensity = sfi$intensity
#' )
#' @export
find_peaks_low_res <- function(mz,
                                rt,
                                intensity,
                                deltart = 5,
                                snr = 3.0) {

  # Round mz to the nearest integer
  mz_int <- round(mz)

  # Create a combined key for aggregation
  keys <- paste(mz_int, rt, sep = "_")

  # Aggregate intensity by summing up intensities for the same integer mass and retention time
  # This is much faster than aggregate() for large datasets
  aggregated_intensities <- tapply(intensity, keys, sum)

  # Extract unique mz and rt combinations
  key_list <- strsplit(names(aggregated_intensities), "_")
  aggregated_mz <- as.numeric(vapply(key_list, function(x) x[1],
                                     FUN.VALUE = character(1)))
  aggregated_rt <- as.numeric(vapply(key_list, function(x) x[2],
                                    FUN.VALUE = character(1)))
  aggregated_intensity <- as.numeric(aggregated_intensities)

  # Calculate low_res_mz_bins
  if (length(aggregated_mz) > 0) {
    min_mz_int <- min(aggregated_mz)
    max_mz_int <- max(aggregated_mz)
    low_res_mz_bins <- as.integer(max_mz_int - min_mz_int + 1)
    low_res_mz_bins <- max(1L, low_res_mz_bins) # Ensure at least 1 bin
  } else {
    low_res_mz_bins <- 1L # Fallback for empty aggregated_mz
  }
  
  # Dynamically calculate ppm to keep unit masses separate
  if (length(aggregated_mz) > 0) {
      max_mz <- max(aggregated_mz)
      if (max_mz > 0) {
          # Set ppm so that the mz_window at max_mz is < 0.5
          # max_mz * ppm * 1e-6 < 0.5  => ppm < 0.5 / (max_mz * 1e-6)
          # We use 0.4 for a safety margin.
          dynamic_ppm <- 0.4 / (max_mz * 1e-6)
      } else {
          dynamic_ppm <- 1 # Fallback for edge case where max_mz is 0
      }
  } else {
      dynamic_ppm <- 1 # Fallback for empty data
  }

  # Call the fast 2D peak finder on the aggregated data
  peaks <- find_2d_peaks(
    mz = aggregated_mz,
    rt = aggregated_rt,
    intensity = aggregated_intensity,
    ppm = dynamic_ppm,  # Use dynamically calculated ppm
    deltart = deltart,
    snr = snr,
    mz_bins = low_res_mz_bins,  # Use calculated mz_bins
    rt_bins = NULL  # Let find_2d_peaks calculate automatically
  )

  return(peaks)
}
