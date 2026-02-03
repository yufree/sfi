#' @useDynLib sfi, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL
#' Demo sfi data
#' @docType data
#' @usage data(sfi)
#' @format A data.frame object with mass to charge ratio, intensity and retention time from sfi mode.
"sfi"
#' Read mzML File and Extract m/z, Retention Time, and Intensity
#' @param path path of SFI mzML file.
#' @return A data frame containing m/z, retention time and intensity.
#' @examples
#' # Load demo data
#' data(sfi)
#' head(sfi)
#' # In practice, you would use a real mzML file path:
#' # peak <- getmzml("path/to/your/file.mzML")
#' # The function returns a data frame with m/z, retention time, and intensity columns
#' @export
getmzml <- function(path) {
  mzml_file <- mzR::openMSfile(path)
  # read meta data
  tt <- mzR::header(mzml_file)
  # generate retention time vector
  rt <- rep(tt$retentionTime, tt$peaksCount)
  # extract peaks
  peaks <- mzR::peaks(mzml_file)
  peak <- do.call(rbind, peaks)
  peak <- cbind(peak, rt)
  peak <- as.data.frame(peak)
  colnames(peak) <- c("mz", "intensity", "rt")
  return(peak)
}
#' Cluster and Pair m/z and Retention Time Features
#'
#' This function clusters m/z values based on Manhattan distance and pairs features within clusters.
#'
#' @param mz Numeric vector of m/z values.
#' @param rt Numeric vector of retention times corresponding to m/z values.
#' @param ppm Numeric. Parts per million tolerance for m/z matching. Default is 5.
#' @param minn Integer. Minimum number of features in a cluster to be retained. Default is 2.
#' @param refmz Optional numeric vector of reference m/z values for alignment. Default is NULL.
#'
#' @return A data frame containing paired m/z and retention time values with their differences:
#' \itemize{
#'   \item mz1: m/z of the first feature in the pair.
#'   \item rt1: retention time of the first feature in the pair.
#'   \item mz2: m/z of the second feature in the pair.
#'   \item rt2: retention time of the second feature in the pair.
#'   \item pmr: absolute difference in retention time (Pair Mass Retention).
#'   \item pmd: absolute difference in m/z (Pair Mass Difference).
#' }
#' @examples
#' data(sfi)
#' peak <- find_2d_peaks(mz=sfi$mz,rt=sfi$rt,intensity=sfi$intensity)
#' sff_features <- getsff(peak$mz, peak$rt)
#' @param mz Numeric vector of m/z values or an object of class `sfi_peaks`.
#' @param rt Numeric vector of retention times corresponding to m/z values.
#' @param ... Additional arguments passed to methods.
#' @export
getsff <- function(mz, rt, ...) {
  UseMethod("getsff")
}

#' @describeIn getsff Method for sfi_peaks object
#' @export
getsff.sfi_peaks <- function(mz, rt = NULL, ...) {
  getsff.default(mz$mz, mz$rt, ...)
}

#' @describeIn getsff Default method for vectors
#' @export
getsff.default <- function(mz,
                   rt,
                   ppm = 5,
                   minn = 2,
                   refmz = NULL, ...) {
  if (!is.unsorted(rt)) {
    warning("The input retention time vector appears to be sorted. This suggests raw data is being used instead of peak-picked data. Please run find_2d_peaks() first.")
  }
  # Fast replacement for hclust
  if (length(mz) < 2) {
    return(data.frame())
  }
  ord <- order(mz)
  mz_sorted <- mz[ord]
  
  group_ids_sorted <- integer(length(mz))
  group_id <- 1
  group_ids_sorted[1] <- group_id
  
  if (length(mz) > 1) {
    for (i in 2:length(mz)) {
      if ((mz_sorted[i] - mz_sorted[i-1]) * 1e6 / mz_sorted[i-1] > ppm) {
        group_id <- group_id + 1
      }
      group_ids_sorted[i] <- group_id
    }
  }
  
  # Create mzcluster vector in original order
  mzcluster <- integer(length(mz))
  mzcluster[ord] <- group_ids_sorted

  # Identify clusters with at least 'minn' members
  valid_clusters <- names(which(table(mzcluster) >= minn))
  idx <- mzcluster %in% as.numeric(valid_clusters)

  # Filter m/z, rt, and cluster assignments
  mz <- mz[idx]
  rt <- rt[idx]
  mzcluster <- mzcluster[idx]

  # If reference m/z values are provided, align to reference
  if (!is.null(refmz)) {
    merge <- enviGCMS::getalign(mz, refmz, ppm = ppm)
    idx <- unique(merge$xid)
    mz <- mz[idx]
    rt <- rt[idx]
    mzcluster <- mzcluster[idx]
  }

  # Helper function to process each cluster
  getmzbin <- function(cluster_id) {
    # Extract m/z and rt for the current cluster
    bin <- data.frame(mz = mz[mzcluster == cluster_id], rt = rt[mzcluster == cluster_id])

    # Use combn to get all pairs of indices
    if (nrow(bin) < 2) return(NULL)
    pairs <- utils::combn(nrow(bin), 2)

    # Pair retention times and m/z values
    rt1 <- bin$rt[pairs[1,]]
    rt2 <- bin$rt[pairs[2,]]
    mz1 <- bin$mz[pairs[1,]]
    mz2 <- bin$mz[pairs[2,]]

    # Calculate differences
    pmr <- abs(rt1 - rt2)
    pmd <- abs(mz1 - mz2)

    # Apply ppm tolerance filter
    idx <- pmd / rowMeans(cbind(mz1, mz2)) < (ppm * 1e-6)

    # Return filtered paired data
    data.frame(
      mz1 = mz1[idx],
      rt1 = rt1[idx],
      mz2 = mz2[idx],
      rt2 = rt2[idx],
      pmr = pmr[idx],
      pmd = pmd[idx]
    )
  }

  # Apply the helper function to all unique clusters
  paired_list <- lapply(unique(mzcluster), getmzbin)

  # Combine all paired data into a single data frame
  paired_df <- do.call(rbind, paired_list)

  return(paired_df)
}

#' Determine Optimal Retention Time Window
#'
#' This function calculates the optimal retention time window based on QC sequences and m/z/rt data.
#'
#' @param mz Numeric vector of m/z values.
#' @param rt Numeric vector of retention times.
#' @param lower Numeric. Lower bound for the retention time window. Default is 620.
#' @param upper Numeric. Upper bound for the retention time window. Default is 650.
#' @param ppm Numeric. Parts per million tolerance for m/z matching. Default is 5.
#' @param minn Integer. Minimum number of features in a QC cluster. Default is 1.
#' @param qcseq Integer vector. QC sequence indicating which samples are QC. Default is c(1, 1, 0, 1, 1, 0, 1, 1, 0).
#'
#' @return Numeric value representing the optimal retention time window.
#' @examples
#' data(sfi)
#' peak <- find_2d_peaks(mz=sfi$mz,rt=sfi$rt,intensity=sfi$intensity)
#' window_opt <- getwindow(peak$mz, peak$rt)
#' @param mz Numeric vector of m/z values or an object of class `sfi_peaks`.
#' @param rt Numeric vector of retention times.
#' @param ... Additional arguments passed to methods.
#' @export
getwindow <- function(mz, rt, ...) {
  UseMethod("getwindow")
}

#' @describeIn getwindow Method for sfi_peaks object
#' @export
getwindow.sfi_peaks <- function(mz, rt = NULL, ...) {
  getwindow.default(mz$mz, mz$rt, ...)
}

#' @describeIn getwindow Default method for vectors
#' @export
getwindow.default <- function(mz,
                      rt,
                      lower = 620,
                      upper = 650,
                      ppm = 5,
                      minn = 1,
                      qcseq = c(1, 1, 0, 1, 1, 0, 1, 1, 0), ...) {
  if (!is.unsorted(rt)) {
    warning("The input retention time vector appears to be sorted. This suggests raw data is being used instead of peak-picked data. Please run find_2d_peaks() first.")
  }
  # Calculate the cutoff retention time based on QC sequence length and upper window
  recut <- length(qcseq) * upper

  # Filter m/z and rt values within the cutoff
  valid_idx <- rt < recut
  mz_filtered <- mz[valid_idx]
  rt_filtered <- rt[valid_idx]

  # Generate paired m/z and rt data
  sff <- getsff(mz_filtered, rt_filtered, ppm = ppm, minn = minn)

  # Determine the most frequent rounded pmr values corresponding to QC samples
  rtc <- as.numeric(names(table(round(sff$pmr, 2)))[order(table(round(sff$pmr, 2)), decreasing = TRUE)])[seq_along(qcseq)]

  # Filter rtc within specified window bounds
  rtcx <- rtc[rtc < upper & rtc > lower]

  # Determine the window as the mean of rtcx or handle edge cases
  if (length(rtcx) > 1) {
    window <- mean(rtcx)
  } else if (length(rtcx) == 0) {
    stop("Matrix effects found. You might increase qcseq's length")
  } else {
    window <- rtcx
  }

  message("Window is ", window)
  return(window)
}

#' Optimize Delta Retention Time
#'
#' This function optimizes the delta retention time (idelta) using a binary search approach.
#'
#' @param mz Numeric vector of m/z values.
#' @param rt Numeric vector of retention times.
#' @param qcmz Numeric vector of QC m/z values.
#' @param qcrt Numeric vector of QC retention times.
#' @param idelta Initial delta retention time guess. Default is 60.
#' @param shift Numeric. Shift applied to idelta. Default is 0.
#' @param ppm Numeric. Parts per million tolerance for m/z matching. Default is 5.
#' @param deltart Numeric. Tolerance for retention time matching. Default is 5.
#' @param window Numeric. Retention time window. Default is 600.
#' @param n Integer. Number of iterations or samples. Default is 160.
#' @param tol Numeric. Tolerance for binary search convergence. Default is 0.03.
#' @param max_iter Integer. Maximum number of binary search iterations. Default is 100.
#'
#' @return Optimized delta retention time (idelta).
#' @examples
#' data(sfi)
#' peak <- find_2d_peaks(mz=sfi$mz,rt=sfi$rt,intensity=sfi$intensity)
#' delta_opt <- getidelta(peak$mz, peak$rt,qcmz=195.0876,qcrt=74,window=632,idelta=90)
#' @param mz Numeric vector of m/z values or an object of class `sfi_peaks`.
#' @param rt Numeric vector of retention times.
#' @param ... Additional arguments passed to methods.
#' @export
getidelta <- function(mz, rt, ...) {
  UseMethod("getidelta")
}

#' @describeIn getidelta Method for sfi_peaks object
#' @export
getidelta.sfi_peaks <- function(mz, rt = NULL, ...) {
  getidelta.default(mz$mz, mz$rt, ...)
}

#' @describeIn getidelta Default method for vectors
#' @export
getidelta.default <- function(mz,
                      rt,
                      qcmz,
                      qcrt,
                      idelta = 60,
                      shift = 0,
                      ppm = 5,
                      deltart = 5,
                      window = 600,
                      n = 160,
                      tol = 0.03,
                      max_iter = 100, ...) {
  if (!is.unsorted(rt)) {
    warning("The input retention time vector appears to be sorted. This suggests raw data is being used instead of peak-picked data. Please run find_2d_peaks() first.")
  }
  # Align QC m/z with sample m/z values
  merge2 <- enviGCMS::getalign(qcmz, mz, ppm = ppm)
  qcmz_aligned <- qcmz[unique(merge2$xid)]
  qcrt_aligned <- qcrt[unique(merge2$xid)]

  # Initialize retention time QC vector
  rtqc <- rep(NA, length(mz))

  # Assign QC retention times based on m/z alignment within ppm tolerance
  for (i in seq_along(qcmz_aligned)) {
    idx1 <- which(abs(mz - qcmz_aligned[i]) / qcmz_aligned[i] < ppm * 1e-6)
    rtqc[idx1] <- qcrt_aligned[i]
  }

  # Calculate intermediate retention time
  irt <- rt - window - rtqc

  # Define search bounds for binary search
  ideltax <- idelta + shift
  upper <- ideltax + tol * 10
  lower <- ideltax - tol * 10

  # Function to evaluate the number of valid matches for a given idelta
  evaluate_idelta <- function(current_idelta, irt, n, deltart) {
    sampleidx <- round(irt / current_idelta)
    shiftrt <- rt - window - sampleidx * current_idelta - rtqc
    sum(abs(shiftrt) < deltart &
          sampleidx %in% seq_len(n), na.rm = TRUE)
  }

  # Binary search optimization to find the best idelta
  binary_search_optimization <- function(f,
                                         lower,
                                         upper,
                                         irt,
                                         n,
                                         deltart,
                                         tol = 0.01,
                                         max_iter = 100) {
    iter <- 0
    converged <- FALSE

    while (!converged && iter < max_iter) {
      mid <- (lower + upper) / 2

      f_mid <- f(mid, irt, n, deltart)
      f_mid_tol <- f(mid + tol, irt, n, deltart)

      if (f_mid < f_mid_tol) {
        lower <- mid
      } else {
        upper <- mid
      }

      if (abs(upper - lower) < tol) {
        converged <- TRUE
      }

      iter <- iter + 1
    }

    if (converged) {
      return((lower + upper) / 2)
    } else {
      warning("Binary search did not converge within the maximum number of iterations.")
      return(NULL)
    }
  }

  # Perform binary search to optimize idelta
  optimized_idelta <- binary_search_optimization(
    f = evaluate_idelta,
    lower = lower,
    upper = upper,
    irt = irt,
    n = n,
    deltart = deltart,
    tol = tol,
    max_iter = max_iter
  )

  if (is.null(optimized_idelta)) {
    stop("Failed to optimize idelta.")
  }

  message("Delta retention time is ", optimized_idelta)
  return(optimized_idelta)
}

#' Quality Control for Mass Spectrometry Data
#'
#' This function performs quality control (QC) on mass spectrometry data by aligning QC and sample features.
#'
#' @param mz Numeric vector of m/z values.
#' @param rt Numeric vector of retention times.
#' @param intensity Numeric vector of intensities corresponding to m/z and rt values.
#' @param idelta Numeric. Initial delta retention time. Default is 60.
#' @param window Numeric. Retention time window. Default is 600.
#' @param qcseq Integer vector indicating QC samples. Default is c(1, 1, 0, 1, 1, 0, 1, 1, 0).
#' @param deltart Numeric. Tolerance for retention time matching. Default is 5.
#' @param ppm Numeric. Parts per million tolerance for m/z matching. Default is 5.
#' @param minn Integer. Minimum number of QC samples required. Default is 1.
#' @param n Integer. Number of samples for delta optimization. Default is 160.
#' @param tol Numeric. Tolerance for binary search in delta optimization. Default is 0.03.
#' @param max_iter Integer. Maximum iterations for binary search. Default is 100.
#' @param wlower Numeric. Lower bound for window determination. Default is 620.
#' @param wupper Numeric. Upper bound for window determination. Default is 650.
#'
#' @return A named numeric vector containing the optimal window and delta retention time:
#' \itemize{
#'   \item window: The optimized retention time window.
#'   \item idelta: The optimized delta retention time.
#' }
#' @examples
#' data(sfi)
#' peak <- find_2d_peaks(mz=sfi$mz,rt=sfi$rt,intensity=sfi$intensity)
#' sfi_params <- get_sfi_params(peak$mz, peak$rt, peak$intensity, deltart=10)
#' @param mz Numeric vector of m/z values or an object of class `sfi_peaks`.
#' @param rt Numeric vector of retention times.
#' @param intensity Numeric vector of intensities corresponding to m/z and rt values.
#' @param ... Additional arguments passed to methods.
#' @export
get_sfi_params <- function(mz, rt, intensity, ...) {
  UseMethod("get_sfi_params")
}

#' @describeIn get_sfi_params Method for sfi_peaks object
#' @export
get_sfi_params.sfi_peaks <- function(mz, rt = NULL, intensity = NULL, ...) {
  get_sfi_params.default(mz$mz, mz$rt, mz$intensity, ...)
}

#' @describeIn get_sfi_params Default method for vectors
#' @export
get_sfi_params.default <- function(mz,
                  rt,
                  intensity,
                  idelta = 60,
                  window = 600,
                  qcseq = c(1, 1, 0, 1, 1, 0, 1, 1, 0),
                  deltart = 5,
                  ppm = 5,
                  minn = 1,
                  n = 160,
                  tol = 0.03,
                  max_iter = 100,
                  wlower = 620,
                  wupper = 650, ...) {
  if (!is.unsorted(rt)) {
    warning("The input retention time vector appears to be sorted. This suggests raw data is being used instead of peak-picked data. Please run find_2d_peaks() first.")
  }
  # Determine the optimal retention time window
  windows <- getwindow(
    mz = mz,
    rt = rt,
    lower = wlower,
    upper = wupper,
    ppm = ppm,
    minn = minn,
    qcseq = qcseq
  )

  # Calculate the shift based on the determined window
  shift <- windows - window

  # Calculate cutoff retention time based on QC sequence and window
  recut <- length(qcseq) * windows
  valid_idx <- rt < recut
  rtz <- rt[valid_idx]
  mzz <- mz[valid_idx]

  # Calculate sample index and modulo retention time for QC
  qcrt <- rtz %% windows
  sampleidx <- floor(rtz / windows) + 1
  df <- data.frame(
    mz = mzz,
    rt = qcrt,
    rti = rtz,
    sampleidx = sampleidx
  )

  # Separate QC and matrix samples based on qcseq
  dfqc <- subset(df, sampleidx %in% which(qcseq == 1))
  dfqc <- dfqc[stats::complete.cases(dfqc), ]
  dfm <- subset(df, sampleidx %in% which(qcseq == 0))
  dfm <- dfm[stats::complete.cases(dfm), ]

  # Align QC and matrix samples
  xxx <- enviGCMS::getalign2(dfqc$mz, dfqc$rt, ppm = ppm, deltart = deltart)
  qcx <- dfqc[xxx, ]
  yyy <- enviGCMS::getalign2(dfm$mz, dfm$rt, ppm = ppm, deltart = deltart)
  mtx <- dfm[yyy, ]

  # Align QC features
  t_qc <- enviGCMS::getalign(qcx$mz,
                             dfqc$mz,
                             qcx$rt,
                             dfqc$rt,
                             ppm = ppm,
                             deltart = deltart)
  t_qc$idx <- paste(t_qc$mz1, t_qc$rt1)
  t_qc$idx2 <- paste(t_qc$mz2, t_qc$rt2)
  dfqc$idx <- paste(dfqc$mz, dfqc$rt)
  dfqc$mzqc <- t_qc$mz1[match(dfqc$idx, t_qc$idx2)]
  dfqc$rtqc <- t_qc$rt1[match(dfqc$idx, t_qc$idx2)]
  dfqct <- subset(dfqc, stats::complete.cases(dfqc))
  dfqct$idxq <- paste(dfqct$mzqc, dfqct$rtqc)

  # Filter QC features based on minimum sample index occurrences
  z_qc <- stats::aggregate(
    dfqct$sampleidx,
    by = list(dfqct$idxq),
    FUN = function(x)
      length(table(x)) >= minn
  )
  dfqct2 <- dfqct[dfqct$idxq %in% z_qc$Group.1[z_qc$x], ]

  # Align matrix features
  t_mtx <- enviGCMS::getalign(mtx$mz,
                              dfm$mz,
                              mtx$rt,
                              dfm$rt,
                              ppm = ppm,
                              deltart = deltart)
  t_mtx$idx <- paste(t_mtx$mz1, t_mtx$rt1)
  t_mtx$idx2 <- paste(t_mtx$mz2, t_mtx$rt2)
  dfm$idx <- paste(dfm$mz, dfm$rt)
  dfm$mzqc <- t_mtx$mz1[match(dfm$idx, t_mtx$idx2)]
  dfm$rtqc <- t_mtx$rt1[match(dfm$idx, t_mtx$idx2)]
  dfmt <- subset(dfm, stats::complete.cases(dfm))
  dfmt$idxq <- paste(dfmt$mzqc, dfmt$rtqc)

  # Filter matrix features based on minimum QC sample index occurrences
  z_mtx <- stats::aggregate(
    dfmt$sampleidx,
    by = list(dfmt$idxq),
    FUN = function(x)
      length(table(x)) >= sum(qcseq == 0)
  )
  dfmt2 <- dfmt[dfmt$idxq %in% z_mtx$Group.1[z_mtx$x], ]

  # Align QC and matrix data frames to identify duplicates
  re <- enviGCMS::getalign(dfqct2$mz,
                           dfmt2$mz,
                           dfqct2$rt,
                           dfmt2$rt,
                           ppm = ppm,
                           deltart = deltart)

  # Extract QC features that are not aligned with matrix features
  qc <- dfqct2[-re$xid, c('mzqc', 'rtqc', 'sampleidx', 'idxq')]

  # Optimize delta retention time
  optimized_idelta <- getidelta(
    mz = mz,
    rt = rt,
    qcmz = qc$mzqc,
    qcrt = qc$rtqc,
    idelta = idelta,
    shift = shift,
    ppm = ppm,
    deltart = deltart,
    window = recut,
    n = n,
    tol = tol,
    max_iter = max_iter
  )

  return(c(window = windows, idelta = optimized_idelta))
}

#' Generate Quality Control Feature List
#'
#' This function generates a list of features found in Quality Control (QC) samples by aligning QC and matrix samples
#' and filtering based on detection frequency criteria.
#'
#' @param mz Numeric vector of m/z values.
#' @param rt Numeric vector of retention times.
#' @param intensity Numeric vector of intensities corresponding to m/z and rt values.
#' @param idelta Numeric. Initial delta retention time. Default is 60.
#' @param windows Numeric. Retention time window. Default is 600.
#' @param qcseq Integer vector indicating QC samples. Default is c(1, 1, 0, 1, 1, 0, 1, 1, 0).
#' @param deltart Numeric. Tolerance for retention time matching. Default is 5.
#' @param ppm Numeric. Parts per million tolerance for m/z matching. Default is 5.
#' @param minn Integer. Minimum number of QC samples required. Default is 6.
#'
#' @return A data frame containing filtered QC features with the following columns:
#' \itemize{
#'   \item mzqc: aligned m/z of the QC feature.
#'   \item rtqc: aligned retention time of the QC feature.
#'   \item intensity: intensity of the feature in the specific QC sample.
#'   \item sampleidx: index of the QC sample injection.
#'   \item idxq: unique identifier for the QC feature group (mz rt).
#' }
#' The row names of the data frame are set to the sample index (injection number), with suffixes to ensure uniqueness.
#' @examples
#' data(sfi)
#' peak <- find_2d_peaks(mz=sfi$mz,rt=sfi$rt,intensity=sfi$intensity)
#' qc_features <- get_qc_features(peak$mz, peak$rt, peak$intensity,
#'                                idelta=92.25,windows=632.11,minn=6,deltart=10)
#' @param mz Numeric vector of m/z values or an object of class `sfi_peaks`.
#' @param rt Numeric vector of retention times.
#' @param intensity Numeric vector of intensities corresponding to m/z and rt values.
#' @param ... Additional arguments passed to methods.
#' @export
get_qc_features <- function(mz, rt, intensity, ...) {
  UseMethod("get_qc_features")
}

#' @describeIn get_qc_features Method for sfi_peaks object
#' @export
get_qc_features.sfi_peaks <- function(mz, rt = NULL, intensity = NULL, ...) {
  get_qc_features.default(mz$mz, mz$rt, mz$intensity, ...)
}

#' @describeIn get_qc_features Default method for vectors
#' @export
get_qc_features.default <- function(mz,
                    rt,
                    intensity,
                    idelta = 60,
                    windows = 600,
                    qcseq = c(1, 1, 0, 1, 1, 0, 1, 1, 0),
                    deltart = 5,
                    ppm = 5,
                    minn = 6, ...) {
  if (!is.unsorted(rt)) {
    warning("The input retention time vector appears to be sorted. This suggests raw data is being used instead of peak-picked data. Please run find_2d_peaks() first.")
  }
  # Calculate the retention time cutoff
  recut <- length(qcseq) * windows

  # Filter m/z, rt, and intensity within the cutoff
  valid_idx <- rt < recut
  rtz <- rt[valid_idx]
  mzz <- mz[valid_idx]
  intensityz <- intensity[valid_idx]

  # Calculate sample index and modulo retention time for QC
  qcrt <- rtz %% windows
  sampleidx <- floor(rtz / windows) + 1
  df <- data.frame(
    mz = mzz,
    rt = qcrt,
    rti = rtz,
    sampleidx = sampleidx,
    intensity = intensityz
  )

  # Separate QC and matrix samples based on qcseq
  dfqc <- subset(df, sampleidx %in% which(qcseq == 1))
  dfqc <- dfqc[stats::complete.cases(dfqc), ]
  dfm <- subset(df, sampleidx %in% which(qcseq == 0))
  dfm <- dfm[stats::complete.cases(dfm), ]

  # Align QC and matrix samples
  xxx <- enviGCMS::getalign2(dfqc$mz, dfqc$rt, ppm = ppm, deltart = deltart)
  qcx <- dfqc[xxx, ]
  yyy <- enviGCMS::getalign2(dfm$mz, dfm$rt, ppm = ppm, deltart = deltart)
  mtx <- dfm[yyy, ]

  # Align QC features
  t_qc <- enviGCMS::getalign(qcx$mz,
                             dfqc$mz,
                             qcx$rt,
                             dfqc$rt,
                             ppm = ppm,
                             deltart = deltart)
  t_qc$idx <- paste(t_qc$mz1, t_qc$rt1)
  t_qc$idx2 <- paste(t_qc$mz2, t_qc$rt2)
  dfqc$idx <- paste(dfqc$mz, dfqc$rt)
  dfqc$mzqc <- t_qc$mz1[match(dfqc$idx, t_qc$idx2)]
  dfqc$rtqc <- t_qc$rt1[match(dfqc$idx, t_qc$idx2)]
  dfqct <- subset(dfqc, stats::complete.cases(dfqc))
  dfqct$idxq <- paste(dfqct$mzqc, dfqct$rtqc)

  # Filter QC features based on minimum occurrences
  z_qc <- stats::aggregate(
    dfqct$sampleidx,
    by = list(dfqct$idxq),
    FUN = function(x)
      length(table(x)) >= minn
  )
  dfqct2 <- dfqct[dfqct$idxq %in% z_qc$Group.1[z_qc$x], ]

  # Align matrix features
  t_mtx <- enviGCMS::getalign(mtx$mz,
                              dfm$mz,
                              mtx$rt,
                              dfm$rt,
                              ppm = ppm,
                              deltart = deltart)
  t_mtx$idx <- paste(t_mtx$mz1, t_mtx$rt1)
  t_mtx$idx2 <- paste(t_mtx$mz2, t_mtx$rt2)
  dfm$idx <- paste(dfm$mz, dfm$rt)
  dfm$mzqc <- t_mtx$mz1[match(dfm$idx, t_mtx$idx2)]
  dfm$rtqc <- t_mtx$rt1[match(dfm$idx, t_mtx$idx2)]
  dfmt <- subset(dfm, stats::complete.cases(dfm))
  dfmt$idxq <- paste(dfmt$mzqc, dfmt$rtqc)

  # Filter matrix features based on minimum QC sample index occurrences
  z_mtx <- stats::aggregate(
    dfmt$sampleidx,
    by = list(dfmt$idxq),
    FUN = function(x)
      length(table(x)) >= sum(qcseq == 0)
  )
  dfmt2 <- dfmt[dfmt$idxq %in% z_mtx$Group.1[z_mtx$x], ]

  # Align QC and matrix data frames to identify duplicates
  re <- enviGCMS::getalign(dfqct2$mz,
                           dfmt2$mz,
                           dfqct2$rt,
                           dfmt2$rt,
                           ppm = ppm,
                           deltart = deltart)

  # Extract QC features that are not aligned with matrix features
  qc <- dfqct2[-re$xid, c('mzqc', 'rtqc', 'intensity', 'sampleidx', 'idxq')]
  rownames(qc) <- make.unique(as.character(qc$sampleidx))

  # Return the QC data frame
  return(qc)
}

#' Generate Sample Feature Matrix (SFM)
#'
#' This function generates a Sample Feature Matrix (SFM) by aligning and filtering sample peaks against QC peaks.
#' The SFM contains features extracted from individual samples within the single file injection.
#'
#' @param mz Numeric vector of m/z values.
#' @param rt Numeric vector of retention times.
#' @param intensity Numeric vector of intensities corresponding to m/z and rt values.
#' @param idelta Numeric. Initial delta retention time. Default is 60.
#' @param windows Numeric. Retention time window. Default is 600.
#' @param qcseq Integer vector indicating QC samples. Default is c(1, 1, 0, 1, 1, 0, 1, 1, 0).
#' @param deltart Numeric. Tolerance for retention time matching. Default is 5.
#' @param ppm Numeric. Parts per million tolerance for m/z matching. Default is 5.
#' @param minn Integer. Minimum number of QC samples required. Default is 1.
#' @param n Integer. Number of samples for delta optimization. Default is 160.
#'
#' @return A data frame containing the aligned and filtered sample features with the following columns:
#' \itemize{
#'   \item mz: m/z of the feature in the sample.
#'   \item rt: retention time of the feature in the sample (global).
#'   \item srt: relative retention time of the feature within the sample injection window.
#'   \item sampleidx: index of the sample injection.
#'   \item intensity: intensity of the feature.
#'   \item qcmz: m/z of the matching reference QC feature.
#'   \item qcrt: retention time of the matching reference QC feature.
#'   \item shiftrt: absolute difference between sample srt and QC reference retention time.
#'   \item ppmshift: absolute difference in ppm between sample m/z and QC reference m/z.
#' }
#' The row names of the data frame are set to the sample index (injection number).
#' @examples
#' data(sfi)
#' peak <- find_2d_peaks(mz=sfi$mz,rt=sfi$rt,intensity=sfi$intensity)
#' sfm_df <- getsfm(peak$mz, peak$rt, peak$intensity,idelta=92,windows=632,minn=6,n=158,deltart=10)
#' @param mz Numeric vector of m/z values or an object of class `sfi_peaks`.
#' @param rt Numeric vector of retention times.
#' @param intensity Numeric vector of intensities corresponding to m/z and rt values.
#' @param ... Additional arguments passed to methods.
#' @export
getsfm <- function(mz, rt, intensity, ...) {
  UseMethod("getsfm")
}

#' @describeIn getsfm Method for sfi_peaks object
#' @export
getsfm.sfi_peaks <- function(mz, rt = NULL, intensity = NULL, ...) {
  getsfm.default(mz$mz, mz$rt, mz$intensity, ...)
}

#' @describeIn getsfm Default method for vectors
#' @export
getsfm.default <- function(mz,
                   rt,
                   intensity,
                   idelta = 60,
                   windows = 600,
                   qcseq = c(1, 1, 0, 1, 1, 0, 1, 1, 0),
                   deltart = 5,
                   ppm = 5,
                   minn = 1,
                   n = 160, ...) {
  if (!is.unsorted(rt)) {
    warning("The input retention time vector appears to be sorted. This suggests raw data is being used instead of peak-picked data. Please run find_2d_peaks() first.")
  }
  # Calculate the retention time cutoff
  recut <- length(qcseq) * windows

  # Filter m/z, rt, and intensity within the cutoff
  valid_idx <- rt < recut
  rtz <- rt[valid_idx]
  mzz <- mz[valid_idx]
  intensityz <- intensity[valid_idx]

  # Calculate sample index and modulo retention time for QC
  qcrt <- rtz %% windows
  sampleidx <- floor(rtz / windows) + 1
  df <- data.frame(
    mz = mzz,
    rt = qcrt,
    rti = rtz,
    sampleidx = sampleidx,
    intensity = intensityz
  )

  # Separate QC and matrix samples based on qcseq
  dfqc <- subset(df, sampleidx %in% which(qcseq == 1))
  dfqc <- dfqc[stats::complete.cases(dfqc), ]
  dfm <- subset(df, sampleidx %in% which(qcseq == 0))
  dfm <- dfm[stats::complete.cases(dfm), ]

  # Align QC and matrix samples
  xxx <- enviGCMS::getalign2(dfqc$mz, dfqc$rt, ppm = ppm, deltart = deltart)
  qcx <- dfqc[xxx, ]
  yyy <- enviGCMS::getalign2(dfm$mz, dfm$rt, ppm = ppm, deltart = deltart)
  mtx <- dfm[yyy, ]

  # Align QC features
  t_qc <- enviGCMS::getalign(qcx$mz,
                             dfqc$mz,
                             qcx$rt,
                             dfqc$rt,
                             ppm = ppm,
                             deltart = deltart)
  t_qc$idx <- paste(t_qc$mz1, t_qc$rt1)
  t_qc$idx2 <- paste(t_qc$mz2, t_qc$rt2)
  dfqc$idx <- paste(dfqc$mz, dfqc$rt)
  dfqc$mzqc <- t_qc$mz1[match(dfqc$idx, t_qc$idx2)]
  dfqc$rtqc <- t_qc$rt1[match(dfqc$idx, t_qc$idx2)]
  dfqct <- subset(dfqc, stats::complete.cases(dfqc))
  dfqct$idxq <- paste(dfqct$mzqc, dfqct$rtqc)

  # Filter QC features based on minimum occurrences
  z_qc <- stats::aggregate(
    dfqct$sampleidx,
    by = list(dfqct$idxq),
    FUN = function(x)
      length(table(x)) >= minn
  )
  dfqct2 <- dfqct[dfqct$idxq %in% z_qc$Group.1[z_qc$x], ]

  # Align matrix features
  t_mtx <- enviGCMS::getalign(mtx$mz,
                              dfm$mz,
                              mtx$rt,
                              dfm$rt,
                              ppm = ppm,
                              deltart = deltart)
  t_mtx$idx <- paste(t_mtx$mz1, t_mtx$rt1)
  t_mtx$idx2 <- paste(t_mtx$mz2, t_mtx$rt2)
  dfm$idx <- paste(dfm$mz, dfm$rt)
  dfm$mzqc <- t_mtx$mz1[match(dfm$idx, t_mtx$idx2)]
  dfm$rtqc <- t_mtx$rt1[match(dfm$idx, t_mtx$idx2)]
  dfmt <- subset(dfm, stats::complete.cases(dfm))
  dfmt$idxq <- paste(dfmt$mzqc, dfmt$rtqc)

  # Filter matrix features based on minimum QC sample index occurrences
  z_mtx <- stats::aggregate(
    dfmt$sampleidx,
    by = list(dfmt$idxq),
    FUN = function(x)
      length(table(x)) >= sum(qcseq == 0)
  )
  dfmt2 <- dfmt[dfmt$idxq %in% z_mtx$Group.1[z_mtx$x], ]

  # Align QC and matrix data frames to identify duplicates
  re <- enviGCMS::getalign(dfqct2$mz,
                           dfmt2$mz,
                           dfqct2$rt,
                           dfmt2$rt,
                           ppm = ppm,
                           deltart = deltart)

  # Extract QC features that are not aligned with matrix features
  qc <- dfqct2[-re$xid, c('mzqc', 'rtqc', 'intensity', 'sampleidx', 'idxq')]

  # Print summary statistics
  zx <- length(unique(dfmt2$idxq))
  xn <- length(unique(qc$idxq))
  message(zx, " peaks found in all blank samples")
  message(xn, " peaks found in at least ", minn, " QC samples")

  # Extract sample and matrix peaks after the cutoff
  rts <- rt[rt >= recut]
  mzs <- mz[rt >= recut]

  # Initialize data frames for storing results
  mdf <- ndf <- data.frame()

  # Iterate through each sample window to align and filter peaks
  for (i in seq_len(n)) {
    lower_rt <- recut + (i - 1) * idelta
    upper_rt <- lower_rt + windows
    idx <- rt >= lower_rt & rt < upper_rt
    mzi <- mz[idx]
    rti <- rt[idx]
    intensityi <- intensity[idx]
    srt <- rti - (i - 1) * idelta - recut

    # Create data frames for current sample
    dfj <- data.frame(
      mz = mzi,
      rt = rti,
      srt = srt,
      sampleidx = i,
      intensity = intensityi
    )
    dfi <- dfj

    # Align sample with QC
    merge_qc <- suppressMessages(enviGCMS::getalign(mzi,
                                   qc$mzqc,
                                   srt,
                                   qc$rtqc,
                                   ppm = ppm,
                                   deltart = deltart))
    dfi <- dfi[unique(merge_qc$xid), ]
    dfi$qcmz <- merge_qc$mz2[!duplicated(merge_qc$xid)]
    dfi$qcrt <- merge_qc$rt2[!duplicated(merge_qc$xid)]
    dfi <- dfi[!duplicated(paste(dfi$mz, dfi$rt)), ]
    ndf <- rbind(ndf, dfi)

    # Align sample with matrix
    merge_mtx <- suppressMessages(enviGCMS::getalign(mzi,
                                    dfmt2$mzqc,
                                    srt,
                                    dfmt2$rtqc,
                                    ppm = ppm,
                                    deltart = deltart))
    dfj <- dfj[unique(merge_mtx$xid), ]
    dfj$mtmz <- merge_mtx$mz2[!duplicated(merge_mtx$xid)]
    dfj$mtrt <- merge_mtx$rt2[!duplicated(merge_mtx$xid)]
    dfj <- dfj[!duplicated(paste(dfj$mz, dfj$rt)), ]
    mdf <- rbind(mdf, dfj)
  }

  # Calculate shifts and ppm shifts for QC and matrix data
  ndf$shiftrt <- abs(ndf$qcrt - ndf$srt)
  ndf$ppmshift <- abs(ndf$qcmz - ndf$mz) / ndf$qcmz * 1e6

  mdf$shiftrt <- abs(mdf$mtrt - mdf$srt)
  mdf$ppmshift <- abs(mdf$mtmz - mdf$mz) / mdf$mtmz * 1e6

  # Calculate duplicate peaks
  nn <- sum(duplicated(paste(ndf$mz, ndf$rt)) |
              duplicated(paste(ndf$mz, ndf$rt), fromLast = TRUE))
  nnn <- sum(duplicated(paste(mdf$mz, mdf$rt)) |
               duplicated(paste(mdf$mz, mdf$rt), fromLast = TRUE))

  # Calculate unique peaks
  xx <- length(unique(paste(ndf$qcmz, ndf$qcrt)))
  xxx <- length(unique(paste(mdf$mtmz, mdf$mtrt)))

  # Print summary statistics
  message(length(mzs), " sample peaks found")
  message(nrow(ndf), " aligned QC peaks found")
  message("Recover ", (nrow(ndf) - nn) / length(mzs), " peaks")
  message(nrow(mdf), " aligned matrix peaks found")
  message("Recover ", (nrow(mdf) - nnn) / length(mzs), " peaks")

  message(xx, " peaks found in samples")
  message(xxx, " matrix peaks found in samples")
  message(nn, " isomer peaks found in samples")
  message(nnn, " isomer matrix peaks found in samples")

  rownames(ndf) <- make.unique(as.character(ndf$sampleidx))
  return(ndf)
}

#' Feature extraction core function
#'
#' This function finds local max peaks on m/z-retention time 2D plane.
#'
#' @param mz Numeric vector of m/z values.
#' @param rt Numeric vector of retention times.
#' @param intensity Numeric vector of intensities corresponding to m/z and rt values.
#' @param deltart Numeric. Tolerance for retention time matching. Default is 5.
#' @param snr Numeric. signal to ratio to find peaks.
#' @param ppm Numeric. Parts per million tolerance for m/z matching. Default is 5.
#' @param mz_bins Numeric. m/z bins. Default 50000.
#' @param rt_bins Numeric. retention time bins. Default 100.
#'
#' @return A data frame containing m/z, retention time, and intensity of identified peaks.
#' @examples
#' data(sfi)
#' peak <- find_2d_peaks(mz=sfi$mz,rt=sfi$rt,intensity=sfi$intensity)
#' @export
find_2d_peaks <- function(mz,
                          rt,
                          intensity,
                          ppm = 5,
                          deltart = 5,
                          snr = 3.0,
                          mz_bins = NULL,
                          rt_bins = NULL) {
  if (!is.numeric(mz) || !is.numeric(rt) || !is.numeric(intensity)) {
    stop("All inputs must be numeric vectors")
  }
  if (length(mz) != length(rt) || length(mz) != length(intensity)) {
    stop("All input vectors must have the same length")
  }

  if (length(mz) == 0) {
    return(data.frame(mz=numeric(), rt=numeric(), intensity=numeric()))
  }
  
  if (is.null(mz_bins)) {
    max_mz <- max(mz)
    desired_mz_bin_width <- max_mz * ppm * 1e-6
    mz_range <- max(mz) - min(mz)
    if (desired_mz_bin_width > 0 && mz_range > 0) {
      mz_bins <- ceiling(mz_range / desired_mz_bin_width)
    } else {
      mz_bins <- 1
    }
  }
  
  if (is.null(rt_bins)) {
    rt_range <- max(rt) - min(rt)
    if (deltart > 0 && rt_range > 0) {
      rt_bins <- ceiling(rt_range / deltart)
    } else {
      rt_bins <- 1
    }
  }

  result <- find_2d_peaks_c(
    mz,
    rt,
    intensity,
    mz_ppm = ppm,
    rt_window = deltart,
    snr_threshold = snr,
    mz_bins = mz_bins,
    rt_bins = rt_bins
  )
  class(result) <- c("sfi_peaks", "data.frame")
  return(result)
}
