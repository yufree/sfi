# Load necessary libraries
library(proxy)
library(enviGCMS)
library(stats)
library(Rcpp)

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
#' @return A data frame containing paired m/z and retention time values with their differences.
#' @export
getsff <- function(mz, rt, ppm = 5, minn = 2, refmz = NULL) {
  # Calculate Manhattan distance matrix for m/z values
  dis <- proxy::dist(mz, method = "manhattan")

  # Perform hierarchical clustering
  fit <- stats::hclust(dis)

  # Cut the dendrogram into clusters based on the specified height
  mzcluster <- stats::cutree(fit, h = 0.001 * ppm)

  # Identify clusters with at least 'minn' members
  valid_clusters <- names(which(table(mzcluster) > minn))
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

    # Calculate Manhattan distance for retention times within the cluster
    rtd <- proxy::dist(bin$rt, method = "manhattan")

    # Extract lower triangle indices
    lower_indices <- which(lower.tri(as.matrix(rtd)), arr.ind = TRUE)

    # Pair retention times and m/z values
    rt1 <- bin$rt[lower_indices[, 1]]
    rt2 <- bin$rt[lower_indices[, 2]]
    mz1 <- bin$mz[lower_indices[, 1]]
    mz2 <- bin$mz[lower_indices[, 2]]

    # Calculate differences
    pmr <- as.numeric(rtd)
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
#' @export
getwindow <- function(mz,
                      rt,
                      lower = 620,
                      upper = 650,
                      ppm = 5,
                      minn = 1,
                      qcseq = c(1, 1, 0, 1, 1, 0, 1, 1, 0)) {
  # Calculate the cutoff retention time based on QC sequence length and upper window
  recut <- length(qcseq) * upper

  # Filter m/z and rt values within the cutoff
  valid_idx <- rt < recut
  mz_filtered <- mz[valid_idx]
  rt_filtered <- rt[valid_idx]

  # Generate paired m/z and rt data
  sff <- getsff(mz_filtered, rt_filtered, ppm = ppm, minn = minn)

  # Determine the most frequent rounded pmr values corresponding to QC samples
  rtc <- as.numeric(names(table(round(sff$pmr, 2)))[
    order(table(round(sff$pmr, 2)), decreasing = TRUE)
  ])[1:length(qcseq)]

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

  message(paste("Window is", window))
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
#' @export
getidelta <- function(mz,
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
                      max_iter = 100) {
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
    sum(abs(shiftrt) < deltart & sampleidx %in% seq_len(n), na.rm = TRUE)
  }

  # Binary search optimization to find the best idelta
  binary_search_optimization <- function(f, lower, upper, irt, n, deltart, tol = 0.01, max_iter = 100) {
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

  message(paste("Delta retention time is", optimized_idelta))
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
#' @return A numeric vector containing the optimal window and delta retention time.
#' @export
getqc <- function(mz,
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
                  wupper = 650) {
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
  t_qc <- enviGCMS::getalign(qcx$mz, dfqc$mz, qcx$rt, dfqc$rt, ppm = ppm, deltart = deltart)
  t_qc$idx <- paste(t_qc$mz1, t_qc$rt1)
  t_qc$idx2 <- paste(t_qc$mz2, t_qc$rt2)
  dfqc$idx <- paste(dfqc$mz, dfqc$rt)
  dfqc$mzqc <- t_qc$mz1[match(dfqc$idx, t_qc$idx2)]
  dfqc$rtqc <- t_qc$rt1[match(dfqc$idx, t_qc$idx2)]
  dfqct <- subset(dfqc, stats::complete.cases(dfqc))
  dfqct$idxq <- paste(dfqct$mzqc, dfqct$rtqc)

  # Filter QC features based on minimum sample index occurrences
  z_qc <- stats::aggregate(dfqct$sampleidx, by = list(dfqct$idxq), FUN = function(x) length(table(x)) >= minn)
  dfqct2 <- dfqct[dfqct$idxq %in% z_qc$Group.1[z_qc$x], ]

  # Align matrix features
  t_mtx <- enviGCMS::getalign(mtx$mz, dfm$mz, mtx$rt, dfm$rt, ppm = ppm, deltart = deltart)
  t_mtx$idx <- paste(t_mtx$mz1, t_mtx$rt1)
  t_mtx$idx2 <- paste(t_mtx$mz2, t_mtx$rt2)
  dfm$idx <- paste(dfm$mz, dfm$rt)
  dfm$mzqc <- t_mtx$mz1[match(dfm$idx, t_mtx$idx2)]
  dfm$rtqc <- t_mtx$rt1[match(dfm$idx, t_mtx$idx2)]
  dfmt <- subset(dfm, stats::complete.cases(dfm))
  dfmt$idxq <- paste(dfmt$mzqc, dfmt$rtqc)

  # Filter matrix features based on minimum QC sample index occurrences
  z_mtx <- stats::aggregate(dfmt$sampleidx, by = list(dfmt$idxq), FUN = function(x) length(table(x)) >= sum(qcseq == 0))
  dfmt2 <- dfmt[dfmt$idxq %in% z_mtx$Group.1[z_mtx$x], ]

  # Align QC and matrix data frames to identify duplicates
  re <- enviGCMS::getalign(dfqct2$mz, dfmt2$mz, dfqct2$rt, dfmt2$rt, ppm = ppm, deltart = deltart)

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

  return(c(windows, optimized_idelta))
}

#' Generate Quality Control Data Frame
#'
#' This function generates a QC data frame by aligning QC and matrix samples and filtering based on criteria.
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
#' @return A data frame containing QC information after alignment and filtering.
#' @export
getqcdf <- function(mz,
                    rt,
                    intensity,
                    idelta = 60,
                    windows = 600,
                    qcseq = c(1, 1, 0, 1, 1, 0, 1, 1, 0),
                    deltart = 5,
                    ppm = 5,
                    minn = 6) {
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
  t_qc <- enviGCMS::getalign(qcx$mz, dfqc$mz, qcx$rt, dfqc$rt, ppm = ppm, deltart = deltart)
  t_qc$idx <- paste(t_qc$mz1, t_qc$rt1)
  t_qc$idx2 <- paste(t_qc$mz2, t_qc$rt2)
  dfqc$idx <- paste(dfqc$mz, dfqc$rt)
  dfqc$mzqc <- t_qc$mz1[match(dfqc$idx, t_qc$idx2)]
  dfqc$rtqc <- t_qc$rt1[match(dfqc$idx, t_qc$idx2)]
  dfqct <- subset(dfqc, stats::complete.cases(dfqc))
  dfqct$idxq <- paste(dfqct$mzqc, dfqct$rtqc)

  # Filter QC features based on minimum occurrences
  z_qc <- stats::aggregate(dfqct$sampleidx, by = list(dfqct$idxq), FUN = function(x) length(table(x)) >= minn)
  dfqct2 <- dfqct[dfqct$idxq %in% z_qc$Group.1[z_qc$x], ]

  # Align matrix features
  t_mtx <- enviGCMS::getalign(mtx$mz, dfm$mz, mtx$rt, dfm$rt, ppm = ppm, deltart = deltart)
  t_mtx$idx <- paste(t_mtx$mz1, t_mtx$rt1)
  t_mtx$idx2 <- paste(t_mtx$mz2, t_mtx$rt2)
  dfm$idx <- paste(dfm$mz, dfm$rt)
  dfm$mzqc <- t_mtx$mz1[match(dfm$idx, t_mtx$idx2)]
  dfm$rtqc <- t_mtx$rt1[match(dfm$idx, t_mtx$idx2)]
  dfmt <- subset(dfm, stats::complete.cases(dfm))
  dfmt$idxq <- paste(dfmt$mzqc, dfmt$rtqc)

  # Filter matrix features based on minimum QC sample index occurrences
  z_mtx <- stats::aggregate(dfmt$sampleidx, by = list(dfmt$idxq), FUN = function(x) length(table(x)) >= sum(qcseq == 0))
  dfmt2 <- dfmt[dfmt$idxq %in% z_mtx$Group.1[z_mtx$x], ]

  # Align QC and matrix data frames to identify duplicates
  re <- enviGCMS::getalign(dfqct2$mz, dfmt2$mz, dfqct2$rt, dfmt2$rt, ppm = ppm, deltart = deltart)

  # Extract QC features that are not aligned with matrix features
  qc <- dfqct2[-re$xid, c('mzqc', 'rtqc', 'intensity', 'sampleidx', 'idxq')]

  # Return the QC data frame
  return(qc)
}

#' Generate Quality Control Matrix (SFM)
#'
#' This function generates a Quality Control Matrix by aligning and filtering sample and matrix peaks.
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
#' @return A data frame containing the aligned and filtered Quality Control Matrix.
#' @export
getsfm <- function(mz,
                   rt,
                   intensity,
                   idelta = 60,
                   windows = 600,
                   qcseq = c(1, 1, 0, 1, 1, 0, 1, 1, 0),
                   deltart = 5,
                   ppm = 5,
                   minn = 1,
                   n = 160) {
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
  t_qc <- enviGCMS::getalign(qcx$mz, dfqc$mz, qcx$rt, dfqc$rt, ppm = ppm, deltart = deltart)
  t_qc$idx <- paste(t_qc$mz1, t_qc$rt1)
  t_qc$idx2 <- paste(t_qc$mz2, t_qc$rt2)
  dfqc$idx <- paste(dfqc$mz, dfqc$rt)
  dfqc$mzqc <- t_qc$mz1[match(dfqc$idx, t_qc$idx2)]
  dfqc$rtqc <- t_qc$rt1[match(dfqc$idx, t_qc$idx2)]
  dfqct <- subset(dfqc, stats::complete.cases(dfqc))
  dfqct$idxq <- paste(dfqct$mzqc, dfqct$rtqc)

  # Filter QC features based on minimum occurrences
  z_qc <- stats::aggregate(dfqct$sampleidx, by = list(dfqct$idxq), FUN = function(x) length(table(x)) >= minn)
  dfqct2 <- dfqct[dfqct$idxq %in% z_qc$Group.1[z_qc$x], ]

  # Align matrix features
  t_mtx <- enviGCMS::getalign(mtx$mz, dfm$mz, mtx$rt, dfm$rt, ppm = ppm, deltart = deltart)
  t_mtx$idx <- paste(t_mtx$mz1, t_mtx$rt1)
  t_mtx$idx2 <- paste(t_mtx$mz2, t_mtx$rt2)
  dfm$idx <- paste(dfm$mz, dfm$rt)
  dfm$mzqc <- t_mtx$mz1[match(dfm$idx, t_mtx$idx2)]
  dfm$rtqc <- t_mtx$rt1[match(dfm$idx, t_mtx$idx2)]
  dfmt <- subset(dfm, stats::complete.cases(dfm))
  dfmt$idxq <- paste(dfmt$mzqc, dfmt$rtqc)

  # Filter matrix features based on minimum QC sample index occurrences
  z_mtx <- stats::aggregate(dfmt$sampleidx, by = list(dfmt$idxq), FUN = function(x) length(table(x)) >= sum(qcseq == 0))
  dfmt2 <- dfmt[dfmt$idxq %in% z_mtx$Group.1[z_mtx$x], ]

  # Align QC and matrix data frames to identify duplicates
  re <- enviGCMS::getalign(dfqct2$mz, dfmt2$mz, dfqct2$rt, dfmt2$rt, ppm = ppm, deltart = deltart)

  # Extract QC features that are not aligned with matrix features
  qc <- dfqct2[-re$xid, c('mzqc', 'rtqc', 'intensity', 'sampleidx', 'idxq')]

  # Print summary statistics
  zx <- length(unique(dfmt2$idxq))
  xn <- length(unique(qc$idxq))
  message(paste(zx, "peaks found in all blank samples"))
  message(paste(xn, "peaks found in at least", minn, "QC samples"))

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
    merge_qc <- enviGCMS::getalign(mzi, qc$mzqc, srt, qc$rtqc, ppm = ppm, deltart = deltart)
    dfi <- dfi[unique(merge_qc$xid), ]
    dfi$qcmz <- merge_qc$mz2[!duplicated(merge_qc$xid)]
    dfi$qcrt <- merge_qc$rt2[!duplicated(merge_qc$xid)]
    dfi <- dfi[!duplicated(paste(dfi$mz, dfi$rt)), ]
    ndf <- rbind(ndf, dfi)

    # Align sample with matrix
    merge_mtx <- enviGCMS::getalign(mzi, dfmt2$mzqc, srt, dfmt2$rtqc, ppm = ppm, deltart = deltart)
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
  nn <- sum(duplicated(paste(ndf$mz, ndf$rt)) | duplicated(paste(ndf$mz, ndf$rt), fromLast = TRUE))
  nnn <- sum(duplicated(paste(mdf$mz, mdf$rt)) | duplicated(paste(mdf$mz, mdf$rt), fromLast = TRUE))

  # Calculate unique peaks
  xx <- length(unique(paste(ndf$qcmz, ndf$qcrt)))
  xxx <- length(unique(paste(mdf$mtmz, mdf$mtrt)))

  # Print summary statistics
  message(paste(length(mzs), 'sample peaks found'))
  message(paste(nrow(ndf), "aligned QC peaks found"))
  message(paste("Recover", (nrow(ndf) - nn) / length(mzs), "peaks"))
  message(paste(nrow(mdf), "aligned matrix peaks found"))
  message(paste("Recover", (nrow(mdf) - nnn) / length(mzs), "peaks"))

  message(paste(xx, "peaks found in samples"))
  message(paste(xxx, "matrix peaks found in samples"))
  message(paste(nn, "isomer peaks found in samples"))
  message(paste(nnn, "isomer matrix peaks found in samples"))

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
#' @return A data frame containing the aligned and filtered Quality Control Matrix.
#' @export
find_2d_peaks <- function(mz, rt, intensity,
                          ppm = 5,
                          deltart = 5,
                          snr = 3.0,
                          mz_bins = 50000,
                          rt_bins = 100) {

  if (!is.numeric(mz) || !is.numeric(rt) || !is.numeric(intensity)) {
    stop("All inputs must be numeric vectors")
  }
  if (length(mz) != length(rt) || length(mz) != length(intensity)) {
    stop("All input vectors must have the same length")
  }

  find_2d_peaks_c(mz, rt, intensity,
                  mz_ppm=ppm, rt_window=deltart, snr_threshold = snr,mz_bins = mz_bins,rt_bins = rt_bins)
}
