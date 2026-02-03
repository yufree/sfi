library(sfi)

test_that("find_2d_peaks works with demo data", {
  data(sfi, package = "sfi")
  peak <- find_2d_peaks(mz = sfi$mz, rt = sfi$rt, intensity = sfi$intensity)
  
  expect_s3_class(peak, "data.frame")
  expect_true(nrow(peak) > 0)
  expect_named(peak, c("mz", "rt", "intensity"))
  
  # Basic sanity check on values
  expect_true(all(peak$mz > 0))
  expect_true(all(peak$rt > 0))
  expect_true(all(peak$intensity > 0))
})

test_that("getsfm returns correctly formatted SFM", {
  data(sfi, package = "sfi")
  # Reduced parameters for faster testing if needed, 
  # but using standard parameters ensures it actually finds something if data is good.
  # We will use the parameters from the vignette.
  peak <- find_2d_peaks(mz = sfi$mz, rt = sfi$rt, intensity = sfi$intensity)
  
  # Vignette params
  idelta <- 92
  windows <- 632
  n <- 100
  deltart <- 10
  minn <- 6
  
  sfm <- getsfm(peak$mz, peak$rt, peak$intensity,
                idelta = idelta, windows = windows, n = n, 
                deltart = deltart, minn = minn)
  
  expect_s3_class(sfm, "data.frame")
  
  # Check expected columns (as defined in documentation update)
  expected_cols <- c("mz", "rt", "srt", "sampleidx", "intensity", 
                     "qcmz", "qcrt", "shiftrt", "ppmshift")
  # Note: getsfm might return extra columns, but should contain these.
  # Based on code: c('mz', 'rt', 'srt', 'sampleidx', 'intensity', 'qcmz', 'qcrt', 'shiftrt', 'ppmshift')
  # Checking intersection to allow for potential extras if any (though usually strict is better).
  expect_true(all(expected_cols %in% colnames(sfm)))
  
  # Check types
  expect_type(sfm$mz, "double")
  expect_true(is.numeric(sfm$sampleidx)) # Accepts integer or double
})

test_that("get_sfi_params returns expected window and delta", {
  data(sfi, package = "sfi")
  peak <- find_2d_peaks(mz = sfi$mz, rt = sfi$rt, intensity = sfi$intensity)
  
  # Use a subset of high intensity peaks for speed/stability in test as per vignette "Advanced Usage"
  sfmsub <- peak[peak$intensity > 1e4, ]
  
  # The test might fail if optimization doesn't converge, so we wrap safely or check for non-null
  # However, with demo data it should work.
  sfi_params <- get_sfi_params(sfmsub$mz, sfmsub$rt, sfmsub$intensity, 
                     n = 158, deltart = 5)
  
  expect_type(sfi_params, "double")
  expect_length(sfi_params, 2)
  expect_named(sfi_params, c("window", "idelta"))
  expect_true(sfi_params["window"] > 0) # Window
  expect_true(sfi_params["idelta"] > 0) # Delta
})

test_that("get_qc_features returns QC feature list", {
  data(sfi, package = "sfi")
  peak <- find_2d_peaks(mz = sfi$mz, rt = sfi$rt, intensity = sfi$intensity)
  
  qcdf <- get_qc_features(peak$mz, peak$rt, peak$intensity,
                  idelta = 92.25, windows = 632.11, minn = 6, deltart = 10)
  
  expect_s3_class(qcdf, "data.frame")
  expect_true(all(c("mzqc", "rtqc", "intensity", "sampleidx", "idxq") %in% colnames(qcdf)))
})

test_that("getsff clusters and pairs features", {
  data(sfi, package = "sfi")
  peak <- find_2d_peaks(mz = sfi$mz, rt = sfi$rt, intensity = sfi$intensity)
  
  # Use a small subset for getsff to avoid combinatorial explosion/slow test
  # Select a small slice of RT
  idx <- peak$rt > 600 & peak$rt < 700
  sub_peak <- peak[idx, ]
  
  if(nrow(sub_peak) > 1) {
    sff <- getsff(sub_peak$mz, sub_peak$rt)
    expect_s3_class(sff, "data.frame")
    if(nrow(sff) > 0) {
      expect_named(sff, c("mz1", "rt1", "mz2", "rt2", "pmr", "pmd"))
    }
  }
})
