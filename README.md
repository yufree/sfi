# sfi

## Introduction

The `sfi` package provides features extraction tools for single file injections(sfi) mode Gas/liquid chromatography-mass spectrometry (GC/LC-MS) data. 

## Installation

```
# Install from CRAN (when available)
install.packages("sfi")

# Or install the development version from GitHub
remotes::install_github("yufree/sfi")
```

## Basic Usage

### Loading the Package

```
library(sfi)
```

### Loading raw data

```
library(mzR)
path <- 'sfi.mzML'
mzml_file <- openMSfile(path)
```

### Feature extraction

```
# read meta data
tt <- header(mzml_file)
# generate retention time vector
rt <- rep(tt$retentionTime,tt$peaksCount)
# extract peaks
peaks <- peaks(mzml_file)
peak <- do.call(rbind,peaks)
peak <- cbind(peak,rt)
# perform peaks picking
peaklist <- find_2d_peaks(peak[,1],peak[,3],peak[,2],rt_window = 30,mz_bins = 40000,rt_bins = 2000)
```

## Feature alignment

The `getsfm()` function will extract sample features from one injections file.

```
mz <- peaklist$mz
rt <- peaklist$rt
intensity <- peaklist$intensity
# injection interval
idelta <- 80
# time windows for a full seperation
windows <- 900
# sample numbers in the files
n <- 100
# retention time shift in seconds
deltart <- 10
# min peak number in pooled qc samples 
minn <- 6
# align peaks from sfi
ndf <- getsfm(mz,rt,intensity,idelta=idelta,windows=windows,n=100,deltart = 10,minn = 6)
# generate feature list
ndf$peakid <- paste0(round(ndf$qcmz,4),'@',round(ndf$qcrt))
nndf <- ndf[,c('sampleidx','intensity','peakid')]
library(data.table)
feature <- dcast(setDT(nndf), ... ~ sampleidx, value.var = "intensity",fun.aggregate = sum)
```

Save feature list as csv file.

```
# save feature as csv file
fwrite(feature,'featurelist.csv')
```
