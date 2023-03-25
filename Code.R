rm(list=ls())

# Install and load packages

if (!require("pacman")) install.packages("pacman")

pacman::p_load("tidyverse", "stats", "DescTools", "xcms", "rlang", "stringr",
               install = TRUE)

# 1. Load User Data and Parameters----

# Import the user supplied table of metabolites and conditions

mass.df <- readxl::read_excel("Mass List and Parameters.xlsx") %>%
  as.data.frame()

is.df <- readxl::read_excel("Mass List and Parameters.xlsx", sheet = 2) %>%
  as.data.frame()

parameters.df <- readxl::read_excel("Mass List and Parameters.xlsx", sheet = 3) %>%
  as.data.frame()

# Determine the number of metabolites and internal standards

num.of.metabolites <- nrow(mass.df)

num.of.is <- nrow(is.df)

# Get the number of injections

num.of.injections <- parameters.df$number.of.injections[1]

# Create a "Plots" folder to store figures

dir.create(path = "Plots",
           showWarnings = FALSE)

# Generate plot sub-folders for each internal standard and metabolite

name.vec <- c(is.df$name, mass.df$name)

for(n in 1:length(name.vec)){
  path = paste("Plots/",name.vec[n], sep = "")
  dir.create(path = path,
             showWarnings = FALSE)
  rm(list = c("path", "n"))
}

# Create vector of data file names 

data.files <- list.files(path = "mzML Files",
                         full.names = TRUE)

data.file.names <- list.files(path = "mzML Files")

# Add a progress bar

pb <- winProgressBar(title = "Peak Seeker", 
                     label = paste("Number of Files Completed: 0 /", length(data.files)), 
                     min = 0,      
                     max = length(data.files), 
                     initial = 0,  
                     width = 300L) 

## Initiate data file for loop ----

for (d in 1:length(data.files)){
  
  print(paste(d, ". ", "Analyzing Data File: ", data.file.names[d], sep = ""))

  # 2. Prepare Data File ----

  # Make a copy of the data file
  
  file.copy(data.files[d], to = paste(data.files[d], "temp", sep = "_"))
  
  # Read in the copied data file
  
  print("Reading Data File")

  run.data <- readMSData(
    file = paste(data.files[d], "temp", sep = "_"),
    pdata = NULL,
    msLevel = 1,
    verbose = isMSnbaseVerbose(),
    centroided. = FALSE,
    smoothed. = FALSE,
    cache. = 0,
    mode =  "inMemory"
  )
  
  print("File Reading Complete")
  
  # Unlock "assayData" environment 
  
  env_binding_unlock(run.data@assayData)
  
  # 3. Perform Mass Calibration ----
  
  print("Performing Mass Calibration")
  
  ## Find mz closest to expected reference masses with greatest intensity to use as the experimental accurate mass
  
  # Read reference mass values, window to search in, minimum counts to be included
  
  ref.mass <- c(parameters.df$ref.mass.one[1], parameters.df$ref.mass.two[1]) 
  mass.window <- parameters.df$ref.mass.window.ppm[1]
  minimum.counts <- parameters.df$ref.mass.minimum.counts[1]
  
  # Define reference mass ranges to be searched
  
  min.ref.mass <- ref.mass - ref.mass * mass.window/1000000
  max.ref.mass <- ref.mass + ref.mass * mass.window/1000000
  
  # Make a vector of spectrum numbers to be used to reference spectrums 
  
  assaydata.names <- str_pad(1:end(rtime(run.data))[1], 4, pad = "0")
  
  # The correction needs to be applied to each spectrum
  
  for (s in 1:end(rtime(run.data))[1]){
    
    spectrum <- paste("F1.S", assaydata.names[s], sep = "")
    
    # Create vector of mass-to-charges
    
    mz <- run.data@assayData[[spectrum]]@mz
    
    # Create vector of intensities
    
    intensity <- run.data@assayData[[spectrum]]@intensity
    
    ## Lower reference mass correction factor ----
    
    mz.indices <- mz %>%
      between(min.ref.mass[1], max.ref.mass[1]) %>%
      which()
    
    # Ignore spectrums where no reference mass signals are being detected (no spray)
    
    if(length(mz.indices) < 2){
      next
    }
    
    # find mz with max intensity
    
    max.intensity.index <- mz.indices[which.max(intensity[mz.indices])]
    
    # if reference mass intensity is below threshold, do not apply a correction
    
    if(intensity[max.intensity.index] < minimum.counts){
      next
    }
    
    # Build two linear models to predict the true apex of the reference mass
    # Build a model for the left two points
    
    left_line_points <- data.frame("mz" = c(mz[max.intensity.index - 1], mz[max.intensity.index - 2]),
                                   "intensity" = c(intensity[max.intensity.index - 1], intensity[max.intensity.index - 2]))
    
    model_left <- lm(formula = left_line_points$intensity ~ left_line_points$mz)
    
    # Build a model for the right two points
    
    right_line_points <- data.frame("mz" = c(mz[max.intensity.index + 1], mz[max.intensity.index + 2]),
                                    "intensity" = c(intensity[max.intensity.index + 1], intensity[max.intensity.index + 2]))
    
    model_right <- lm(formula = right_line_points$intensity ~ right_line_points$mz)
    
    # Solve for where the two models are equal to each other
    
    slope <- model_left[["coefficients"]][["left_line_points$mz"]] - model_right[["coefficients"]][["right_line_points$mz"]]
    intercept <- model_right[["coefficients"]][["(Intercept)"]] - model_left[["coefficients"]][["(Intercept)"]]
    
    experimental.mz <- solve(slope, intercept)
    
    # determine mass difference and mz index for model
    
    experimental.mass.diff <- experimental.mz - ref.mass[1]
    
    index.of.lower.ref.mass <- max.intensity.index
    
    ## Upper reference mass correction factor ----
    
    mz.indices <- mz %>%
      between(min.ref.mass[2], max.ref.mass[2]) %>%
      which()
    
    # Ignore spectrums where no reference mass signals are being detected (no spray)
    
    if(length(mz.indices) < 2){
      next
    }
    
    # find mz with max intensity
    
    max.intensity.index <- mz.indices[which.max(intensity[mz.indices])]
    
    # if reference mass intensity is below threshold, do not apply a correction
    
    if(intensity[max.intensity.index] < minimum.counts){
      next
    }
    
    # Build two linear models to predict the true apex of the reference mass
    # Build a model for the left two points
    
    left_line_points <- data.frame("mz" = c(mz[max.intensity.index - 1], mz[max.intensity.index - 2]),
                                   "intensity" = c(intensity[max.intensity.index - 1], intensity[max.intensity.index - 2]))
    
    model_left <- lm(formula = left_line_points$intensity ~ left_line_points$mz)
    
    # Build a model for the right two points
    
    right_line_points <- data.frame("mz" = c(mz[max.intensity.index + 1], mz[max.intensity.index + 2]),
                                    "intensity" = c(intensity[max.intensity.index + 1], intensity[max.intensity.index + 2]))
    
    model_right <- lm(formula = right_line_points$intensity ~ right_line_points$mz)
    
    # Solve for where the two models are equal to each other
    
    slope <- model_left[["coefficients"]][["left_line_points$mz"]] - model_right[["coefficients"]][["right_line_points$mz"]]
    intercept <- model_right[["coefficients"]][["(Intercept)"]] - model_left[["coefficients"]][["(Intercept)"]]
    
    experimental.mz <- solve(slope, intercept) 
    
    # determine mass difference and mz index for model
    
    experimental.mass.diff[2] <- experimental.mz - ref.mass[2]
    
    index.of.upper.ref.mass <- max.intensity.index
    
    ## Develop correction model ----
    
    model.data <- data.frame("x" = c(run.data@assayData[[spectrum]]@mz[index.of.lower.ref.mass], run.data@assayData[[spectrum]]@mz[index.of.upper.ref.mass]), 
                             "y" = c(experimental.mass.diff[1], experimental.mass.diff[2]))
    
    model <- lm(y ~ x, model.data)
    
    ## Apply correction model----
    
    correction.vector <- c(model[["coefficients"]][["x"]] * mz[1:length(mz)]) + model[["coefficients"]][["(Intercept)"]]
    
    run.data@assayData[[spectrum]]@mz <- run.data@assayData[[spectrum]]@mz - correction.vector 
    
  }
  
  # clean-up environment 
  
  rm(list = c("model", "model.data", "assaydata.names", "correction.vector", "mz.indices",
              "experimental.mass.diff", "experimental.mz", "index.of.lower.ref.mass", "minimum.counts",
              "index.of.upper.ref.mass", "intensity", "mass.window", "max.intensity.index",
              "max.ref.mass", "min.ref.mass", "mz", "ref.mass", "s", "spectrum", "slope", "intercept",
              "left_line_points", "model_left", "model_right", "right_line_points"))
  
  print("Mass Calibration Complete")
  
  # 4. Extract Electropherograms ----
  
  print("Extracting Electropherograms")
  
  # Define mass error in ppm
  
  mass.error.vec <- c(is.df$extraction.window.ppm, mass.df$extraction.window.ppm)
  
  # Create a matrix of minimum and maximum m/z values for each internal standard and metabolite
  
  mz <- c(is.df$mz, mass.df$mz)
  
  min <- mz - mz * mass.error.vec/1000000
  max <- mz + mz * mass.error.vec/1000000
  mzr <- matrix(c(min, max), ncol = 2)
  
  # Extract electropherograms
  
  electropherograms <- chromatogram(run.data,
                                    mz = mzr,
                                    rt = c(0,end(rtime(run.data))),
                                    aggregationFun = "mean",
                                    missing = 0,
                                    msLevel = 1)
  
  # Create a data frame of migration times and intensities with electropherograms data
  
  eie.df <- data.frame("mt.seconds" = electropherograms[1]@rtime)
  
  for (n in 1:length(name.vec)){
    temp.df <- data.frame(electropherograms[n]@intensity)
    colnames(temp.df) <- paste(name.vec[n], "intensity", sep = " ")
    eie.df <- cbind(eie.df,temp.df)
  }
  
  print("Extraction Complete")
  
  # 5. Smooth Intensity Vectors ----
  
  print("Smoothing Electropherograms")
  
  smoothing.kernal.vec <- c(is.df$smoothing.kernal, mass.df$smoothing.kernal)
  smoothing.strength.vec <- c(is.df$smoothing.strength, mass.df$smoothing.strength)
  
  for (n in 1:length(name.vec)){ 
    Smooth <- with(eie.df, 
                   ksmooth(x = mt.seconds, 
                           y = eie.df[,n + 1], 
                           kernel = smoothing.kernal.vec[n], 
                           bandwidth = smoothing.strength.vec[n]))
    eie.df[,n + 1] <- Smooth[["y"]]
  }
  
  # Clean-up environment
  
  rm(list = c("electropherograms", "mzr", "Smooth", "temp.df", "max", "min", 
              "n", "mass.error.vec", "run.data"))
  
  print("Electropherograms Smoothing Complete")
  
  # 6. Internal Standard Peak Detection, Integration, and Filtering ----
  
  print("Performing Peak Picking and Filtering for Internal Standards")
  
  ## Peak detection ---- 
  
  n <- parameters.df$required.points.for.peak.picking[1]
  
  for (s in 1:num.of.is){
    
    rle.output <- eie.df[,s+1] %>%
      diff() %>%
      sign() %>%
      rle()
    
    consecutive.runs <- which(rle.output$lengths > n & rle.output$values == 1)
    consecutive.runs <- subset(consecutive.runs, (consecutive.runs + 1) %in% (which(rle.output$lengths > n)) == TRUE)
    
    run.lengths <- cumsum(rle.output$lengths) + 1
    
    start <- eie.df$mt.seconds[run.lengths[consecutive.runs - 1]]
    apex <- eie.df$mt.seconds[run.lengths[consecutive.runs]]
    end <- eie.df$mt.seconds[run.lengths[consecutive.runs + 1]]
    
    # For FWHM calculations I will also add intensity values here as well
    
    start.intensity <- eie.df[run.lengths[consecutive.runs - 1], (s+1)]
    apex.intensity <- eie.df[run.lengths[consecutive.runs], (s+1)]
    end.intensity <- eie.df[run.lengths[consecutive.runs + 1], (s+1)]
    
    # Account for peaks that start immediately during the analysis
    
    if(length(start) != length(apex)){
      start <- append(start, 0, 0)
      start.intensity <- append(start.intensity, 0, 0)
    }
    
    # Create a data frame containing the start, apex, and end migration times of each peak
    
    peak.df <- data.frame(start,
                          apex,
                          end,
                          start.intensity,
                          apex.intensity,
                          end.intensity)
    
    ## Migration time filtering ----
    
    # Filter peaks that are outside migration time limits
    
    peak.df <- subset(peak.df, peak.df$start >= is.df$min.rt.min[s] * 60 & peak.df$end <= is.df$max.rt.min[s] * 60)
    
    ## Integrate peaks ----
    
    peak.area.vector = c(1:nrow(peak.df))
    
    for (p in 1:nrow(peak.df)){
      
      peak.area.vector[p] <- AUC(eie.df$mt.seconds,
                                 eie.df[,s+1],
                                 method = "trapezoid",
                                 from = peak.df[p,1],
                                 to = peak.df[p,3],
                                 absolutearea = FALSE,
                                 na.rm = FALSE)
      
      peak.area.vector[p] <- peak.area.vector[p] - (peak.df[p,3] - peak.df[p,1]) * min(peak.df[p,4], peak.df[p,6])
    }
    
    peak.df <- cbind(peak.df, peak.area.vector)
    
    # rename peak.df columns
    
    colnames.vector = c(paste(name.vec[s], "start.seconds", sep = "."),
                        paste(name.vec[s], "apex.seconds", sep = "."),
                        paste(name.vec[s], "end.seconds", sep = "."),
                        paste(name.vec[s], "start.intensity", sep = "."),
                        paste(name.vec[s], "apex.intensity", sep = "."),
                        paste(name.vec[s], "end.intensity", sep = "."),
                        paste(name.vec[s], "peak.area", sep = "."))
    
    colnames(peak.df) <- colnames.vector
    
    # Retain peak.df for future filtering steps
    
    peak.df.fill <- peak.df
    
    ## FWHM filtering ----
    
    # Find the peak intensity at half the peak height
    
    intensity.fwhm <- peak.df[,4] + (peak.df[,5] - peak.df[,4])/2 
    
    # Find the migration times closest to these intensities within each peak
    
    fwhm.vec <- vector()
    single.eie <- eie.df[,c(1,s+1)]
    
    for(p in 1:nrow(peak.df)){
      
      single.eie.temp <- subset(single.eie, single.eie$mt.seconds >= peak.df[p,1] & single.eie$mt.seconds <= peak.df[p,2])
      fwhm.mt.left <- single.eie.temp$mt.seconds[which.min(abs(single.eie.temp[,2] - intensity.fwhm[p]))]
      
      single.eie.temp <- subset(single.eie, single.eie$mt.seconds <= peak.df[p,3] & single.eie$mt.seconds >= peak.df[p,2])
      fwhm.mt.right <- single.eie.temp$mt.seconds[which.min(abs(single.eie.temp[,2] - intensity.fwhm[p]))]
      
      fwhm.vec <- append(fwhm.vec, fwhm.mt.right - fwhm.mt.left)
      
    }
    
    peak.df$fwhm <- fwhm.vec
    
    # Determine the fwhm of the peaks (n = num.of.injections) with the greatest area
    
    df.temp <- peak.df[order(-peak.df[,8]),]
    df.temp <- df.temp[1:num.of.injections,]
    
    fwhm.cutoff <- median(df.temp$fwhm)
    
    peak.df <- subset(peak.df, peak.df$fwhm <= (fwhm.cutoff * is.df$peak.fwhm.tolerance.multiplier[s]))
    peak.df <- peak.df[,c(1:7)]
    
    # subset peak.df so that only the peaks (n = number.of.injections) with the greatest area are kept
    
    cut.off <- sort(peak.df[,7], decreasing = TRUE)[num.of.injections]
    
    peak.df <- subset(peak.df, peak.df[,7] >= cut.off)
    
    ## Peak space filtering ----
    
    # Determine the upper and lower migration time limits for space between peaks
    
    median.space <- peak.df[,2] %>%
      diff() %>%
      median()
    
    median.space.tol <- is.df$peak.space.tolerance.percent[s] / 100
    
    median.space.lower.lim <- median.space - median.space * median.space.tol
    median.space.upper.lim <- median.space + median.space * median.space.tol
    
    # Check if peaks migrate within the tolerance limits
    
    peak.space.tol.check <- between(diff(peak.df[,2]), median.space.lower.lim, median.space.upper.lim)
    
    if(all(peak.space.tol.check) != TRUE){
      bad.space <- which(peak.space.tol.check == FALSE)
    }else{
      bad.space <- NA
    }
    
    ## Scenario 1 ----
    # Only one bad space is detected
    
    if(length(bad.space) == 1 & is.na(bad.space[1]) == FALSE){
      
      false.peak.diff <- peak.df[num.of.injections, 2] - peak.df[(num.of.injections - 1), 2]
      
      # This algorithm always assumes the final peak is false - likely due to carryover
      # It is possible the first peak is false but this seems less likely
      # Final peak only removed if case 1 does not produce a duplicate
      
      # Case 1- An interior peak is missing (usually a blank)
      
      if(bad.space != (num.of.injections - 1)){
        
        # Since the we know the bad space is not at the end, use the space after the bad space to find the expected apex
        
        expected.peak.apex <- peak.df[bad.space[1],2] + peak.df[(bad.space[1] + 2),2] - peak.df[(bad.space[1] + 1),2]
        
      }
      
      # Case 2 -  Final space is false and less than median - suspect that true final peak was missed
      
      if(bad.space == (num.of.injections - 1) & false.peak.diff < median.space.upper.lim){
        expected.peak.apex <- peak.df[(num.of.injections - 1 ),2] + median.space
      }
      
      # Case 3 - Final space is false and greater than median - suspect that peak 1 was missed
      
      if(bad.space == (num.of.injections - 1) & false.peak.diff > median.space.lower.lim){
        expected.peak.apex <- peak.df[1,2] - median.space
      }
      
      peak <- which.min(abs(peak.df.fill[,2] - expected.peak.apex))
      
      # If the nearest peak is too far from the expected migration time use a place holder
      
      if (abs(peak.df.fill[peak,2] - expected.peak.apex) > (median.space / 2)){
        
        nearest.mt <- (eie.df$mt.seconds - expected.peak.apex) %>%
          abs() %>%
          which.min(.)
        
        peaks <- data.frame(eie.df[nearest.mt,1],
                            eie.df[nearest.mt,1],
                            eie.df[nearest.mt,1],
                            eie.df[nearest.mt,s + 1],
                            eie.df[nearest.mt,s + 1],
                            eie.df[nearest.mt,s + 1],
                            0)
        
        colnames(peaks) <- colnames(peak.df)
        
        peak.df <- rbind(peak.df, peaks)
        
        peak.df <- peak.df[order(peak.df[,2]),]
        
      }else{
        
        peak.df <- rbind(peak.df, peak.df.fill[peak,])
        
        peak.df <- peak.df[order(peak.df[,2]),]
        
      }
      
      # Remove bad peak
      
      if (any(duplicated(peak.df[,2]))){
        
        peak.df <- peak.df[-c(which(duplicated(peak.df[,2]))),]
        
      }else{
        
        peak.df <- peak.df[1:(num.of.injections),]
        
      }
    }
    
    ## Scenario 2 ---- 
    #Two or three bad spaces are detected 
    
    if(length(bad.space) == 2 | length(bad.space) == 3){
      
      peaks.to.check <- c(bad.space, (bad.space + 1)) %>%
        unique() %>%
        sort()
      
      # Outside peaks must be checked last
      
      peaks.to.check <- c(peaks.to.check[-1], peaks.to.check[1])
      
      # Loop through each suspect bad peak and remove them iteratively until the number of bad spaces reaches 1
      
      for (p in 1:length(peaks.to.check)){
        
        num.bad.space <- peak.df[,2] %>%
          .[-c(peaks.to.check[p])] %>%
          diff(.) %>%
          between (median.space.lower.lim, median.space.upper.lim) 
        num.bad.space <- length(which(num.bad.space == FALSE))
        
        if (num.bad.space == 1){
          peak.df <- peak.df[-c(peaks.to.check[p]),]
          break
        }
      }
      
      # To avoid errors where removing a peak results in 0 bad spaces only fill
      # in gap if nrow(peak.df) == number of injections - 1
      
      if(nrow(peak.df) == (num.of.injections - 1)){
        
        # Find the peak gap and calculate the expected migration time for the missing peak
        
        gap <- which(between(diff(peak.df[,2]), median.space.lower.lim, median.space.upper.lim) == FALSE)
        
        expected.peak.apex <- (peak.df[(gap[1] + 1),2] - peak.df[gap[1],2])/2 + peak.df[gap[1],2]
        
        # Find the nearest peak in the peak.df.fill data frame to the expected migration time
        # Avoid duplicate peaks by not using exisitng peaks in peak.df
        
        peak.df.fill <- subset(peak.df.fill, !(peak.df.fill[,2] %in% peak.df[,2]))
        
        peak <- which.min(abs(peak.df.fill[,2] - expected.peak.apex))
        peak.df <- rbind(peak.df, peak.df.fill[peak,])
        
        peak.df <- peak.df[order(peak.df[,2]),]
        
      }
    }
    
    # Summarize peak.df data in is.peak.df
    
    if(s == 1){
      is.peaks.df = peak.df
    }else{
      is.peaks.df = cbind(is.peaks.df, peak.df)  
    }
  }
  
  # Make a data frame containing the apex migration times of the internal standards
  # to be used to filter metabolite peaks
  
  is.mt.df <- is.peaks.df[,seq(from = 2,
                               to = ncol(is.peaks.df),
                               by = 7)]
  
  print("Peak Picking and Filtering for Internal Standards Complete")
  
  # 7. Build Relative Migration Time Correction Models ----
  
  # Only use internal standards
  
  is.names.vec <- subset(is.df$name, is.df$class == "Internal Standard") %>%
    paste(., ".apex.seconds", sep = "")
  
  correction.df <- is.mt.df[,is.names.vec]
  
  # Reorder columns in order of peak elution
  
  correction.df <- correction.df[,order(correction.df[1,])] %>%
    suppressWarnings()
  
  ## Model 1 - For analytes with rmts <= 1 ----
  
  # Create two data frames with same dimensions and names as correction.df 
  
  is.rmt.df.1 <- correction.df
  is.mt.diff.df <-  correction.df
  
  # Compute correction values and time ranges
  
  for (c in 2:ncol(correction.df)){
    
    is.rmt.df.1[,c] <- correction.df[,(c - 1)]/correction.df[,c]
    is.mt.diff.df[,c] <- correction.df[,c] - correction.df[,(c - 1)]
    
  }
  
  # The data frame is.rmts.1 is stored for reference during the first data file
  
  if (d == 1){
    
    reference.rmt.df.1 <-  is.rmt.df.1
    
  }
  
  # Compute correction values 
  
  correction.values.df.1 <- (is.rmt.df.1 - reference.rmt.df.1) / is.mt.diff.df
  
  # Assign 0 to first column since no correction is applies to analytes without an internal standard eluting before it
  
  correction.values.df.1[,1] <- 0
  
  # Remove non finite values
  
  is.na(correction.values.df.1)<-sapply(correction.values.df.1, is.infinite)
  
  correction.values.df.1[is.na(correction.values.df.1)] <- 0
  
  ## Model 2 - For analytes with rmts > 1 ----
  
  # Create data frames with same dimensions and names as correction.df
  
  is.rmt.df.2 <- correction.df
  
  # Compute correction values and time ranges
  
  for (c in 1:(ncol(correction.df) - 1)){
    
    is.rmt.df.2[,c] <- correction.df[,(c + 1)]/correction.df[,c]
    is.mt.diff.df[,c] <- correction.df[,(c + 1)] - correction.df[,c]
    
  }
  
  # This is stored for reference during the first data file
  
  if (d == 1){
    
    reference.rmt.df.2 <-  is.rmt.df.2
    
  }
  
  # Compute correction values for rmts > 1
  
  correction.values.df.2 <- (is.rmt.df.2 - reference.rmt.df.2) / is.mt.diff.df
  
  # Assign 0 to last column since no correction is applies to analytes without an internal standard eluting after it
  
  correction.values.df.2[,ncol(correction.values.df.2)] <- 0
  
  # Remove non finite values
  
  is.na(correction.values.df.2)<-sapply(correction.values.df.2, is.infinite)
  
  correction.values.df.2[is.na(correction.values.df.2)] <- 0
  
  # 8. Metabolite  Peak Detection, Integration, and Filtering ----
  
  print("Performing Peak Picking and Filtering for Analytes")
  
  ## Peak detection ---- 
  
  for (m in (num.of.is + 1):length(name.vec)){
    
    peak.df <- data.frame()
      
    # Determine the start, apex, and end of peaks. Use the user defined value "n" to detect peaks.
    # If n results in fewer peaks then injection, decrease n by 1 and repeat
    
    while (nrow(peak.df) < num.of.injections){
      
      n <- parameters.df$required.points.for.peak.picking[1]
      
      rle.output <- eie.df[,m + 1] %>%
        diff() %>%
        sign() %>%
        rle()
      
      consecutive.runs <- which(rle.output$lengths > n & rle.output$values == 1)
      consecutive.runs <- subset(consecutive.runs, (consecutive.runs + 1) %in% (which(rle.output$lengths > n)) == TRUE)
      
      run.lengths <- cumsum(rle.output$lengths) + 1
      
      start <- eie.df$mt.seconds[run.lengths[consecutive.runs - 1]]
      apex <- eie.df$mt.seconds[run.lengths[consecutive.runs]]
      end <- eie.df$mt.seconds[run.lengths[consecutive.runs + 1]]
      
      # I will also add intensity values here as well
      
      start.intensity <- eie.df[run.lengths[consecutive.runs - 1], (m+1)]
      apex.intensity <- eie.df[run.lengths[consecutive.runs], (m+1)]
      end.intensity <- eie.df[run.lengths[consecutive.runs + 1], (m+1)]
      
      # Account for peaks that start immediately during the analysis
      
      if(length(start) != length(apex)){
        start <- append(start, 0, 0)
        start.intensity <- append(start.intensity, 0, 0)
      }
      
      # Create a data frame containing the start, apex, and end migration times of each 
      # peak in addition to required intensities for FWHM calculations
      
      peak.df <- data.frame(start,
                            apex,
                            end,
                            start.intensity,
                            apex.intensity,
                            end.intensity)
      
      n <- n - 1
      
    }
    
    ## Filter peaks by peak width ----
    
    # Define a minimum peak width cut off in seconds. Remove peaks with a width <= cutoff
    # If the cutoff results in fewer peaks than injections, decrease cutoff by 1 and repeat
    
    min.width.cut.off <- mass.df$minimim.peak.width.seconds[m - num.of.is]
    
    peak.df.trim <- subset(peak.df, (peak.df$end - peak.df$start) >= min.width.cut.off)
    
    while (nrow(peak.df.trim) < num.of.injections){
      
      min.width.cut.off <- min.width.cut.off - 1
      
      peak.df.trim <- subset(peak.df, (peak.df$end - peak.df$start) >= min.width.cut.off)
      
    }
    
    peak.df <- peak.df.trim
    
    ## Integrate peaks ----
    
    peak.area.vector = c(1:nrow(peak.df))
    
    for (p in 1:nrow(peak.df)){
      
      peak.area.vector[p] <- AUC(eie.df$mt.seconds,
                                 eie.df[,m+1],
                                 method = "trapezoid",
                                 from = peak.df[p,1],
                                 to = peak.df[p,3],
                                 absolutearea = FALSE,
                                 na.rm = FALSE)
      
      peak.area.vector[p] <- peak.area.vector[p] - (peak.df[p,3] - peak.df[p,1]) * min(peak.df[p,4], peak.df[p,6])
    }
    
    peak.df <- cbind(peak.df, peak.area.vector)
    
    # rename peak.df columns
    
    colnames.vector = c(paste(name.vec[m], "start.seconds", sep = "."),
                        paste(name.vec[m], "apex.seconds", sep = "."),
                        paste(name.vec[m], "end.seconds", sep = "."),
                        paste(name.vec[m], "start.intensity", sep = "."),
                        paste(name.vec[m], "apex.intensity", sep = "."),
                        paste(name.vec[m], "end.intensity", sep = "."),
                        paste(name.vec[m], "peak.area", sep = "."))
    
    colnames(peak.df) <- colnames.vector
    
    ## Filter peaks ----
    
    ### Filter peaks based on smallest rmt difference ----
    
    # Determine the expected migration times of the metabolites
    
    rmt.internal.standard <- paste(mass.df$rmt.internal.standard[m - num.of.is], ".apex.seconds", sep = "")
    
    is.mt.vec <- is.mt.df[,rmt.internal.standard]
    
    rmts <- mass.df[m - num.of.is,(ncol(mass.df) - num.of.injections + 1):(ncol(mass.df))] %>%
      t() %>%
      as.vector()
    
    expected.mt <- rmts * is.mt.vec
    
    #### Apply relative migration time correction ----
    
    for(i in 1:num.of.injections){
      
      # correction only applies to compounds between first and last internal standard
      
      column.index <- which(colnames(correction.values.df.1) == rmt.internal.standard)
      
      if(rmts[i] <= 1 & column.index > 1){
        
        correction.value <- correction.values.df.1[,rmt.internal.standard][i] * (correction.df[,rmt.internal.standard][i] - expected.mt[i])
        
        expected.mt[i] <- (rmts[i] + correction.value) * is.mt.vec[i]
        
      }
      
      column.index <- which(colnames(correction.values.df.2) == rmt.internal.standard)
      
      if(rmts[i] > 1 & column.index < ncol(correction.values.df.2)){
        
        correction.value <- correction.values.df.2[,(column.index + 1)][i] * (expected.mt[i] - correction.df[,rmt.internal.standard][i])
        
        expected.mt[i] <- (rmts[i] + correction.value) * is.mt.vec[i]
      }
      
    }
    
    # Filter peak.df for peaks within rmt tolerance
    
    rmt.tolerance <- mass.df$rmt.tolerance.percent[m - num.of.is] / 100
    
    for (i in 1:num.of.injections){
      
      peaks <- peak.df %>%
        filter(., peak.df[,2] <= (1 + rmt.tolerance) * expected.mt[i] & 
                 peak.df[,2] >= (1 - rmt.tolerance) * expected.mt[i])
      
      # If more than one peak is found choose the nearest one
      
      if(nrow(peaks) > 1){
        peaks <- (peak.df[,2] - expected.mt[i]) %>%
          abs() %>%
          which.min(.)
        peaks <- peak.df[peaks,]
      }
      
      # If no peak is found, generate a place holder
      
      if(nrow(peaks) == 0){
        
        nearest.mt <- (eie.df$mt.seconds - expected.mt[i]) %>%
          abs() %>%
          which.min(.)
        
        peaks <- data.frame(eie.df[nearest.mt,1],
                            eie.df[nearest.mt,1],
                            eie.df[nearest.mt,1],
                            eie.df[nearest.mt,m + 1],
                            eie.df[nearest.mt,m + 1],
                            eie.df[nearest.mt,m + 1],
                            0)
        
        colnames(peaks) <- colnames(peak.df)
      }

      if(i == 1){
        filtered.peaks.df <- peaks
      }else{
        filtered.peaks.df <- rbind(filtered.peaks.df, peaks)
      }
    }
    
    ### Filter peaks outside of run time limits ----
    
    # Set a place holder for peaks where the expected migration time > total run time
    
    total.run.time <- eie.df$mt.seconds[nrow(eie.df)]
    
    late.peaks <- (expected.mt > total.run.time) %>%
      which()
    
    # Find migration times to use as placeholders that do not belong to other identified peaks
    
    mt <- tail(eie.df$mt.seconds, n = 15)
    
    mt <- mt[!(mt %in% filtered.peaks.df[,2])]
    
    for (i in late.peaks){
      mt.temp <- mt[i]
      intensity.temp <- mt
      filtered.peaks.df[i,] <- c(mt.temp,
                                 mt.temp,
                                 mt.temp,
                                 eie.df[which(eie.df$mt.seconds == mt.temp) ,m + 1],
                                 eie.df[which(eie.df$mt.seconds == mt.temp) ,m + 1],
                                 eie.df[which(eie.df$mt.seconds == mt.temp) ,m + 1],
                                 0)
    }
    
    
    ### Filter duplicated peaks ----
    
    # If the same peak is assigned to multiple injection numbers, reapply rmt filter with more austere rmt tolerances
    # New rmt tolerance will be the original / count, which starts at 2 and increases by 1 each iteration
    
    count = 2
    
    while (any(duplicated(filtered.peaks.df[,2]))){
      
      strict.rmt.tolerance <- rmt.tolerance/count
      
      # find rows with duplicated values 
      
      duplicate.location <- filtered.peaks.df[,2] %>%
        duplicated() %>%
        which()
      
      duplicate.rows <- which(filtered.peaks.df[,2] %in% filtered.peaks.df[duplicate.location,2])
      
      # reapply filtering for these peaks with the more strict rmt tolerance
      
      for (r in duplicate.rows){
        
        peaks <- peak.df %>%
          filter(., peak.df[,2] <= (1 + strict.rmt.tolerance) * expected.mt[r] & 
                   peak.df[,2] >= (1 - strict.rmt.tolerance) * expected.mt[r])
        
        # If more than one peak is found choose the nearest one
        
        if(nrow(peaks) > 1){
          peaks <- (peak.df[,2] - expected.mt[r]) %>%
            abs() %>%
            which.min(.)
          peaks <- peak.df[peaks,]
        }
        
        # If no peak is found, generate a place holder
        
        if(nrow(peaks) == 0){
          
          nearest.mt <- (eie.df$mt.seconds - expected.mt[r]) %>%
            abs() %>%
            which.min(.)
          
          peaks <- data.frame(eie.df[nearest.mt,1],
                              eie.df[nearest.mt,1],
                              eie.df[nearest.mt,1],
                              eie.df[nearest.mt,m + 1],
                              eie.df[nearest.mt,m + 1],
                              eie.df[nearest.mt,m + 1],
                              0)
          
          colnames(peaks) <- colnames(peak.df)
        }
        
        filtered.peaks.df[r,] <- peaks
        
      }
      
      # Designate a more strict cut off
      
      count = count + 1
      
      strict.rmt.tolerance <- rmt.tolerance/count
    }
    
    ### Filter using peak spaces ----
    
    # Get peak space tolerance
    
    median.space.tol <- mass.df$peak.space.tolerance.percent[m - num.of.is] / 100
    
    # Calculate median peak space
    
    median.space <- filtered.peaks.df[,2] %>%
      diff() %>%
      median()
    
    # Define upper and lower peak space limits
    
    median.space.lower.lim <- median.space - median.space * median.space.tol
    median.space.upper.lim <- median.space + median.space * median.space.tol
    
    # Check if peaks migrate within the tolerance limits
    
    peak.space.tol.check <- between(diff(filtered.peaks.df[,2]), median.space.lower.lim, median.space.upper.lim)
    
    # Do not include peaks where expected migration time > total run time
    
    peak.space.tol.check <- peak.space.tol.check[-c(late.peaks - 1)]
    
    # Check if peaks migrate within the tolerance limits
    
    peak.space.tol.check <- between(diff(filtered.peaks.df[,2]), median.space.lower.lim, median.space.upper.lim)
    
    if(all(peak.space.tol.check) != TRUE){
      bad.space <- which(peak.space.tol.check == FALSE)
    }else{
      bad.space <- NA
    }
    
    # identify which peaks are potentially incorrectly assigned (bad peaks)
    # these are peaks before and after each bad space
    
    bad.peaks <- c(bad.space, bad.space + 1) %>%
      unique() %>%
      sort()
    
    # Define a count that will be used to modify the peak space tolerance
    
    count = 4
    
    # Identify bad peaks, and replace them with peaks meeting peak space criteria
    # If the number of bad peaks is equal to the number of injections, do not apply this filter
    
    while(length(bad.peaks) > 0 & length(bad.peaks) < num.of.injections){
      
      # define remaining peaks which are correctly assigned (good peaks)
      
      good.peaks <- c(1:num.of.injections) %>%
        setdiff(., c(bad.peaks))
      
      # Use a quarter of the mt difference between good peaks as a tolerance to find the new peaks
      
      peak.tolerance <- filtered.peaks.df[good.peaks,2] %>%
        diff() %>%
        median () / count
      
      # find the nearest good peak neighbor for each bad peak
      
      for (b in 1:length(bad.peaks)){
        
        # find the nearest good peak neighbor for each bad peak
        
        nearest.good.peak <- (good.peaks - bad.peaks[b]) %>%
          abs() %>%
          which.min()
        
        # calculate the expected migration time 
        
        expected.mt <- filtered.peaks.df[good.peaks[nearest.good.peak], 2] - 
          (good.peaks[nearest.good.peak] - bad.peaks[b]) * median.space
        
        # find peaks nearest to the expected migration time within the tolerance
        
        peaks <- peak.df %>%
          filter(., peak.df[,2] <= expected.mt + peak.tolerance & peak.df[,2] >= expected.mt - peak.tolerance)
        
        # if more than one peak is found, select the closest one
        
        if(nrow(peaks) > 1){
          peaks <- (peak.df[,2] - expected.mt) %>%
            abs() %>%
            which.min(.)
          peaks <- peak.df[peaks,]
        }
        
        # if no peaks are found, define a place holder
        
        if(nrow(peaks) == 0){
          
          nearest.mt <- (eie.df$mt.seconds - expected.mt) %>%
            abs() %>%
            which.min(.)
          
          peaks <- data.frame(eie.df[nearest.mt,1],
                              eie.df[nearest.mt,1],
                              eie.df[nearest.mt,1],
                              eie.df[nearest.mt,m + 1],
                              eie.df[nearest.mt,m + 1],
                              eie.df[nearest.mt,m + 1],
                              0)
          
          colnames(peaks) <- colnames(peak.df)
          
          filtered.peaks.df[bad.peaks[b],] <- peaks
        }
        
        # if only one peak is found
        
        filtered.peaks.df[bad.peaks[b],] <- peaks

      }
      
      count = count + 1
      
      duplicate.location <- filtered.peaks.df[,2] %>%
        duplicated() %>%
        which()
      
      # bad peaks correspond to any rows that are not unique
      
      bad.peaks <- which(filtered.peaks.df[,2] %in% filtered.peaks.df[duplicate.location,2])
      
    }
    
    # Summarize filtered.peak.df data in metabolite.peak.df
    
    if(m == (num.of.is + 1)){
      metabolite.peaks.df = filtered.peaks.df
    }else{
      metabolite.peaks.df = cbind(metabolite.peaks.df, filtered.peaks.df)  
    }
  }
  
  ### Filter peaks below LOD ----
  
  # Build a data frame to store comments for each metabolite peak
  
  comment.df <- matrix(nrow = num.of.injections, ncol = num.of.metabolites, "") %>%
    as.data.frame
  
  colnames(comment.df) <- mass.df$name
  
  # Loop through each metabolite and see if its area is below the LOD threshold
  
  for (m in 1:num.of.metabolites){
    
    # Determine the noise of the electropherogram
    # Fine the noise levels in 60 seconds intervals
    
    region.start <- seq(1, nrow(eie.df), 60)
    region.end <- seq(60, nrow(eie.df), 60)
    length(region.start) <- length(region.end)
    
    # Generate a vector to store noise data
    
    noise.vec <- rep(NA, length(region.start))
    
    # Define a function to calculate noise
    
    noise_calculation <- function(temp.noise) {
      mean(temp.noise) + mass.df$snr.threshold[m] * sd(temp.noise)
    }
    
    # Calculate the noise in each region
    
    for (r in 1:length(region.start)){
      temp.noise <- eie.df[region.start[r]:region.end[r], m + num.of.is + 1]
      noise.vec[r] <- noise_calculation(temp.noise)
    }

    # Define the noise as the 20th percentile noise region
    
    noise <- noise.vec %>%
      sort()
    
    noise <- noise[as.integer(length(noise)/5)]
    
    peak.area.df <- metabolite.peaks.df[,seq(7, ncol(metabolite.peaks.df), 7)]
    
    comment.df[,m] <- ifelse(peak.area.df[,m] < noise, "<LOD", comment.df[,m])
    
    ### Annotate injections that are not detected
    
    comment.df[,m] <- ifelse(peak.area.df[,m] == 0, "NPD", comment.df[,m])
    
  }
  
  ### Filter interfered peaks ----
  
  # Build an interference data frame since some are metabolites and some are internal standards
  
  interference.df <- cbind(is.peaks.df[,seq(2, ncol(is.peaks.df), 7)],
                           metabolite.peaks.df[,seq(2, ncol(metabolite.peaks.df), 7)])
  
  for (m in 1:num.of.metabolites){
    
    # Skip metabolites with no reported interference
    
    if(is.na(mass.df$interference[m])){
      next
    }
    
    # Get the names of the interferences from mass.df
    
    interferences <- strsplit(mass.df$interference[m], ", ") %>%
      unlist()
    
    for (k in 1:length(interferences)){
      
      interference <- paste(interferences[k], ".apex.seconds", sep = "")
      
      # See if there is any overlap between the metabolite peak and its interference 
      
      metabolite.name <- paste(mass.df$name[m], ".apex.seconds", sep = "")
      
      for (i in 1:num.of.injections){
        for (j in 1:num.of.injections){
          
          diff.temp <- (metabolite.peaks.df[i,metabolite.name] - interference.df[j,interference]) %>%
            abs() 
          
          comment.df[i,mass.df$name[m]] <- ifelse(diff.temp < mass.df$interference.comigration.threshold.seconds[m], "Interfered", comment.df[i,mass.df$name[m]])
        }
      }
    }
  }
  
  # Combine internal standard and metabolite data frames for plotting
  
  peaks.df <- cbind(is.peaks.df, metabolite.peaks.df)
  
  # update comment data frame account for internal standards
  
  is.comment.df <- matrix(nrow = num.of.injections, ncol = nrow(is.df), "") %>%
    as.data.frame()
  
  colnames(is.comment.df) <- is.df$name
  
  comment.df <- cbind(is.comment.df, comment.df)
  
  print("Peak Picking and Filtering for Analytes Complete")
  
  # 9. Plotting ----
  
  print("Plotting Electropherograms")
  
  for (n in 1:length(name.vec)){
    
    ## Create annotation data frame ----
    
    peak.mt.df <- peaks.df[,seq(from = 2, to = ncol(peaks.df), by = 7)]
      
      ann.df <- data.frame("peak.number" = c(1:num.of.injections),
                           "comment" = comment.df[,n],
                           "peak.apex.seconds" = peak.mt.df[,n],
                           "peak.height.counts" = eie.df[which(eie.df$mt.seconds %in% (peak.mt.df[,n])),n+1])
    
    max.peak.height = max(ann.df$peak.height.counts)
    
    ## Create peak fill data frame ----
    
    pf.df <- data.frame("peak.number" = 1,
                        "mt.seconds" = eie.df[,1],
                        "intensity" = eie.df[,n+1])
    
    mt.vec <- vector()
    start.df <- peaks.df[,seq(from = 1, to = ncol(peaks.df), by = 7)]
    end.df <- peaks.df[,seq(from = 3, to = ncol(peaks.df), by = 7)]
    
    for (i in 1:num.of.injections){
      
      # Create a migration time vector to track where peaks elute
      
      if(comment.df[i,n] == ""){
        mt.vec.temp <- eie.df$mt.seconds[between(eie.df$mt.seconds, start.df[i,n], end.df[i,n])]
        mt.vec <- append(mt.vec, mt.vec.temp)
        
        # Update peak.number in pf.df
        
        pf.df$peak.number <- ifelse(pf.df$mt.seconds >= start.df[i,n], i, pf.df$peak.number)
        pf.df$peak.number <- as.factor(pf.df$peak.number)
      }else{
        next
      }
    }
    
    pf.df$intensity <- ifelse(pf.df$mt.seconds %in% mt.vec == TRUE, pf.df$intensity , 0)
    
    ## Add baseline intensity
    
    pf.df$baseline <- 0
    
    for (i in 1:num.of.injections){
      if(comment.df[i,n] == ""){
        lower.intensity <- min(c(peaks.df[i, n * 7 - 3], peaks.df[i, n * 7 - 1]))
        pf.df$baseline <- ifelse(pf.df$mt.seconds >= start.df[i,n] & pf.df$mt.seconds <= end.df[i,n], lower.intensity, pf.df$baseline)
      }else{
        next
      }
    }
    
    ## Plot ----
    
    mz <- c(is.df$mz, mass.df$mz)
    
    name <- name.vec[n]
    mz <- mz[n]
    
    ggplot(data = eie.df) +
      geom_line(aes(x = mt.seconds/60, y = eie.df[,n+1]), colour = "grey50") +
      theme_classic() +
      coord_cartesian(xlim = c(start.df[1,n]/60-1, end.df[num.of.injections,n]/60+1),
                      ylim = c(min(eie.df[(which(eie.df$mt.seconds == start.df[1,n]) - 60) : (which(eie.df$mt.seconds == end.df[1,n]) + 60),n + 1]) / 3,
                               1.2 * max.peak.height)) +
      scale_y_continuous(name = "Ion Counts",
                         labels = function(x) format(x, scientific = TRUE),
                         expand = c(0,0),
                         breaks = scales::pretty_breaks(n = 10)) +
      scale_x_continuous(name = "Migration Time (Minutes)",
                         breaks = scales::pretty_breaks(n = 10))+
      ggtitle(paste(name, " EIE", " (m/z = ", mz,")",sep = ""),
              subtitle = paste("Data File: ", data.files[d])) +
      geom_ribbon(data = pf.df,
                  aes(x = mt.seconds/60, ymax = intensity, ymin = baseline, fill = peak.number),
                  alpha =0.4) +
      geom_text(data = ann.df,
                label = ann.df$peak.number,
                size  = 5,
                family = "sans",
                aes(x = peak.apex.seconds/60,
                    y = peak.height.counts + 0.07 * max.peak.height)) +
      geom_text(data = ann.df,
                label = ann.df$comment,
                size  = 5,
                family = "sans",
                aes(x = peak.apex.seconds/60,
                    y = peak.height.counts + 0.11 * max.peak.height)) +
      geom_segment(data = ann.df,
                   aes(x = peak.apex.seconds/60,
                       y = peak.height.counts + 0.04 * max.peak.height,
                       xend = peak.apex.seconds/60,
                       yend = peak.height.counts + 0.01 * max.peak.height),
                   arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
      theme(legend.position = "none",
            text = element_text(size = 15, family = "sans"))
    
    # Save plots to their respective folders within the "Plots" folder
    
    data.files.name <- list.files(path = "mzML Files")
    data.files.name <- gsub(".mzML", "", data.files.name, fixed = TRUE)
    
    ggsave(filename=paste(name,"_",data.files.name[d],".png",sep=""),
           width = 16,
           height = 9,
           plot = last_plot(),
           path = paste("Plots/", name.vec[n], sep = ""))
    
  }
  
  print("Plotting Complete")
  
  # 10. Export Data ----
  
  ## Generate peak area data frame ----
  
  peak.area.df <- cbind("file.name" = c(data.files.name[d],2:num.of.injections),
                        "peak.number" = c(1:num.of.injections),
                        peaks.df[,seq(from = 7, to = ncol(peaks.df), by = 7)])
  peak.area.df$file.name[2:num.of.injections] <- ""
  colnames(peak.area.df)[3:(length(name.vec) + 2)] <- name.vec
  
  # Update values to include <LOD and Interfered
  
  for (i in 1:num.of.injections){
    peak.area.df[i,3:ncol(peak.area.df)] <- ifelse(comment.df[i,] == "", peak.area.df[i, 3:ncol(peak.area.df)], comment.df[i,])
  }
  
  if(d == 1){
    peak.area.report = peak.area.df
  }else{
    peak.area.report = rbind(peak.area.report, peak.area.df)
  }
  
  ## Generate peak migration time data frame ----
  
  peak.mt.df <- cbind("file.name" = c(data.files.name[d],2:num.of.injections),
                      "peak.number" = c(1:num.of.injections),
                      peaks.df[,seq(from = 2, to = ncol(peaks.df), by = 7)] / 60)
  peak.mt.df$file.name[2:num.of.injections] <- ""
  colnames(peak.mt.df)[3:(length(name.vec) + 2)] <- name.vec
  
  # Update values to include <LOD and Interfered
  
  for (i in 1:num.of.injections){
    peak.mt.df[i,3:ncol(peak.mt.df)] <- ifelse(comment.df[i,] == "", peak.mt.df[i, 3:ncol(peak.mt.df)], comment.df[i,])
  }
  
  if(d == 1){
    peak.mt.report = peak.mt.df
  }else{
    peak.mt.report = rbind(peak.mt.report, peak.mt.df)
  }
  
  ## Update progress bar ----
  
  progress <- paste(d,"/", length(data.files), "Files Completed")
  setWinProgressBar(pb, d, label = progress) 
  
  # Delete temporary mzml file
  
  file.remove(paste(data.files[d], "temp", sep = "_"))
  
  print(paste("Completed Analysis of Data File: ", data.file.names[d], sep = ""))
  print("")
  
}

## Export data ----

write.csv(peak.area.report,
          file = "csv Outputs/Metabolite Peak Areas.csv",
          row.names = FALSE)

write.csv(peak.mt.report,
          file = "csv Outputs/Metabolite Migration Times.csv",
          row.names = FALSE)

# close progress bar

close(pb)
