rm(list=ls())

pacman::p_load("tidyverse", "stats", "DescTools", "xcms", "rlang", "stringr",
               install = TRUE)

# 1. Load User Data and Parameters----

# Import the user supplied table of metabolites and conditions

mass.df <- readxl::read_excel("Inputs/Targeted Mass List.xlsx") %>%
  as.data.frame()

parameters.df <- readxl::read_excel("Inputs/Parameters.xlsx") %>%
  as.data.frame()

# Determine the number of metabolites

num.of.metabolites <- nrow(mass.df)

# Get the number of injections

num.of.injections <- parameters.df$number.of.injections[1]

# Generate plot folders for each metabolite

for(m in 1:num.of.metabolites){
  path = paste("Metabolite_Plots/", mass.df$name[m], sep = "")
  dir.create(path = path,
             showWarnings = FALSE)
  rm(list = c("path", "m"))
}

# Create vector of data file names 

data.files <- list.files(path = "mzML Files",
                         full.names = TRUE)

# Add a progress bar

pb <- winProgressBar(title = "Progress Bar", 
                     label = paste("Number of Files Completed: 0 /", length(data.files)), 
                     min = 0,      
                     max = length(data.files), 
                     initial = 0,  
                     width = 300L) 

## Initiate data file for loop ----

#for (d in 1:length(data.files)){
  
  d= 1

  # 2. Prepare Data File ----

  # Make a copy of the data file
  
  file.copy(data.files[d], to = paste(data.files[d], "temp", sep = "_"))
  
  # Read in the copied data file

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
  
  # Unlock "assayData" environment 
  
  env_binding_unlock(run.data@assayData)
  
  # 3. Perform Mass Calibration ----
  
  ## Find mz closest to expected reference masses with greatest intensity to use as the experimental accurate mass ----
  
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
              "max.ref.mass", "min.ref.mass", "mz", "ref.mass", "s", "spectrum"))
  
  # 4. Extract Electropherograms ----
  
  # Define mass error in ppm
  
  mass.error = parameters.df$extraction.window.ppm[1]
  
  # Create a matrix of minimum and maximum m/z values for each metabolite
  
  min <- mass.df$mz - mass.df$mz * mass.error/1000000
  max <- mass.df$mz + mass.df$mz * mass.error/1000000
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
  
  for (m in 1:num.of.metabolites){
    temp.df <- data.frame(electropherograms[m]@intensity)
    colnames(temp.df) <- paste(mass.df[m,1], "intensity", sep = ".")
    eie.df <- cbind(eie.df,temp.df)
  }
  
  # 5. Smooth Intensity Vectors ----
  
  for (m in 2:ncol(eie.df)){ 
    Smooth <- with(eie.df, 
                   ksmooth(x = mt.seconds, 
                           y = eie.df[,m], 
                           kernel = "normal", 
                           bandwidth = parameters.df$eie.smoothing.strength[1]))
    eie.df[,m] <- Smooth[["y"]]
  }
  
  # Clean-up environment
  
  rm(list = c("electropherograms", "mzr", "Smooth", "temp.df", "m", "max", "min"))
  
  # 6. Peak Detection and Integration, and Filtering ----
  
  ## Peak detection ---- 
  
  # Define the minimum number of consecutive increasing and decreasing points a peak must have
  # Define a minimum peak width in seconds
  
  n <- parameters.df$required.points.for.peak.picking[1]
  min.width.cut.off <- mass.df$minimim.peak.width.seconds[1]
  
  for (m in 1:num.of.metabolites){
    
    # Determine the start, apex, and end of peaks.
    # Peaks are intensity features that have n consecutive increases followed by n consecutive decreases
    
    rle.output <- eie.df[,m+1] %>%
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
    
    # Remove peaks that do not meet the minimum peak width cutoff
    
    peak.df <- subset(peak.df, (peak.df$end - peak.df$start) >= min.width.cut.off)
    
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
    
    colnames.vector = c(paste(mass.df[m,1], "start.seconds", sep = "."),
                        paste(mass.df[m,1], "apex.seconds", sep = "."),
                        paste(mass.df[m,1], "end.seconds", sep = "."),
                        paste(mass.df[m,1], "start.intensity", sep = "."),
                        paste(mass.df[m,1], "apex.intensity", sep = "."),
                        paste(mass.df[m,1], "end.intensity", sep = "."),
                        paste(mass.df[m,1], "peak.area", sep = "."))
    
    colnames(peak.df) <- colnames.vector
    
    ## Filter peaks ----
    
    ### Filter peaks based on smallest rmt difference ----
    
    # Determine the expected migration times of the metabolites
    
    rmt.internal.standard <- mass.df$rmt.internal.standard[m]
    
    is.mt.df <- read.table("csv Outputs/Internal Standard Migration Times.csv",
                           header = TRUE,
                           sep = ",",
                           check.names = FALSE) %>%
      as.data.frame()
    
    is.mt.df.start <- seq(from = 1, to = num.of.injections * length(data.files), by = num.of.injections)
    is.mt.df.end <- seq(from = num.of.injections, to = num.of.injections * length(data.files), by = num.of.injections)
    
    is.mt.vec <- is.mt.df[c(is.mt.df.start[d]:is.mt.df.end[d]),rmt.internal.standard]
    
    expected.mt <- mass.df[m,(ncol(mass.df) - num.of.injections + 1):(ncol(mass.df))] %>%
      t() %>%
      as.vector()
    
    expected.mt <- expected.mt * is.mt.vec
    
    # Filter peak.df for peaks within rmt tolerance
    
    rmt.tolerance <- mass.df$rmt.tolerance.percent[m] / 100
    
    for (i in 1:num.of.injections){
      
      peaks <- peak.df %>%
        filter(., peak.df[,2] <= (1 + rmt.tolerance) * expected.mt[i] * 60 & 
                 peak.df[,2] >= (1 - rmt.tolerance) * expected.mt[i] * 60)
      
      # If more than one peak is found choose the nearest one
      
      if(nrow(peaks) > 1){
        peaks <- (peak.df[,2] - expected.mt[i] * 60) %>%
          abs() %>%
          which.min(.)
        peaks <- peak.df[peaks,]
      }
      
      # If no peak is found, generate a place holder
      
      if(nrow(peaks) == 0){
        
        nearest.mt <- (eie.df$mt.seconds - expected.mt[i] * 60) %>%
          abs() %>%
          which.min(.)
        
        peaks <- data.frame(eie.df[nearest.mt,1],
                            eie.df[nearest.mt,1],
                            eie.df[nearest.mt,1],
                            eie.df[nearest.mt,i + 1],
                            eie.df[nearest.mt,i + 1],
                            eie.df[nearest.mt,i + 1],
                            0)
        
        colnames(peaks) <- colnames(peak.df)
      }

      if(i == 1){
        filtered.peaks.df <- peaks
      }else{
        filtered.peaks.df <- rbind(filtered.peaks.df, peaks)
      }
    }
    
    ### Filter duplicated peaks ----
    
    # If the same peak is assigned to multiple injection numbers, reapply rmt filter with more austere rmt tolerances
    # New rmt tolerance will be the original / count, which starts at 2 and increases by 1 each iteration
    
    count = 2
    
    if (any(duplicated(filtered.peaks.df[,2]))){
      
      strict.rmt.tolerance <- rmt.tolerance/count
      
      # find rows with duplicated values 
      
      duplicate.location <- filtered.peaks.df[,2] %>%
        duplicated() %>%
        which()
      
      duplicate.rows <- which(filtered.peaks.df[,2] %in% filtered.peaks.df[duplicate.location,2])
      
      # reapply filtering for these peaks with the more strict rmt tolerance
      
      for (r in duplicate.rows){
        
        peaks <- peak.df %>%
          filter(., peak.df[,2] <= (1 + strict.rmt.tolerance) * expected.mt[r] * 60 & 
                   peak.df[,2] >= (1 - strict.rmt.tolerance) * expected.mt[r] * 60)
        
        # If more than one peak is found choose the nearest one
        
        if(nrow(peaks) > 1){
          peaks <- (peak.df[,2] - expected.mt[r] * 60) %>%
            abs() %>%
            which.min(.)
          peaks <- peak.df[peaks,]
        }
        
        # If no peak is found, generate a place holder
        
        if(nrow(peaks) == 0){
          
          nearest.mt <- (eie.df$mt.seconds - expected.mt[r] * 60) %>%
            abs() %>%
            which.min(.)
          
          peaks <- data.frame(eie.df[nearest.mt,1],
                              eie.df[nearest.mt,1],
                              eie.df[nearest.mt,1],
                              eie.df[nearest.mt,r + 1],
                              eie.df[nearest.mt,r + 1],
                              eie.df[nearest.mt,r + 1],
                              0)
          
          colnames(peaks) <- colnames(peak.df)
        }
        
        filtered.peaks.df[r,] <- peaks
      }
      
      # check if duplicate rows still appear
      
      if (any(duplicated(filtered.peaks.df[,2]))){
        
        duplicate.location <- filtered.peaks.df[,2] %>%
          duplicated() %>%
          which()
        
        duplicate.rows <- which(filtered.peaks.df[,2] %in% filtered.peaks.df[duplicate.location,2])
        
        # restart loop with remaining duplicate rows
        
        r <- duplicate.rows[1]
        
        # Designate a more strict cut off
        
        count = count + 1
        
        strict.rmt.tolerance <- rmt.tolerance/count
    }
    
    ### Filter using peak spaces ----
    
    # Get peak space tolerance
    
    median.space.tol <- mass.df$peak.space.tolerance.percent[m] / 100
    
    # Calculate median peak space
    
    median.space <- filtered.peaks.df[,2] %>%
      diff() %>%
      median()
    
    # Define upper and lower peak space limits
    
    median.space.lower.lim <- median.space - median.space * median.space.tol
    median.space.upper.lim <- median.space + median.space * median.space.tol
    
    # Check if peaks migrate within the tolerance limits
    
    peak.space.tol.check <- between(diff(filtered.peaks.df[,2]), median.space.lower.lim, median.space.upper.lim)
    
    if(all(peak.space.tol.check) != TRUE){
      bad.space <- which(peak.space.tol.check == FALSE)
    }else{
      bad.space <- NA
    }
    
    ## Identify bad peaks, and replace them with peaks meeting peak space criteria
    
    if(is.na(bad.space[1]) == FALSE){
      
      # identify which peaks are potentially incorrectly assigned (bad peaks)
      # these are peaks before and after each bad space
      
      bad.peaks <- c(bad.space, bad.space + 1) %>%
        unique() %>%
        sort()
      
      # define remaining peaks which are correctly assigned (good peaks)
      
      good.peaks <- c(1:num.of.injections) %>%
        setdiff(., c(bad.peaks))
      
      # find the nearest good peak neighbor for each bad peak
      
      for (b in 1:length(bad.peaks)){
        
        # find the nearest good peak neighbor for each bad peak
        
        nearest.good.peak <- (good.peaks - bad.peaks[b]) %>%
          abs() %>%
          which.min()
        
        # calculate the expected migration time 
        
        expected.mt <- filtered.peaks.df[good.peaks[nearest.good.peak], 2] - 
          (good.peaks[nearest.good.peak] - bad.peaks[b]) * median.space
        
        # find peaks nearest to the expected migration time within rmt limits
        
        peaks <- peak.df %>%
          filter(., peak.df[,2] <= (1 + rmt.tolerance) * expected.mt & peak.df[,2] >= (1 - rmt.tolerance) * expected.mt)
        
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
                              eie.df[nearest.mt,i + 1],
                              eie.df[nearest.mt,i + 1],
                              eie.df[nearest.mt,i + 1],
                              0)
          
          colnames(peaks) <- colnames(peak.df)
          
          filtered.peaks.df[bad.peaks[b],] <- peaks
        }
        
        filtered.peaks.df[bad.peaks[b],] <- peaks
      }
    }
    
    # Summarize filtered.peak.df data in metabolite.peak.df
    
    if(m == 1){
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
      mean(temp.noise) + sd(temp.noise)
    }
    
    # Calculate the noise in each region
    
    for (r in 1:length(region.start)){
      temp.noise <- eie.df[region.start[r]:region.end[r], m + 1]
      noise.vec[r] <- noise_calculation(temp.noise)
    }

    # Define the noise as the median noise region
    
    noise <- noise.vec %>%
      sort()
    
    noise <- noise[as.integer(length(noise)/5)]
    
    peak.area.df <- metabolite.peaks.df[,seq(7, ncol(metabolite.peaks.df), 7)]
    
    comment.df[,m] <- ifelse(peak.area.df[,m] < (mass.df$snr.threshold[m] * noise), "<LOD", comment.df[,m])
    
  }
  
  ### Filter interfered peaks ----
  
  for (m in 1:num.of.metabolites){
    
    # Skip metabolites with no reported interference
    
    if(is.na(mass.df$interference[m])){
      next
    }
    
    # Get the name of the interference from mass.df
    
    interference <- mass.df$interference[m] 
    
    # Build an interference data frame since some are metabolites and some are internal standards
    
    interference.df <- cbind(metabolite.peaks.df[,seq(2, ncol(metabolite.peaks.df), 7)],
                             is.mt.df[is.mt.df.start[d]:is.mt.df.end[d],2:ncol(is.mt.df)])
    
    colnames(interference.df)[1:num.of.metabolites] <- mass.df$name
    
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

  
  # 7. Plotting ----
  
  for (m in 1:num.of.metabolites){
    
    ## Create annotation data frame ----
    
    peak.mt.df <- metabolite.peaks.df[,seq(from = 2, to = ncol(metabolite.peaks.df), by = 7)]
    
    # Avoid error caused by the same peak being allocated to two different injections
    
    ann.df <- data.frame("peak.number" = c(1:num.of.injections),
                         "comment" = comment.df[,m],
                         "peak.apex.seconds" = peak.mt.df[,m],
                         "peak.height.counts" = eie.df[which(eie.df$mt.seconds %in% (peak.mt.df[,m])),m+1])
    
    max.peak.height = max(ann.df$peak.height.counts)
    
    ## Create peak fill data frame ----
    
    pf.df <- data.frame("peak.number" = 1,
                        "mt.seconds" = eie.df[,1],
                        "intensity" = eie.df[,m+1])
    
    mt.vec <- vector()
    start.df <- metabolite.peaks.df[,seq(from = 1, to = ncol(metabolite.peaks.df), by = 7)]
    end.df <- metabolite.peaks.df[,seq(from = 3, to = ncol(metabolite.peaks.df), by = 7)]
    
    for (i in 1:num.of.injections){
      
      # Create a migration time vector to track where peaks elute
      
      if(comment.df[i,m] == ""){
        mt.vec.temp <- eie.df$mt.seconds[between(eie.df$mt.seconds, start.df[i,m], end.df[i,m])]
        mt.vec <- append(mt.vec, mt.vec.temp)
        
        # Update peak.number in pf.df
        
        pf.df$peak.number <- ifelse(pf.df$mt.seconds >= start.df[i,m], i, pf.df$peak.number)
        pf.df$peak.number <- as.factor(pf.df$peak.number)
      }else{
        next
      }
    }
    
    pf.df$intensity <- ifelse(pf.df$mt.seconds %in% mt.vec == TRUE, pf.df$intensity , 0)
    
    ## Add baseline intensity
    
    pf.df$baseline <- 0
    
    for (i in 1:num.of.injections){
      if(comment.df[i,m] == ""){
        lower.intensity <- min(c(metabolite.peaks.df[i, m * 7 - 3], metabolite.peaks.df[i, m * 7 - 1]))
        pf.df$baseline <- ifelse(pf.df$mt.seconds >= start.df[i,m] & pf.df$mt.seconds <= end.df[i,m], lower.intensity, pf.df$baseline)
      }else{
        next
      }
    }
    
    ## Plot ----
    
    name <- mass.df$name[m]
    mz <- mass.df$mz[m]
    
    ggplot(data = eie.df) +
      geom_line(aes(x = mt.seconds/60, y = eie.df[,m+1]), colour = "grey50") +
      theme_classic() +
      coord_cartesian(xlim = c(start.df[1,m]/60-2, end.df[num.of.injections,m]/60+2),
                      ylim = c(0, 1.2 * max.peak.height)) +
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
                aes(x = peak.apex.seconds/60,
                    y = peak.height.counts + 0.07 * max.peak.height)) +
      geom_text(data = ann.df,
                label = ann.df$comment,
                aes(x = peak.apex.seconds/60,
                    y = peak.height.counts + 0.11 * max.peak.height)) +
      geom_segment(data = ann.df,
                   aes(x = peak.apex.seconds/60,
                       y = peak.height.counts + 0.04 * max.peak.height,
                       xend = peak.apex.seconds/60,
                       yend = peak.height.counts + 0.01 * max.peak.height),
                   arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
      theme(legend.position = "none")
    
    # Save plots to their respective folders within the "Plots" folder
    
    data.files.name <- list.files(path = "mzML Files")
    data.files.name <- gsub(".mzML", "", data.files.name, fixed = TRUE)
    
    ggsave(filename=paste(name,"_",data.files.name[d],".png",sep=""),
           plot = last_plot(),
           path = paste("Metabolite_Plots/", mass.df$name[m], sep = ""))
    
  }
  
  # 8. Export Data ----
  
  ## Generate peak area data frame ----
  
  peak.area.df <- cbind("file.name" = c(data.files.name[d],2:num.of.injections),
                        "peak.number" = c(1:num.of.injections),
                        metabolite.peaks.df[,seq(from = 7, to = ncol(metabolite.peaks.df), by = 7)])
  peak.area.df$file.name[2:num.of.injections] <- ""
  colnames(peak.area.df)[3:(num.of.metabolites+2)] <- mass.df$name
  
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
                      metabolite.peaks.df[,seq(from = 2, to = ncol(metabolite.peaks.df), by = 7)] / 60)
  peak.mt.df$file.name[2:num.of.injections] <- ""
  colnames(peak.mt.df)[3:(num.of.metabolites+2)] <- mass.df$name
  
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
