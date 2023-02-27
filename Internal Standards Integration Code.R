rm(list=ls())

# Install and load packages

if (!require("pacman")) install.packages("pacman")

pacman::p_load("tidyverse", "stats", "DescTools", "xcms", "rlang", "stringr",
               install = TRUE)

# 1. Load User Data and Parameters----

# Import the user supplied table of internal standards and conditions

is.df <- readxl::read_excel("Inputs/Internal Standard Mass List.xlsx") %>%
  as.data.frame()

parameters.df <- readxl::read_excel("Inputs/Parameters.xlsx") %>%
  as.data.frame()

# Determine the number of internal standards

num.of.iss <- nrow(is.df)

# Get the number of injections

num.of.injections <- parameters.df$number.of.injections[1] 

# Generate plot folders for each internal standard

for(i in 1:num.of.iss){
  path = paste("Internal Standard Plots/", is.df$name[i], sep = "")
  dir.create(path = path,
             showWarnings = FALSE)
}

# Create vector of mzML file names

data.files <- list.files(path = "mzML Files",
                         full.names = TRUE)

# Add a progress bar

pb <- winProgressBar(title = "Peak Seeker", 
                     label = paste("Number of Files Completed: 0 /", length(data.files)), 
                     min = 0,      
                     max = length(data.files), 
                     initial = 0,  
                     width = 300L) 

## Initiate data file for loop ----

for (d in 1:length(data.files)){

  # 2. Prepare Data File----
  
  file.copy(data.files[d], to = paste(data.files[d], "temp", sep = "_"))
  
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
  
  # Unlock assayData environment 
  
  env_binding_unlock(run.data@assayData)
  
  # 3. Perform Mass Calibration ----
  
  # Find mz closest to expected reference masses with greatest intensity to use as the experimental accurate mass
  
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
  
  # Create a matrix of minimum and maximum m/z values for each internal standard
  
  min <- is.df$mz - is.df$mz * mass.error/1000000
  max <- is.df$mz + is.df$mz * mass.error/1000000
  mzr <- matrix(c(min, max), ncol = 2)
  
  electropherograms <- chromatogram(run.data,
                                    mz = mzr,
                                    rt = c(0,end(rtime(run.data))),
                                    aggregationFun = "mean",
                                    missing = 0,
                                    msLevel = 1)
  
  # Create a data frame of migration times and intensities (electropherogram data)
  
  eie.df <- data.frame("mt.seconds" = electropherograms[1]@rtime)
  
  for (i in 1:num.of.iss){
    temp.df <- data.frame(electropherograms[i]@intensity)
    colnames(temp.df) <- paste(is.df$name[i], "intensity", sep = ".")
    eie.df <- cbind(eie.df,temp.df)
  }
  
  # 5. Smooth intensity vectors ----
  
  for (j in 2:ncol(eie.df)){ 
    Smooth <- with(eie.df, 
                   ksmooth(x = mt.seconds, 
                           y = eie.df[,j], 
                           kernel = "normal", 
                           bandwidth = parameters.df$eie.smoothing.strength[1]))
    eie.df[,j] <- Smooth[["y"]]
  }
  
  # 6. Peak Detection, Integration, and Filtering ----
  
  ## Peak detection ---- 
  
  n <- parameters.df$required.points.for.peak.picking[1]
  
  for (i in 1:num.of.iss){
    
    rle.output <- eie.df[,i+1] %>%
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
    
    start.intensity <- eie.df[run.lengths[consecutive.runs - 1], (i+1)]
    apex.intensity <- eie.df[run.lengths[consecutive.runs], (i+1)]
    end.intensity <- eie.df[run.lengths[consecutive.runs + 1], (i+1)]
    
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
    
    peak.df <- subset(peak.df, peak.df$start >= is.df$min.rt.min[i] * 60 & peak.df$end <= is.df$max.rt.min[i] * 60)
    
    ## Integrate peaks ----
    
    peak.area.vector = c(1:nrow(peak.df))
    
    for (p in 1:nrow(peak.df)){
      
      peak.area.vector[p] <- AUC(eie.df$mt.seconds,
                                 eie.df[,i+1],
                                 method = "trapezoid",
                                 from = peak.df[p,1],
                                 to = peak.df[p,3],
                                 absolutearea = FALSE,
                                 na.rm = FALSE)
      
      peak.area.vector[p] <- peak.area.vector[p] - (peak.df[p,3] - peak.df[p,1]) * min(peak.df[p,4], peak.df[p,6])
    }
    
    peak.df <- cbind(peak.df, peak.area.vector)
    
    # rename peak.df columns
    
    colnames.vector = c(paste(is.df[i,1], "start.seconds", sep = "."),
                        paste(is.df[i,1], "apex.seconds", sep = "."),
                        paste(is.df[i,1], "end.seconds", sep = "."),
                        paste(is.df[i,1], "start.intensity", sep = "."),
                        paste(is.df[i,1], "apex.intensity", sep = "."),
                        paste(is.df[i,1], "end.intensity", sep = "."),
                        paste(is.df[i,1], "peak.area", sep = "."))
    
    colnames(peak.df) <- colnames.vector
    
    # Retain peak.df for future filtering steps
    
    peak.df.fill <- peak.df
    
    ## FWHM filtering ----
    
    # Find the peak intensity at half the peak height
    
    intensity.fwhm <- peak.df[,4] + (peak.df[,5] - peak.df[,4])/2 
    
    # Find the migration times closest to these intensities within each peak
    
    fwhm.vec <- vector()
    single.eie <- eie.df[,c(1,i+1)]
    
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
    
    peak.df <- subset(peak.df, peak.df$fwhm <= (fwhm.cutoff * is.df$peak.fwhm.tolerance.multiplier[i]))
    peak.df <- peak.df[,c(1:7)]
    
    # subset peak.df so that only the peaks (n = number.of.injections) with the greatest area are kept
    
    cut.off <- sort(peak.df[,7], decreasing = TRUE)[num.of.injections]
    
    peak.df <- subset(peak.df, peak.df[,7] >= cut.off)
    
    ## Peak space filtering ----
    
    # Determine the upper and lower migration time limits for space between peaks
    
    median.space <- peak.df[,2] %>%
      diff() %>%
      median()
    
    median.space.tol <- is.df[i,6] / 100
    
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
      
      peak.df <- peak.df[1:(num.of.injections-1),]
      
      # Typical case where an interior peak is missing (usually a blank)
      
      expected.peak.apex <- (peak.df[(bad.space[1] + 1),2] - peak.df[bad.space[1],2])/2 + peak.df[bad.space[1],2]
      
      # Case where final space is false - suspect that peak 1 was missed
      
      if(bad.space == (num.of.injections - 1) & false.peak.diff > median.space.upper.lim){
        expected.peak.apex <- peak.df[(num.of.injections -1 ),2] + median.space
      }
      
      # Case where final space is false - suspect that true final peak was missed
      
      if(bad.space == (num.of.injections - 1) & false.peak.diff < median.space.lower.lim){
        expected.peak.apex <- peak.df[1,2] - median.space
      }
      
      peak <- which.min(abs(peak.df.fill[,2] - expected.peak.apex))
      peak.df <- rbind(peak.df, peak.df.fill[peak,])
      
      peak.df <- peak.df[order(peak.df[,2]),]
    }
    
    ## Scenario 2 ---- 
    #Two or three bad spaces are detected 
    
    if(length(bad.space) == 2 | length(bad.space) == 3){
      
      peaks.to.check <- c(bad.space, (bad.space + 1))
      
      for (p in 1:length(peaks.to.check)){
        
        num.bad.space <- peak.df[,2] %>%
          .[-c(peaks.to.check[p])] %>%
          diff(.) %>%
          between (median.space.lower.lim, median.space.upper.lim) 
        num.bad.space <- length(which(num.bad.space == FALSE) + 1)
        
        if (num.bad.space == 1){
          peak.df <- peak.df[-c(peaks.to.check[p]),]
          break
        }
      }
      
      gap <- which(between(diff(peak.df[,2]), median.space.lower.lim, median.space.upper.lim) == FALSE)
      
      expected.peak.apex <- (peak.df[(gap[1] + 1),2] - peak.df[gap[1],2])/2 + peak.df[gap[1],2]
      
      peak <- which.min(abs(peak.df.fill[,2] - expected.peak.apex))
      peak.df <- rbind(peak.df, peak.df.fill[peak,])
      
      peak.df <- peak.df[order(peak.df[,2]),]
      
    }
    
    
    # Summarize peak.df data in is.peak.df
    
    if(i == 1){
      is.peaks.df = peak.df
    }else{
      is.peaks.df = cbind(is.peaks.df, peak.df)  
    }
  }
  
  # 7. Plotting ----
  
  for (i in 1:num.of.iss){
    
    ## Create annotation data frame ----
    
    peak.mt.df <- is.peaks.df[,seq(from = 2, to = ncol(is.peaks.df), by = 7)]
    
    ann.df <- data.frame("peak.number" = c(1:num.of.injections),
                         "peak.apex.seconds" = peak.mt.df[,i],
                         "peak.height.counts" = eie.df[which(eie.df$mt.seconds %in% (peak.mt.df[,i])),i+1])
    
    max.peak.height = max(ann.df$peak.height.counts)
    
    ## Create peak fill data frame ----
    
    pf.df <- data.frame("peak.number" = 1,
                        "mt.seconds" = eie.df[,1],
                        "intensity" = eie.df[,i+1])
    
    mt.vec <- vector()
    start.df <- is.peaks.df[,seq(from = 1, to = ncol(is.peaks.df), by = 7)]
    end.df <- is.peaks.df[,seq(from = 3, to = ncol(is.peaks.df), by = 7)]
    
    for (j in 1:num.of.injections){
      
      # Create a migration time vector to track where peaks elute
      
      mt.vec.temp <- eie.df$mt.seconds[between(eie.df$mt.seconds, start.df[j,i], end.df[j,i])]
      mt.vec <- append(mt.vec, mt.vec.temp)
      
      pf.df$peak.number <- ifelse(pf.df$mt.seconds >= start.df[j,i], j, pf.df$peak.number)
      pf.df$peak.number <- as.factor(pf.df$peak.number)
    }
    
    pf.df$intensity <- ifelse(pf.df$mt.seconds %in% mt.vec == TRUE, pf.df$intensity , 0)
    
    ## Add baseline intensity
    
    pf.df$baseline <- 0
    
    for (j in 1:num.of.injections){
      
        lower.intensity <- min(c(is.peaks.df[j, i * 7 - 3], is.peaks.df[j, i * 7 - 1]))
        pf.df$baseline <- ifelse(pf.df$mt.seconds >= start.df[j,i] & pf.df$mt.seconds <= end.df[j,i], lower.intensity, pf.df$baseline)
      }
    
    ## Plot ----
    
    name <- is.df$name[i]
    mz <- is.df$mz[i]
    
    ggplot(data = eie.df) +
      geom_line(aes(x = mt.seconds/60, y = eie.df[,i+1]), colour = "grey50") +
      theme_classic() +
      coord_cartesian(xlim = c(start.df[1,i]/60-1, end.df[num.of.injections,i]/60+1),
                      ylim = c(min(eie.df$mt.seconds[(start.df[1,m] + 60):(end.df[num.of.injections,m] + 60)]), 1.2 * max.peak.height)) +
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
           path = paste("Internal Standard Plots/", is.df$name[i], sep = ""))
    
  }
  
  # 8. Data for Export ----
  
  ## Generate peak area data frame ----
  
  peak.area.df <- cbind(file.name = c(data.files.name[d],2:num.of.injections),
                        is.peaks.df[,seq(from = 7, to = ncol(is.peaks.df), by = 7)])
  peak.area.df$file.name[2:num.of.injections] <- ""
  colnames(peak.area.df)[2:(num.of.iss+1)] <- is.df$name
  
  if(d == 1){
    peak.area.report = peak.area.df
  }else{
    peak.area.report = rbind(peak.area.report, peak.area.df)
  }
  
  ## Generate peak migration time data frame ----
  
  peak.mt.df <- cbind(file.name = c(data.files.name[d],2:num.of.injections),
                      is.peaks.df[,seq(from = 2, to = ncol(is.peaks.df), by = 7)] / 60)
  peak.mt.df$file.name[2:num.of.injections] <- ""
  colnames(peak.mt.df)[2:(num.of.iss+1)] <- is.df$name
  
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
          file = "csv Outputs/Internal Standard Peak Areas.csv",
          row.names = FALSE)

write.csv(peak.mt.report,
          file = "csv Outputs/Internal Standard Migration Times.csv",
          row.names = FALSE)

# close progress bar

close(pb)

