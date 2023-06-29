rm(list=ls())

# Install and load packages

if (!require("pacman", quietly = TRUE)) install.packages("pacman")
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

pacman::p_load("tidyverse", "stats", "DescTools", "xcms", "rlang",
               install = TRUE)

# 1. Load User Data and Parameters ----

# Import the user supplied tables of metabolites, internal standards, and general parameters

mass_df <- readxl::read_excel("Mass List and Parameters.xlsx", sheet = 1) %>%
  as.data.frame()

is_df <- readxl::read_excel("Mass List and Parameters.xlsx", sheet = 2) %>%
  as.data.frame()

parameters_df <- readxl::read_excel("Mass List and Parameters.xlsx", sheet = 3) %>%
  as.data.frame()

# Determine the number of metabolites and internal standards

num_of_metabolites <- nrow(mass_df)

num_of_is <- nrow(is_df)

# Get the number of injections

num_of_injections <- parameters_df$number.of.injections[1]

# Create an "Outputs" folder to store csv files

dir.create(path = "csv Outputs",
           showWarnings = FALSE)

# Create a "Plots" folder to store figures

dir.create(path = "Plots",
           showWarnings = FALSE)

# Generate plot sub-folders for plots of each internal standard and metabolite

name_vec <- c(is_df$name, mass_df$name)

for(n in 1:length(name_vec)){
  path = paste("Plots/",name_vec[n], sep = "")
  dir.create(path = path,
             showWarnings = FALSE)
}

rm(list = c("path", "n"))

# Create vector of data file directories and names 

data_files <- list.files(path = "mzML Files",
                         full.names = TRUE)

data_file_names <- list.files(path = "mzML Files")

# Initiate progress bar

pb <- winProgressBar(title = "Peak Seeker", 
                     label = paste("Number of Files Completed: 0 /", length(data_files)), 
                     min = 0,      
                     max = length(data_files), 
                     initial = 0,  
                     width = 300L) 

## Initiate data file for loop ----

for (d in 1:length(data_files)){
  
  print(paste(d, ". ", "Analyzing Data File: ", data_file_names[d], sep = ""))

  # 2. Read Data File ----

  # Make a copy of the data file as data will be written directly to this file during mass calibration
  
  file.copy(data_files[d], to = paste(data_files[d], "temp", sep = "_"))
  
  # Read copied data file
  
  print("Reading Data File")

  run_data <- readMSData(
    file = paste(data_files[d], "temp", sep = "_"),
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
  
  env_binding_unlock(run_data@assayData)
  
  # 3. Perform Mass Calibration ----
  
  print("Performing Mass Calibration")
  
  # Define mass window and minimum lock mass counts
  
  mass_window <- parameters_df$ref.mass.window.ppm[1]
  
  minimum_counts <- parameters_df$ref.mass.minimum.counts[1]
  
  # Define a function to calculate 1. experimental m/z of lock masses, 2. the corresponding mass error, and 3. the index of the lock mass
  # The lock mass value is currently the most intense point in the mass window plus an adjustment to better predict the lock mass value.
  
  calibration_parameters <- function(spectrum, lock_mass, minimum_counts, mass_window) {
    
    lock_mass_range <- c((lock_mass - lock_mass * mass_window/1000000), (lock_mass + lock_mass * mass_window/1000000))
    
    # Find the index of the most intense point in the lock mass window
    
    max_intensity_index <- spectrum %>%
      filter(mz >= lock_mass_range[1] & mz <= lock_mass_range[2]) %>%
      slice_max(intensity) %>%
      pull(index)
    
    if(length(max_intensity_index) == 1){
      if(spectrum$intensity[max_intensity_index] > minimum_counts){
        
        # Create a data frame with all four points
        line_points <- data.frame("mz" = c(spectrum$mz[(max_intensity_index - 2):(max_intensity_index - 1)], 
                                           spectrum$mz[(max_intensity_index + 1):(max_intensity_index + 2)]),
                                  "intensity" = c(spectrum$intensity[(max_intensity_index - 2):(max_intensity_index - 1)], 
                                                  spectrum$intensity[(max_intensity_index + 1):(max_intensity_index + 2)]))
        
        # Create a model for all four points
        model_left <- lm(formula = intensity ~ mz, data = line_points[c(1,2),])
        model_right <- lm(formula = intensity ~ mz, data = line_points[c(3,4),])
        
        # Get the coefficients for both models
        coefficients_left <- c(model_left$coefficients["mz"], model_left$coefficients["(Intercept)"])
        coefficients_right <- c(model_right$coefficients["mz"], model_right$coefficients["(Intercept)"])
        
        # Calculate the slope and intercept
        slope <- coefficients_left[1] - coefficients_right[1]
        intercept <- coefficients_right[2] - coefficients_left[2]
        
        # Solve for the experimental_mz
        experimental_mz <- solve(slope, intercept)
        
        experimental_mass_diff <- lock_mass - experimental_mz
        
        return(c(experimental_mz, experimental_mass_diff, max_intensity_index))
      }
    }
  }
  
  # Check if mass calibration should be applied
  
  if (parameters_df$apply.mass.calibration == "Yes"){
    
    # Loop through each spectrum to build a model and perform the mass calibration
    
    for (s in 1:end(rtime(run_data))[1]){
      
      spectrum_name <- ls(run_data@assayData)[s]
      
      spectrum <- data.frame(index = 1:length(run_data@assayData[[spectrum_name]]@mz),
                             mz = run_data@assayData[[spectrum_name]]@mz,
                             intensity = run_data@assayData[[spectrum_name]]@intensity)
      
      # Lower lock mass
      
      lock_mass <- parameters_df$ref.mass.one[1]
      
      cal_para_1 <- calibration_parameters(spectrum = spectrum,
                                           lock_mass = lock_mass,
                                           minimum_counts = minimum_counts,
                                           mass_window = mass_window)
      
      # Upper lock mass
      
      lock_mass <- parameters_df$ref.mass.two[1]
      
      cal_para_2 <- calibration_parameters(spectrum = spectrum,
                                           lock_mass = lock_mass,
                                           minimum_counts = minimum_counts,
                                           mass_window = mass_window)
      
      if (is.null(cal_para_1) | is.null(cal_para_2)){
        next
      }
      
      ## Develop correction model ----
      
      model_data <- data.frame("x" = c(cal_para_1[1], cal_para_2[1]),
                               "y" = c(cal_para_1[2], cal_para_2[2]))
      
      model <- lm(y ~ x, model_data)
      
      ## Apply correction model----
      
      correction_vector <- c(model[["coefficients"]][["x"]] * spectrum$mz) + model[["coefficients"]][["(Intercept)"]]
      
      correction_vector[1:cal_para_1[3]] <- cal_para_1[2]
      
      correction_vector[cal_para_2[3]:length(correction_vector)] <- cal_para_2[2]
      
      run_data@assayData[[spectrum_name]]@mz <- run_data@assayData[[spectrum_name]]@mz + correction_vector
      
      rm(list = c("cal_para_1", "cal_para_2"))
      
    }
    
    # clean-up global environment
    
    rm(list = c("model", "model_data", "correction_vector", "minimum_counts",
                "mass_window", "s", "spectrum", "lock_mass", "spectrum_name", 
                "calibration_parameters"))
    
    print("Mass Calibration Complete")
    
  }
  
  # 4. Extract Electropherograms ----
  
  print("Extracting Electropherograms")
  
  # Define mass error in ppm
  
  mass_error_vec <- c(is_df$extraction.window.ppm, mass_df$extraction.window.ppm)
  
  # Create a matrix of minimum and maximum m/z values for each internal standard and metabolite
  
  mz_vec <- c(is_df$mz, mass_df$mz)
  
  min <- mz_vec - mz_vec * mass_error_vec/1000000
  max <- mz_vec + mz_vec * mass_error_vec/1000000
  mzr <- matrix(c(min, max), ncol = 2)
  
  # Extract electropherograms
  
  electropherograms <- chromatogram(run_data,
                                    mz = mzr,
                                    rt = c(0,end(rtime(run_data))),
                                    aggregationFun = "mean",
                                    missing = 0,
                                    msLevel = 1)
  
  # Create a data frame of migration times and intensities with electropherograms data
  
  eie_df <- data.frame("mt.seconds" = electropherograms[1]@rtime)
  
  for (n in 1:length(name_vec)){
    temp_df <- data.frame(electropherograms[n]@intensity)
    colnames(temp_df) <- paste(name_vec[n], "intensity", sep = " ")
    eie_df <- cbind(eie_df,temp_df)
  }
  
  print("Extraction Complete")
  
  # 5. Smooth Intensity Vectors ----
  
  print("Smoothing Electropherograms")
  
  smoothing_kernel_vec <- c(is_df$smoothing.kernel, mass_df$smoothing.kernel)
  smoothing_strength_vec <- c(is_df$smoothing.strength, mass_df$smoothing.strength)
  
  for (n in 1:length(name_vec)){ 
    Smooth <- with(eie_df, 
                   ksmooth(x = mt.seconds, 
                           y = eie_df[,n + 1], 
                           kernel = smoothing_kernel_vec[n], 
                           bandwidth = smoothing_strength_vec[n]))
    eie_df[,n + 1] <- Smooth[["y"]]
  }
  
  # Clean-up global environment
  
  rm(list = c("electropherograms", "mzr", "Smooth", "temp_df", "max", "min", 
              "n", "mass_error_vec", "run_data"))
  
  print("Electropherograms Smoothing Complete")
  
  # 6. Migration Factor Calculation ----
  
  if(d == 1){
   
    # Summarize user supplied migration time data
    
    metabolites_mt_df <- data.frame(name = mass_df$name, mass_df[,c((ncol(mass_df) - num_of_injections + 1) : ncol(mass_df))])
    is_mt_df <- subset(is_df, is_df$class == "Internal Standard")
    is_mt_df <- data.frame(name = is_mt_df$name, is_mt_df[,c((ncol(is_mt_df) - num_of_injections + 1) : ncol(is_mt_df))])
    
    # Determine IS on the left
    
    is_left_vec <- c(1:num_of_metabolites)
    
    for (m in 1:nrow(metabolites_mt_df)){
      
      is_temp <- (metabolites_mt_df[m,2] - (is_mt_df[,2]))
      is_temp <- set_names(is_temp, is_mt_df$name)
      is_temp <- is_temp[which(sign(is_temp) == 1 | sign(is_temp) == 0)] %>%
        which.min() %>%
        names()
      
      if(length(is_temp) == 0){
        is_temp <- "none"
      }
      
      is_left_vec[m] <- is_temp
    }
    
    # Determine IS on the Right
    
    is_right_vec <- c(1:num_of_metabolites)
    
    for (m in 1:nrow(metabolites_mt_df)){
      
      is_temp <- (metabolites_mt_df[m,2] - (is_mt_df[,2]))
      is_temp <- set_names(is_temp, is_mt_df$name)
      is_temp <- is_temp[which(sign(is_temp) == -1| sign(is_temp) == 0)] %>%
        which.max() %>%
        names()
      
      if(length(is_temp) == 0){
        is_temp <- "none"
      }
      
      is_right_vec[m] <- is_temp
    }
    
    # Compute migration factors
    
    mf_df <- matrix(NA, ncol = (num_of_injections + 1), nrow = nrow(metabolites_mt_df)) %>%
      as.data.frame()
    mf_df[,1] <- metabolites_mt_df$name
    
    summary_vec <- c(1:nrow(metabolites_mt_df))
    
    for (m in 1:nrow(metabolites_mt_df)){
      for (i in 2:(num_of_injections + 1)){
        
        left <- is_left_vec[m]
        right <- is_right_vec[m]
        
        # If a left or right internal standard does not exist, calculate the relative 
        # migration time to the nearest internal standard instead
        
        if(left == "none"){
          
          mf_df[m,i] <- metabolites_mt_df[m,i] / is_mt_df[which(is_mt_df$name == right),i]
          summary_vec[m] <- right
          next
          
        }
        
        if(right == "none"){
          
          mf_df[m,i] <- metabolites_mt_df[m,i] / is_mt_df[which(is_mt_df$name == left),i]
          summary_vec[m] <- left
          next
          
        }
        
        # If the metabolite migrates near either neighboring internal standards, 
        # calculate the relative migration time to the nearest internal standard instead
        
        if(metabolites_mt_df[m,2] / is_mt_df[which(is_mt_df$name == left),2] < 1.01){
          
          mf_df[m,i] <-metabolites_mt_df[m,i] / is_mt_df[which(is_mt_df$name == left),i]
          summary_vec[m] <- left
          next
          
        }
        
        if(metabolites_mt_df[m,2] / is_mt_df[which(is_mt_df$name == right),2] > 0.99){
          
          mf_df[m,i] <- metabolites_mt_df[m,i] / is_mt_df[which(is_mt_df$name == right),i]
          summary_vec[m] <- right
          next
          
        }
        
        # Calculate the retention factor
        
        mf_df[m,i] <- (metabolites_mt_df[m,i] - is_mt_df[which(is_mt_df$name == left),i]) / (is_mt_df[which(is_mt_df$name == right),i] - is_mt_df[which(is_mt_df$name == left),i])
        summary_vec[m] <- "rf"
        
      }
    }
  }
  
  print("Migration Factor Calculation Complete")
  
  # 7. Internal Standard Peak Detection, Integration, and Filtering ----
  
  print("Performing Peak Picking and Filtering for Internal Standards")
  
  ## Peak detection ---- 
  
  n <- parameters_df$required.points.for.peak.picking[1]
  
  for (s in 1:num_of_is){
    
    rle_output <- eie_df[,s+1] %>%
      diff() %>%
      sign() %>%
      rle()
    
    consecutive_runs <- which(rle_output$lengths > n & rle_output$values == 1)
    consecutive_runs <- subset(consecutive_runs, (consecutive_runs + 1) %in% (which(rle_output$lengths > n)) == TRUE)
    
    run_lengths <- cumsum(rle_output$lengths) + 1
    
    start <- eie_df$mt.seconds[run_lengths[consecutive_runs - 1]]
    apex <- eie_df$mt.seconds[run_lengths[consecutive_runs]]
    end <- eie_df$mt.seconds[run_lengths[consecutive_runs + 1]]
    
    # For FWHM calculations I will also add intensity values here as well
    
    start_intensity <- eie_df[run_lengths[consecutive_runs - 1], (s+1)]
    apex_intensity <- eie_df[run_lengths[consecutive_runs], (s+1)]
    end_intensity <- eie_df[run_lengths[consecutive_runs + 1], (s+1)]
    
    # Account for peaks that start immediately during the analysis
    
    if(length(start) != length(apex)){
      start <- append(start, 0, 0)
      start_intensity <- append(start_intensity, 0, 0)
    }
    
    # Create a data frame containing the start, apex, and end migration times of each peak
    
    peak_df <- data.frame(start,
                          apex,
                          end,
                          start_intensity,
                          apex_intensity,
                          end_intensity)
    
    ## Migration time filtering ----
    
    # Filter peaks that are outside migration time limits
    
    peak_df <- subset(peak_df, peak_df$start >= is_df$min.rt.min[s] * 60 & peak_df$end <= is_df$max.rt.min[s] * 60)
    
    ## Integrate peaks ----
    
    peak_area_vector = c(1:nrow(peak_df))
    
    for (p in 1:nrow(peak_df)){
      
      peak_area_vector[p] <- AUC(eie_df$mt.seconds,
                                 eie_df[,s+1],
                                 method = "trapezoid",
                                 from = peak_df[p,1],
                                 to = peak_df[p,3],
                                 absolutearea = FALSE,
                                 na.rm = FALSE)
      
      # Perform baseline correction
      
      peak_area_vector[p] <- peak_area_vector[p] - (peak_df[p,3] - peak_df[p,1]) * min(peak_df[p,4], peak_df[p,6])
    }
    
    peak_df <- cbind(peak_df, peak_area_vector)
    
    # rename peak_df columns
    
    colnames.vector = c(paste(name_vec[s], "start.seconds", sep = "."),
                        paste(name_vec[s], "apex.seconds", sep = "."),
                        paste(name_vec[s], "end.seconds", sep = "."),
                        paste(name_vec[s], "start_intensity", sep = "."),
                        paste(name_vec[s], "apex_intensity", sep = "."),
                        paste(name_vec[s], "end_intensity", sep = "."),
                        paste(name_vec[s], "peak.area", sep = "."))
    
    colnames(peak_df) <- colnames.vector
    
    # Retain peak_df for future filtering steps
    
    peak_df_fill <- peak_df
    
    ## FWHM filtering ----
    
    # Find the peak intensity at half the peak height
    
    intensity_fwhm <- peak_df[,4] + (peak_df[,5] - peak_df[,4])/2 
    
    # Find the migration times closest to these intensities within each peak
    
    fwhm_vec <- vector()
    single_eie <- eie_df[,c(1,s+1)]
    
    for(p in 1:nrow(peak_df)){
      
      single_eie_temp <- subset(single_eie, single_eie$mt.seconds >= peak_df[p,1] & single_eie$mt.seconds <= peak_df[p,2])
      fwhm_mt_left <- single_eie_temp$mt.seconds[which.min(abs(single_eie_temp[,2] - intensity_fwhm[p]))]
      
      single_eie_temp <- subset(single_eie, single_eie$mt.seconds <= peak_df[p,3] & single_eie$mt.seconds >= peak_df[p,2])
      fwhm_mt_right <- single_eie_temp$mt.seconds[which.min(abs(single_eie_temp[,2] - intensity_fwhm[p]))]
      
      fwhm_vec <- append(fwhm_vec, fwhm_mt_right - fwhm_mt_left)
      
    }
    
    peak_df$fwhm <- fwhm_vec
    
    # Determine the fwhm of the peaks (n = num_of_injections) with the greatest area
    
    df_temp <- peak_df[order(-peak_df[,8]),]
    df_temp <- df_temp[1:num_of_injections,]
    
    fwhm_cutoff <- median(df_temp$fwhm)
    
    peak_df <- subset(peak_df, peak_df$fwhm <= (fwhm_cutoff * is_df$peak.fwhm.tolerance.multiplier[s]))
    peak_df <- peak_df[,c(1:7)]
    
    # subset peak_df so that only the peaks (n = number.of.injections) with the greatest area are kept
    
    cut_off <- sort(peak_df[,7], decreasing = TRUE)[num_of_injections]
    
    peak_df <- subset(peak_df, peak_df[,7] >= cut_off)
    
    ## Peak space filtering ----
    
    # Determine the upper and lower migration time limits for space between peaks
    
    median_space <- peak_df[,2] %>%
      diff() %>%
      median()
    
    median_space_tol <- is_df$peak.space.tolerance.percent[s] / 100
    
    median_space_lower_lim <- median_space - median_space * median_space_tol
    median_space_upper_lim <- median_space + median_space * median_space_tol
    
    # Check if peaks migrate within the tolerance limits
    
    peak_space_tol_check <- between(diff(peak_df[,2]), median_space_lower_lim, median_space_upper_lim)
    
    if(all(peak_space_tol_check) != TRUE){
      bad_space <- which(peak_space_tol_check == FALSE)
    }else{
      bad_space <- NA
    }
    
    ## Scenario 1 ----
    # Only one bad space is detected
    
    if(length(bad_space) == 1 & is.na(bad_space[1]) == FALSE){
      
      false_peak_diff <- peak_df[num_of_injections, 2] - peak_df[(num_of_injections - 1), 2]
      
      # This algorithm always assumes the final peak is false - likely due to carryover
      # It is possible the first peak is false but this seems less likely
      # Final peak only removed if case 1 does not produce a duplicate
      
      # Case 1- An interior peak is missing (usually a blank)
      
      if(bad_space != (num_of_injections - 1)){
        
        # Since the we know the bad space is not at the end, use the space after the bad space to find the expected apex
        
        expected_peak_apex <- peak_df[bad_space[1],2] + peak_df[(bad_space[1] + 2),2] - peak_df[(bad_space[1] + 1),2]
        
      }
      
      # Case 2 -  Final space is false and less than median - suspect that true final peak was missed
      
      if(bad_space == (num_of_injections - 1) & false_peak_diff < median_space_upper_lim){
        expected_peak_apex <- peak_df[(num_of_injections - 1 ),2] + median_space
      }
      
      # Case 3 - Final space is false and greater than median - suspect that peak 1 was missed
      
      if(bad_space == (num_of_injections - 1) & false_peak_diff > median_space_lower_lim){
        expected_peak_apex <- peak_df[1,2] - median_space
      }
      
      peak <- which.min(abs(peak_df_fill[,2] - expected_peak_apex))
      
      # If the nearest peak is too far from the expected migration time use a place holder
      
      if (abs(peak_df_fill[peak,2] - expected_peak_apex) > (median_space / 2)){
        
        nearest_mt <- (eie_df$mt.seconds - expected_peak_apex) %>%
          abs() %>%
          which.min(.)
        
        peaks <- data.frame(eie_df[nearest_mt,1],
                            eie_df[nearest_mt,1],
                            eie_df[nearest_mt,1],
                            eie_df[nearest_mt,s + 1],
                            eie_df[nearest_mt,s + 1],
                            eie_df[nearest_mt,s + 1],
                            0)
        
        colnames(peaks) <- colnames(peak_df)
        
        peak_df <- rbind(peak_df, peaks)
        
        peak_df <- peak_df[order(peak_df[,2]),]
        
      }else{
        
        peak_df <- rbind(peak_df, peak_df_fill[peak,])
        
        peak_df <- peak_df[order(peak_df[,2]),]
        
      }
      
      # Remove bad peak
      
      if (any(duplicated(peak_df[,2]))){
        
        peak_df <- peak_df[-c(which(duplicated(peak_df[,2]))),]
        
      }else{
        
        peak_df <- peak_df[1:(num_of_injections),]
        
      }
    }
    
    ## Scenario 2 ---- 
    #Two or three bad spaces are detected 
    
    if(length(bad_space) == 2 | length(bad_space) == 3){
      
      peaks_to_check <- c(bad_space, (bad_space + 1)) %>%
        unique() %>%
        sort()
      
      # Outside peaks must be checked last
      
      peaks_to_check <- c(peaks_to_check[-1], peaks_to_check[1])
      
      # Loop through each suspect bad peak and remove them iteratively until the number of bad spaces reaches 1
      
      for (p in 1:length(peaks_to_check)){
        
        num_bad_space <- peak_df[,2] %>%
          .[-c(peaks_to_check[p])] %>%
          diff(.) %>%
          between (median_space_lower_lim, median_space_upper_lim) 
        num_bad_space <- length(which(num_bad_space == FALSE))
        
        if (num_bad_space == 1){
          peak_df <- peak_df[-c(peaks_to_check[p]),]
          break
        }
      }
      
      # To avoid errors where removing a peak results in 0 bad spaces only fill
      # in gap if nrow(peak_df) == number of injections - 1
      
      if(nrow(peak_df) == (num_of_injections - 1)){
        
        # Find the peak gap and calculate the expected migration time for the missing peak
        
        gap <- which(between(diff(peak_df[,2]), median_space_lower_lim, median_space_upper_lim) == FALSE)
        
        expected_peak_apex <- (peak_df[(gap[1] + 1),2] - peak_df[gap[1],2])/2 + peak_df[gap[1],2]
        
        # Find the nearest peak in the peak_df_fill data frame to the expected migration time
        # Avoid duplicate peaks by not using exisitng peaks in peak_df
        
        peak_df_fill <- subset(peak_df_fill, !(peak_df_fill[,2] %in% peak_df[,2]))
        
        peak <- which.min(abs(peak_df_fill[,2] - expected_peak_apex))
        peak_df <- rbind(peak_df, peak_df_fill[peak,])
        
        peak_df <- peak_df[order(peak_df[,2]),]
        
      }
    }
    
    # Summarize peak_df data in is.peak_df
    
    if(s == 1){
      is_peaks_df = peak_df
    }else{
      is_peaks_df = cbind(is_peaks_df, peak_df)  
    }
  }
  
  # Make a data frame containing the apex migration times of the internal standards
  # to be used to filter metabolite peaks
  
  is_mt_df <- is_peaks_df[,seq(from = 2,
                               to = ncol(is_peaks_df),
                               by = 7)]
  
  print("Peak Picking and Filtering for Internal Standards Complete")
  
  # 8. Metabolite  Peak Detection, Integration, and Filtering ----
  
  print("Performing Peak Picking and Filtering for Analytes")
  
  ## Peak detection ---- 
  
  for (m in (num_of_is + 1):length(name_vec)){
    
    peak_df <- data.frame()
      
    # Determine the start, apex, and end of peaks. Use the user defined value "n" to detect peaks.
    # If n results in fewer peaks then injection, decrease n by 1 and repeat
    
    n <- parameters_df$required.points.for.peak.picking[1]
    
    while (nrow(peak_df) < num_of_injections & n >= 0){
      
      rle_output <- eie_df[,m + 1] %>%
        diff() %>%
        sign() %>%
        rle()
      
      consecutive_runs <- which(rle_output$lengths > n & rle_output$values == 1)
      consecutive_runs <- subset(consecutive_runs, (consecutive_runs + 1) %in% (which(rle_output$lengths > n)) == TRUE)
      
      run_lengths <- cumsum(rle_output$lengths) + 1
      
      start <- eie_df$mt.seconds[run_lengths[consecutive_runs - 1]]
      apex <- eie_df$mt.seconds[run_lengths[consecutive_runs]]
      end <- eie_df$mt.seconds[run_lengths[consecutive_runs + 1]]
      
      # I will also add intensity values here as well
      
      start_intensity <- eie_df[run_lengths[consecutive_runs - 1], (m+1)]
      apex_intensity <- eie_df[run_lengths[consecutive_runs], (m+1)]
      end_intensity <- eie_df[run_lengths[consecutive_runs + 1], (m+1)]
      
      # Account for peaks that start immediately during the analysis
      
      if(length(start) != length(apex)){
        start <- append(start, 0, 0)
        start_intensity <- append(start_intensity, 0, 0)
      }
      
      # Create a data frame containing the start, apex, and end migration times of each 
      # peak in addition to required intensities for FWHM calculations
      
      peak_df <- data.frame(start,
                            apex,
                            end,
                            start_intensity,
                            apex_intensity,
                            end_intensity)
      
      n <- n - 1
      
    }
    
    ## Filter peaks by peak width ----
    
    # Define a minimum peak width cut off in seconds. Remove peaks with a width <= cutoff
    # If the cutoff results in fewer peaks than injections, decrease cutoff by 1 and repeat
    
    min_width_cut_off <- mass_df$minimim.peak.width.seconds[m - num_of_is]
    
    peak_df_trim <- subset(peak_df, (peak_df$end - peak_df$start) >= min_width_cut_off)
    
    while (nrow(peak_df_trim) < num_of_injections){
      
      min_width_cut_off <- min_width_cut_off - 1
      
      peak_df_trim <- subset(peak_df, (peak_df$end - peak_df$start) >= min_width_cut_off)
      
    }
    
    peak_df <- peak_df_trim
    
    ## Integrate peaks ----
    
    peak_area_vector = c(1:nrow(peak_df))
    
    for (p in 1:nrow(peak_df)){
      
      peak_area_vector[p] <- AUC(eie_df$mt.seconds,
                                 eie_df[,m+1],
                                 method = "trapezoid",
                                 from = peak_df[p,1],
                                 to = peak_df[p,3],
                                 absolutearea = FALSE,
                                 na.rm = FALSE)
      
      peak_area_vector[p] <- peak_area_vector[p] - (peak_df[p,3] - peak_df[p,1]) * min(peak_df[p,4], peak_df[p,6])
    }
    
    peak_df <- cbind(peak_df, peak_area_vector)
    
    # rename peak_df columns
    
    colnames.vector = c(paste(name_vec[m], "start.seconds", sep = "."),
                        paste(name_vec[m], "apex.seconds", sep = "."),
                        paste(name_vec[m], "end.seconds", sep = "."),
                        paste(name_vec[m], "start_intensity", sep = "."),
                        paste(name_vec[m], "apex_intensity", sep = "."),
                        paste(name_vec[m], "end_intensity", sep = "."),
                        paste(name_vec[m], "peak.area", sep = "."))
    
    colnames(peak_df) <- colnames.vector
    
    ## Filter peaks ----
    
    # Build a data frame to store comments for each metabolite peak
    
    comment_df <- matrix(nrow = num_of_injections, ncol = num_of_metabolites, "") %>%
      as.data.frame
    
    colnames(comment_df) <- mass_df$name
    
    ### Filter peaks based on smallest rmt difference ----
    
    # Determine the expected migration times of the metabolites
    
    n <- m - num_of_is
    left_is <- paste(is_left_vec[n], ".apex.seconds", sep ="")
    right_is <- paste(is_right_vec[n], ".apex.seconds", sep ="")
    rmt_is <- paste(summary_vec[n], ".apex.seconds", sep ="")

    if(summary_vec[n] == "rf"){
      
      mf_vec <- mf_df[n,2:ncol(mf_df)] %>%
        unlist() %>%
        unname()

      expected_mt <- mf_vec * (is_mt_df[,which(colnames(is_mt_df) == right_is)] - is_mt_df[,which(colnames(is_mt_df) == left_is)]) + is_mt_df[,which(colnames(is_mt_df) == left_is)]
      
    }else{
    
      mf_vec <- mf_df[n,2:ncol(mf_df)] %>%
        unlist() %>%
        unname()
      
      expected_mt <- mf_vec * is_mt_df[,which(colnames(is_mt_df) == rmt_is)] 
      
    }
    
    # Filter peak_df for peaks within rmt tolerance
    
    rmt_tolerance <- mass_df$rmt.tolerance.percent[m - num_of_is] / 100
    
    for (i in 1:num_of_injections){
      
      peaks <- peak_df %>%
        filter(., peak_df[,2] <= (1 + rmt_tolerance) * expected_mt[i] & 
                 peak_df[,2] >= (1 - rmt_tolerance) * expected_mt[i])
      
      # If more than one peak is found choose the nearest one
      
      if(nrow(peaks) > 1){
        peaks <- (peak_df[,2] - expected_mt[i]) %>%
          abs() %>%
          which.min(.)
        peaks <- peak_df[peaks,]
      }
      
      # If no peak is found, generate a place holder
      
      if(nrow(peaks) == 0){
        
        nearest_mt <- (eie_df$mt.seconds - expected_mt[i]) %>%
          abs() %>%
          which.min(.)
        
        peaks <- data.frame(eie_df[nearest_mt,1],
                            eie_df[nearest_mt,1],
                            eie_df[nearest_mt,1],
                            eie_df[nearest_mt,m + 1],
                            eie_df[nearest_mt,m + 1],
                            eie_df[nearest_mt,m + 1],
                            0)
        
        colnames(peaks) <- colnames(peak_df)
      }

      if(i == 1){
        filtered_peaks_df <- peaks
      }else{
        filtered_peaks_df <- rbind(filtered_peaks_df, peaks)
      }
    }
    
    ### Filter peaks outside of run time limits ----
    
    # Set a place holder for peaks where the expected migration time > total run time
    
    total_run_time <- eie_df$mt.seconds[nrow(eie_df)]
    
    late_peaks <- (expected_mt > total_run_time) %>%
      which()
    
    # Find migration times to use as placeholders that do not belong to other identified peaks
    
    mt <- tail(eie_df$mt.seconds, n = 15)
    
    mt <- mt[!(mt %in% filtered_peaks_df[,2])]
    
    for (i in late_peaks){
      mt_temp <- mt[i]
      filtered_peaks_df[i,] <- c(mt_temp,
                                 mt_temp,
                                 mt_temp,
                                 eie_df[which(eie_df$mt.seconds == mt_temp) ,m + 1],
                                 eie_df[which(eie_df$mt.seconds == mt_temp) ,m + 1],
                                 eie_df[which(eie_df$mt.seconds == mt_temp) ,m + 1],
                                 0)
    }
    
    
    ### Filter duplicated peaks ----
    
    # If the same peak is assigned to multiple injection numbers, reapply rmt filter with more austere rmt tolerances
    # New rmt tolerance will be the original / count, which starts at 2 and increases by 1 each iteration
    
    count = 2
    
    while (any(duplicated(filtered_peaks_df[,2])) & count < 100){
      
      strict_rmt_tolerance <- rmt_tolerance/count
      
      # find rows with duplicated values 
      
      duplicate_location <- filtered_peaks_df[,2] %>%
        duplicated() %>%
        which()
      
      duplicate_rows <- which(filtered_peaks_df[,2] %in% filtered_peaks_df[duplicate_location,2])
      
      # reapply filtering for these peaks with the more strict rmt tolerance
      
      for (r in duplicate_rows){
        
        peaks <- peak_df %>%
          filter(., peak_df[,2] <= (1 + strict_rmt_tolerance) * expected_mt[r] & 
                   peak_df[,2] >= (1 - strict_rmt_tolerance) * expected_mt[r])
        
        # If more than one peak is found choose the nearest one
        
        if(nrow(peaks) > 1){
          peaks <- (peak_df[,2] - expected_mt[r]) %>%
            abs() %>%
            which.min(.)
          peaks <- peak_df[peaks,]
        }
        
        # If no peak is found, generate a place holder
        
        if(nrow(peaks) == 0){
          
          nearest_mt <- (eie_df$mt.seconds - expected_mt[r]) %>%
            abs() %>%
            which.min(.)
          
          peaks <- data.frame(eie_df[nearest_mt,1],
                              eie_df[nearest_mt,1],
                              eie_df[nearest_mt,1],
                              eie_df[nearest_mt,m + 1],
                              eie_df[nearest_mt,m + 1],
                              eie_df[nearest_mt,m + 1],
                              0)
          
          colnames(peaks) <- colnames(peak_df)
        }
        
        filtered_peaks_df[r,] <- peaks
        
      }
      
      count = count + 1
    }
    
    # If the duplicate peak filter fails (count = 100) then generate a place holder peak_df
    # and generate a warning that the filter failed. This filter fails when two or more expected migration 
    # times are too close to each other.
    
    if (count == 100){
      
      mt.temp <- eie_df$mt.seconds[seq(1, num_of_injections * 10, 10)]
      
      filtered_peaks_df[,c(1:3)] <- mt.temp
      filtered_peaks_df[,c(4:7)] <- 0
      
      comment_df[1,name_vec[m]] <- "Peak Filtering Failed to Remove Duplicate Peaks"
      
      print(paste("WARNING: Peak Filtering Failed for", name_vec[m]))
    }
    
    ### Filter using peak spaces ----
    
    # Get peak space tolerance
    
    space_tol <- mass_df$peak.space.tolerance.percent[m - num_of_is] / 100
    
    # Make a vector containing all the expected space lengths
    
    space_vec <- expected_mt %>%
      diff()
    
    # Define upper and lower peak space limits
    
    space_lower_lim <- space_vec - space_vec * space_tol
    space_upper_lim <- space_vec + space_vec * space_tol
    
    # Check if peaks migrate within the tolerance limits
    
    peak_space_tol_check <- between(diff(filtered_peaks_df[,2]), space_lower_lim, space_upper_lim)
    
    if(all(peak_space_tol_check) != TRUE){
      bad_space <- which(peak_space_tol_check == FALSE)
    }else{
      bad_space <- NA
    }
    
    # identify which peaks are potentially incorrectly assigned (bad peaks)
    # these are peaks before and after each bad space
    
    bad_peaks <- c(bad_space, bad_space + 1) %>%
      unique() %>%
      sort()
    
    # Define a count that will be used to modify the peak space tolerance
    
    count = 4
    
    # Keep unaltered filtered peaks data frame for the event that the peak space algorithm fails
    
    filtered_peaks_df_retain <- filtered_peaks_df
    
    # set count limit to determine when the algorithm fails
    
    count_limit = 100
    
    # Identify bad peaks, and replace them with peaks meeting peak space criteria
    # If the number of bad peaks is equal to the number of injections, do not apply this filter
    
    while(length(bad_peaks) > 0 & length(bad_peaks) < (num_of_injections - 1) & count < count_limit){
      
      # define remaining peaks which are correctly assigned (good peaks)
      
      good_peaks <- c(1:num_of_injections) %>%
        setdiff(., c(bad_peaks))
      
      # Use the median of the expected peak space times to find peaks
      
      peak_tolerance <- expected_mt %>%
        diff() %>%
        median () / count
      
      for (b in 1:length(bad_peaks)){
        
        # find the nearest good peak neighbor for each bad peak
        
        nearest_good_peak <- (good_peaks - bad_peaks[b]) %>%
          abs() %>%
          which.min()
        
        # calculate the expected migration time 
        
        expected_mt <- filtered_peaks_df[good_peaks[nearest_good_peak], 2] - 
          (good_peaks[nearest_good_peak] - bad_peaks[b]) * median(space_vec)
        
        # find peaks nearest to the expected migration time within the tolerance
        
        peaks <- peak_df %>%
          filter(., peak_df[,2] <= expected_mt + peak_tolerance & peak_df[,2] >= expected_mt - peak_tolerance)
        
        # if more than one peak is found, select the closest one
        
        if(nrow(peaks) > 1){
          peaks <- (peak_df[,2] - expected_mt) %>%
            abs() %>%
            which.min(.)
          peaks <- peak_df[peaks,]
        }
        
        # if no peaks are found, define a place holder
        
        if(nrow(peaks) == 0){
          
          nearest_mt <- (eie_df$mt.seconds - expected_mt) %>%
            abs() %>%
            which.min(.)
          
          peaks <- data.frame(eie_df[nearest_mt,1],
                              eie_df[nearest_mt,1],
                              eie_df[nearest_mt,1],
                              eie_df[nearest_mt,m + 1],
                              eie_df[nearest_mt,m + 1],
                              eie_df[nearest_mt,m + 1],
                              0)
          
          colnames(peaks) <- colnames(peak_df)
          
          filtered_peaks_df[bad_peaks[b],] <- peaks
        }
        
        # if only one peak is found
        
        filtered_peaks_df[bad_peaks[b],] <- peaks
        
      }
      
      count = count + 1
      
      duplicate_location <- filtered_peaks_df[,2] %>%
        duplicated() %>%
        which()
      
      # bad peaks correspond to any rows that are not unique
      
      bad_peaks <- which(filtered_peaks_df[,2] %in% filtered_peaks_df[duplicate_location,2])
      
    }
    
    # if algorithm failed, revert back to filtered_peaks_df
    
    if (count == count_limit){
      
      filtered_peaks_df <- filtered_peaks_df_retain
      
    }
    
    # Summarize filtered.peak_df data in metabolite_peak_df
    
    if(m == (num_of_is + 1)){
      metabolite_peaks_df = filtered_peaks_df
    }else{
      metabolite_peaks_df = cbind(metabolite_peaks_df, filtered_peaks_df)  
    }
    
 }
  
  ### Filter peaks below LOD ----
  
  # Loop through each metabolite and see if its area is below the LOD threshold
  
  for (m in 1:num_of_metabolites){
    
    # Determine the noise of the electropherogram
    # Fine the noise levels in 60 seconds intervals
    
    region_start <- seq(1, nrow(eie_df), 60)
    region_end <- seq(60, nrow(eie_df), 60)
    length(region_start) <- length(region_end)
    
    # Generate a vector to store noise data
    
    noise_vec <- rep(NA, length(region_start))
    
    # Define a function to calculate noise
    
    noise_calculation <- function(temp_noise) {
      mean(temp_noise) + mass_df$snr.threshold[m] * sd(temp_noise)
    }
    
    # Calculate the noise in each region
    
    for (r in 1:length(region_start)){
      temp_noise <- eie_df[region_start[r]:region_end[r], m + num_of_is + 1]
      noise_vec[r] <- noise_calculation(temp_noise)
    }

    # Define the noise as the 20th percentile noise region
    
    noise <- noise_vec %>%
      sort()
    
    noise <- noise[as.integer(length(noise)/5)]
    
    peak_area_df <- metabolite_peaks_df[,seq(7, ncol(metabolite_peaks_df), 7)]
    
    comment_df[,m] <- ifelse(peak_area_df[,m] < noise, "<LOD", comment_df[,m])
    
    ### Annotate injections that are not detected
    
    comment_df[,m] <- ifelse(peak_area_df[,m] == 0, "NPD", comment_df[,m])
    
  }
  
  ### Filter interfered peaks ----
  
  # Build an interference data frame since some are metabolites and some are internal standards
  
  interference_df <- cbind(is_peaks_df[,seq(2, ncol(is_peaks_df), 7)],
                           metabolite_peaks_df[,seq(2, ncol(metabolite_peaks_df), 7)])
  
  for (m in 1:num_of_metabolites){
    
    # Skip metabolites with no reported interference
    
    if(is.na(mass_df$interference[m])){
      next
    }
    
    # Get the names of the interferences from mass_df
    
    interferences <- strsplit(mass_df$interference[m], ", ") %>%
      unlist()
    
    for (k in 1:length(interferences)){
      
      interference <- paste(interferences[k], ".apex.seconds", sep = "")
      
      # See if there is any overlap between the metabolite peak and its interference 
      
      metabolite_name <- paste(mass_df$name[m], ".apex.seconds", sep = "")
      
      for (i in 1:num_of_injections){
        for (j in 1:num_of_injections){
          
          diff_temp <- (metabolite_peaks_df[i,metabolite_name] - interference_df[j,interference]) %>%
            abs() 
          
          comment_df[i,mass_df$name[m]] <- ifelse(diff_temp < mass_df$interference.comigration.threshold.seconds[m], "Interfered", comment_df[i,mass_df$name[m]])
        }
      }
    }
  }
  
  # Combine internal standard and metabolite data frames for plotting
  
  peaks_df <- cbind(is_peaks_df, metabolite_peaks_df)
  
  # update comment data frame account for internal standards
  
  is_comment_df <- matrix(nrow = num_of_injections, ncol = nrow(is_df), "") %>%
    as.data.frame()
  
  colnames(is_comment_df) <- is_df$name
  
  comment_df <- cbind(is_comment_df, comment_df)
  
  print("Peak Picking and Filtering for Analytes Complete")
  
  # 9. Plotting ----
  
  print("Plotting Electropherograms")
  
  for (n in 1:length(name_vec)){
    
    ## Create annotation data frame ----
    
    peak_mt_df <- peaks_df[,seq(from = 2, to = ncol(peaks_df), by = 7)]
      
      ann_df <- data.frame("peak.number" = c(1:num_of_injections),
                           "comment" = comment_df[,n],
                           "peak.apex.seconds" = peak_mt_df[,n],
                           "peak.height.counts" = eie_df[which(eie_df$mt.seconds %in% (peak_mt_df[,n])),n+1])
    
    max_peak_height = max(ann_df$peak.height.counts)
    
    ## Create peak fill data frame ----
    
    pf_df <- data.frame("peak.number" = 1,
                        "mt.seconds" = eie_df[,1],
                        "intensity" = eie_df[,n+1])
    
    mt_vec <- vector()
    start_df <- peaks_df[,seq(from = 1, to = ncol(peaks_df), by = 7)]
    end_df <- peaks_df[,seq(from = 3, to = ncol(peaks_df), by = 7)]
    
    for (i in 1:num_of_injections){
      
      # Create a migration time vector to track where peaks elute
      
      if(comment_df[i,n] == ""){
        mt_vec_temp <- eie_df$mt.seconds[between(eie_df$mt.seconds, start_df[i,n], end_df[i,n])]
        mt_vec <- append(mt_vec, mt_vec_temp)
        
        # Update peak.number in pf_df
        
        pf_df$peak.number <- ifelse(pf_df$mt.seconds >= start_df[i,n], i, pf_df$peak.number)
        pf_df$peak.number <- as.factor(pf_df$peak.number)
      }else{
        next
      }
    }
    
    pf_df$intensity <- ifelse(pf_df$mt.seconds %in% mt_vec == TRUE, pf_df$intensity , 0)
    
    ## Add baseline intensity
    
    pf_df$baseline <- 0
    
    for (i in 1:num_of_injections){
      if(comment_df[i,n] == ""){
        lower_intensity <- min(c(peaks_df[i, n * 7 - 3], peaks_df[i, n * 7 - 1]))
        pf_df$baseline <- ifelse(pf_df$mt.seconds >= start_df[i,n] & pf_df$mt.seconds <= end_df[i,n], lower_intensity, pf_df$baseline)
      }else{
        next
      }
    }
    
    ## Plot ----
    
    name <- name_vec[n]
    mz <- mz_vec[n]
    start_mt <- start_df[1,n]
    end_mt <- end_df[num_of_injections,n]
    extra_space <- ifelse(start_mt > 70, 1, 0)
    ymin = min(eie_df[(which(eie_df$mt.seconds == start_mt) - extra_space * 60) : (which(eie_df$mt.seconds == end_df[1,n]) + extra_space * 60),n + 1])
    
    ggplot(data = eie_df) +
      geom_line(aes(x = mt.seconds/60, y = eie_df[,n+1]), colour = "grey50") +
      theme_classic() +
      coord_cartesian(xlim = c(start_mt/60 - extra_space, end_mt/60 + extra_space),
                      ylim = c(ymin / 3,
                               1.2 * max_peak_height)) +
      scale_y_continuous(name = "Ion Counts",
                         labels = function(x) format(x, scientific = TRUE),
                         expand = c(0,0),
                         breaks = scales::pretty_breaks(n = 10)) +
      scale_x_continuous(name = "Migration Time (Minutes)",
                         breaks = scales::pretty_breaks(n = 10))+
      ggtitle(paste(name, " EIE", " (m/z = ", mz,")",sep = ""),
              subtitle = paste("Data File: ", data_files[d])) +
      geom_ribbon(data = pf_df,
                  aes(x = mt.seconds/60, ymax = intensity, ymin = baseline, fill = peak.number),
                  alpha =0.4) +
      geom_text(data = ann_df,
                label = ann_df$peak.number,
                size  = 7,
                family = "sans",
                aes(x = peak.apex.seconds/60,
                    y = peak.height.counts + 0.08 * max_peak_height)) +
      geom_text(data = ann_df,
                label = ann_df$comment,
                size  = 7,
                family = "sans",
                aes(x = peak.apex.seconds/60,
                    y = peak.height.counts + 0.12 * max_peak_height)) +
      geom_segment(data = ann_df,
                   aes(x = peak.apex.seconds/60,
                       y = peak.height.counts + 0.04 * max_peak_height,
                       xend = peak.apex.seconds/60,
                       yend = peak.height.counts + 0.01 * max_peak_height),
                   arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
      theme(legend.position = "none",
            text = element_text(size = 25, family = "sans"))
    
    # Save plots to their respective folders within the "Plots" folder
    
    data_files_name <- list.files(path = "mzML Files")
    data_files_name <- gsub(".mzML", "", data_files_name, fixed = TRUE)
    
    ggsave(filename=paste(name,"_",data_files_name[d],".png",sep=""),
           width = 16,
           height = 9,
           plot = last_plot(),
           path = paste("Plots/", name_vec[n], sep = ""))
    
  }
  
  print("Plotting Complete")
  
  # 10. Export Data ----
  
  ## Generate peak area data frame ----
  
  peak_area_df <- cbind("file.name" = c(data_files_name[d],2:num_of_injections),
                        "peak.number" = c(1:num_of_injections),
                        peaks_df[,seq(from = 7, to = ncol(peaks_df), by = 7)])
  peak_area_df$file.name[2:num_of_injections] <- ""
  colnames(peak_area_df)[3:(length(name_vec) + 2)] <- name_vec
  
  # Update values to include <LOD and Interfered
  
  for (i in 1:num_of_injections){
    peak_area_df[i,3:ncol(peak_area_df)] <- ifelse(comment_df[i,] == "", peak_area_df[i, 3:ncol(peak_area_df)], comment_df[i,])
  }
  
  if(d == 1){
    peak_area_report = peak_area_df
  }else{
    peak_area_report = rbind(peak_area_report, peak_area_df)
  }
  
  ## Generate peak migration time data frame ----
  
  peak_mt_df <- cbind("file.name" = c(data_files_name[d],2:num_of_injections),
                      "peak.number" = c(1:num_of_injections),
                      peaks_df[,seq(from = 2, to = ncol(peaks_df), by = 7)] / 60)
  peak_mt_df$file.name[2:num_of_injections] <- ""
  colnames(peak_mt_df)[3:(length(name_vec) + 2)] <- name_vec
  
  if(d == 1){
    peak_mt_report = peak_mt_df
  }else{
    peak_mt_report = rbind(peak_mt_report, peak_mt_df)
  }
  
  ## Update progress bar ----
  
  progress <- paste(d,"/", length(data_files), "Files Completed")
  setWinProgressBar(pb, d, label = progress) 
  
  # Delete temporary mzml file
  
  file.remove(paste(data_files[d], "temp", sep = "_"))
  
  print(paste("Completed Analysis of Data File: ", data_file_names[d], sep = ""))
  print("")
  
  ## Export data ----
  
  write.csv(peak_area_report,
            file = "csv Outputs/Metabolite Peak Areas.csv",
            row.names = FALSE)
  
  write.csv(peak_mt_report,
            file = "csv Outputs/Metabolite Migration Times.csv",
            row.names = FALSE)
  
}

# close progress bar

close(pb)
