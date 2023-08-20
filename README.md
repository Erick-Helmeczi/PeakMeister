# PeakSeeker
## Software for Pre-Processing MSI-CE-MS Datasets

### Description

Multisegment injection-capillary electrophoresis-mass spectrometry (MSI-CE-MS) is a robust analytical technique capable of rapidly acquiring metabolomic data (<4min/sample) via serial injections of samples and quality controls within every analytical run. The high-throughput nature of this platform makes it an attractive technique for large scale metabolomic experiments where keeping costs down and meeting project deadlines are of critical importance. However, until now, data collected by MSI-CE-MS had to manually pre-processed by an experienced analyst due to the high migration time variability of CE-MS and the complexity of multiplexed data sets, which reduce the compatibility of this technique with other pre-processing tools. Unfortunately, manually pre-processing data is slow, expensive, and tedious, ultimately decreasing the merits of this technique. Thus, we are introducing PeakSeeker, an open-source software written in the R statistical environment for the automated pre-processing of targeted full-scan MSI-CE-MS data.

### Current Features

1. The software allows users to process their data with minimal R experience as all required parameters can be adjusted within an excel sheet.
2. Mass calibration to lock masses is performed to minimize mass varibaility which commonly occurs in time-of-flight mass spectrometers, expecially during long studies. 
3. Extracted ion electropherograms (EIEs) are extracted from the calibrated data files using the user supplied mass-to-charge list and then smoothed.
4. Peaks are detected and filtered using user defined parameters.
5. Plots are generated and saved for users to review to ensure correct peak peaking and integration occured. 
6. Peak area and migration time data is further summarized in csv files.
7. A windows progress bar and print statements have been implimented so users can have an estimate of completetion. 
