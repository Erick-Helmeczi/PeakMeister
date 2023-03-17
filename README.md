# Peak Seeker
## Software for Automating MSI-CE-MS Data Preprocessing

### Description

The motivation behind this program was to develop software that can reduce the burden, cost, and time required to manually preprocess multisegment injection-capillary electrophoresis-mass spectrometry (MSI-CE-MS) data. While MSI-CE-MS allows for the rapid acquisiton of metabolomic data, we were unable to find data processing software which could be tailored to the unique output of this technique. Thus, the key objective of this software is to take raw (.mzML) data files and user supplied conditons, such as a targted mass-to-charge list, and output plots and tables of summarized peak areas and migration times. 

### Current Features

1. The software allows users to process their data with minimal R experience as all required parameters can be adjusted within an excel sheet.
2. Mass calibration to lock masses is performed to minimize mass varibaility which commonly occurs in time-of-flight mass spectrometers, expecially during long studies. 
3. Extracted ion electropherograms (EIEs) are extracted from the calibrated data files using the user supplied mass-to-charge list and then smoothed.
4. Peaks are detected and filtered using user defined parameters.
5. Plots are generated and saved for users to review to ensure correct peak peaking and integration occured. 
6. Peak area and migration time data is further summarized in csv files.
7. A widnows progress bar and print statements have been implimented so users can have an estimate of completetion. 
