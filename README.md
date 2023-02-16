# MSI-CE-MS-Data-Preprocessing

### Description

The motivation behind this program was to develop software that can reduce the burden, cost, and time required to manually process multisegment injection-capillary electrophoresis-mass spectrometry (MSI-CE-MS) data. While, MSI-CE-MS allows for the rapid acquisiton of metabolomic data, we were unable to find data processing software which could be tailored to the unique output of this technique. Thus, the key features of this software is to take raw (.mzML) data files and user supplied conditons, such as a targted mass-to-charge list, and output plots and tables of peak areas and migration times. 

### Current Features

1. The software, which is solely written in R, allows user to process their data with minimal R experience as all parameters used within the program can be changed within excel sheets accompaning the project. 
2. Mass calibration to lock masses is performed to minimize mass varibaility which commonly occurs in time-of-flight mass spectrometers expecially during long studies. 
3. Extracted ion electropherograms (EIEs) are extracted from the calibrated data file using the user supplied mass-to-charge list and smoothed.
