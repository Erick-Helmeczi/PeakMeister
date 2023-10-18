# PeakMeister
## Software for Pre-Processing MSI-CE-MS Datasets

### Description

Multisegment injection-capillary electrophoresis-mass spectrometry (MSI-CE-MS) is a robust analytical technique capable of rapidly acquiring metabolomic data (<4min/sample) via serial injections of samples and quality controls within every analytical run. The high-throughput nature of this platform makes it an attractive technique for large scale metabolomic experiments where keeping costs down and meeting project deadlines are of utmost importance. However, until now, data collected by MSI-CE-MS had to be manually pre-processed by an experienced analyst due to the high migration time variability of CE-MS and the complexity of multiplexed data sets, which reduce the compatibility of this technique with other pre-processing tools. Unfortunately, manually pre-processing data is slow, expensive, and tedious, ultimately decreasing the merits of this technique. Thus, we are introducing PeakMeister, an open-source software written in the R statistical environment for the automated pre-processing of targeted full-scan MSI-CE-MS data.

### Current Features

The key differentiating feature of PeakMeister from other currently available software tools for pre-processing metabolomic datasets is its use of migration indexes to predict the elution time of analytes in MSI-CE-MS datasets. Thus, to achieve correspondence between analytical runs, PeakMeister does not perform migration time alignments or converting the time dimension to electrophoretic mobility. Instead, PeakMeister computes migration indexes for each analyte and uses the migration time of internal standards, or any reliable signals with sufficient signal-to-noise ratios, to compute the migration times of analytes which can then be used for peak annotation and integration. Additionally, as PeakMeister was designed for multiplexed datasets, it also uses the spaces between analytical peaks to confirm and adjust peak annotation, as the gaps between peaks are typically consistent in MSI-CE-MS experiments. Results produced by PeakMeister are saved and include:

1. A table containing the migration indexes or relative migration times used to annotating analytes
2. A copy of the parameters used to process the data
3. Plots of each extracted ion electropherogram which can be used to check for proper peak annotation and integration
4. A table containing the migration times of peaks or expected peak positions
5. A table containing the peak areas

### Quick Start

PeakMeister has only two requirements to get you up and running and all required template files and folders are included in the latest releases:
  * Convert all your data files to open-source mzML files and save them in the folder titled "mzML Files"
  * Provide a targeted mass list and the corresponding parameters using the provided "Mass List and Parameters.xlsx" template

Use the project file to open the R script titled "code" and execute the script to begin pre-processing your data

### Detailed Usage

Although PeakMeister is a script based software tool, users will typically require little to no knowledge of coding to pre-process their data as an excel sheet containing all pre-processing parameters is provided. In this detailed usage overview, the purpose of each parameters and how they can be manipulated to accurately annotate MSI-CE-MS data will be explained.

##### Sheet 1: Mass List

1. name - Provide a name for each analyte in your study. All names must be unique, so if you have multiple unknown compounds, use names such as Unknown-1, Unknown-2, Unknown-3, etc.
2. mz - Provide the mass-to-charge to be extracted for each analyte in your study.
3. extraction.window.ppm - This is the extraction window used to extract each mz value provided with units of ppm. We found that a minimum mass window of ~30-35 ppm is required depending on the mass of the analyte, however this will likely be dependent on the mass spectrometer used duing data acquisition.
4. interference - Here you can designate analyte interferences. Write the name of the analyte or internal standard that is a significant interference and PeakMeister will annotate any instances of overlap.
5. interference.comigration.threshold - This parameter only needs to be set when an interference is also provided. This is the window used to determine if two peaks are overlapping enough to be considered an interference. The widths of your peaks will determine how high this value needs to be set.
6. minimum.peak.width.seconds - This parameter is used to filter out noise during peak detection. Only peaks with a width equal to or greater than this value will be considered during annotation.
7. migration.window.seconds - After migration time prediction, PeakMeister will look for analyte peaks at the expected migration time +/- the time set for this parameter.
8. peak.space.tolerance.percent - Peaks are expected to be spaced approximaely equally apart. This parameter designates how differently they can deviate from the median spacing before being subject to reanalysis by PeakMeister. 
9. snr.threshold - After identifiying which signals correspond to the analyte peaks, PeakMeister will compute their S/N to determine if they should be recorded as "<LOD".
10. smoothing.kernal - Smoothing is performed using the "ksmooth" function from the [stats](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/ksmooth.html) package. This parameter sets the kernal to be used.
11. smoothing.strength - Smoothing is performed using the "ksmooth" function from the [stats](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/ksmooth.html) package. This parameter sets the bandwidth to be used.
12. The last columns of of this sheet are where peak migration times for each analyte are to be set. All of these values need to be taken from the same run, including those that will be set for the internal standards in the next sheet.

##### Sheet 2: Reference Mass List

For this sheet, I will focus on the sections not discussed above. Features used in this section can be used as landmark peaks that can help predict the migration times of analytes in the previous sheet. These landmark peaks use a separte algorithm for annotation which requires their peak areas to be reliable and always the most intense within the EIE. Thus, internal standards are typically used here, however analytes that are very concentrated in your sample matrix may also be used. 

1. class - Use "Reference" for any compound you want to be used to calculate migration indexes for analytes on the previous sheet. Typically this should be every compound listed here, however the option is availble to not use them. Simply type "Not Reference" instead if this is the case.
2. min.mt.min - Minimum migration time cutoff. Peaks that elute before this threshold will not be considered as possible reference standard peaks.
3. max.mt.min - Maximum migration time cutoff. Peaks that elute after this threshold will not be considered as possible reference standard peaks.

##### Sheet 3: Parameters

1. number.of.injections - Number of injections used during data acquisition.
2. ref.mass.one - Lower lock-mass used for accurate mass correction.
3. ref.mass.two - Upper lock-mass used for accurate mass correction.
4. ref.mass.window.ppm - mass window to search for lock-mass peak.
5. ref.mass.counts - Minimum peak hight requirement of the lock-mass for correction to be applied.
6. apply.mass.correction - Should the lock mass correction be applied. Correction will only be applied if "Yes" is set for this parameter.
7. apply.smoothing - Should smoothing be applied to the EIEs. Smoothing will only be applied if "Yes" is set for this parameter.
   
### Copyright

PeakMeister is licensed under the [MIT](https://choosealicense.com/licenses/mit/) license

