# PeakSeeker
## Software for Pre-Processing MSI-CE-MS Datasets

### Description

Multisegment injection-capillary electrophoresis-mass spectrometry (MSI-CE-MS) is a robust analytical technique capable of rapidly acquiring metabolomic data (<4min/sample) via serial injections of samples and quality controls within every analytical run. The high-throughput nature of this platform makes it an attractive technique for large scale metabolomic experiments where keeping costs down and meeting project deadlines are of critical importance. However, until now, data collected by MSI-CE-MS had to manually pre-processed by an experienced analyst due to the high migration time variability of CE-MS and the complexity of multiplexed data sets, which reduce the compatibility of this technique with other pre-processing tools. Unfortunately, manually pre-processing data is slow, expensive, and tedious, ultimately decreasing the merits of this technique. Thus, we are introducing PeakSeeker, an open-source software written in the R statistical environment for the automated pre-processing of targeted full-scan MSI-CE-MS data.

### Current Features

The key differentiating feature of PeakSeeker from other currently available software tools for pre-processing metabolomic datasets is its use of migration indexes to predict the elution time of analytes in MSI-CE-MS datasets. Thus, to achieve correspondence between analytical runs, PeakSeeker does not perform migration time alignments or converting the time dimension to electrophoretic mobility. Instead, PeakSeeker computes migration indexes for each analyte and uses the migration time of internal standards, or any reliable signals with sufficient signa-to-noise ratios, to compute the migration times of analytes which can then be used for peak annotation and integration. Additionally, as PeakSeeker was designed for multiplexed datasets, it also uses the spaces between analytical peaks to confirm and adjust peak annotation, as the gaps between peaks is typically consistent in MSI-CE-MS experiments. Results produced by PeakSeeker are saved and include:

1. A table containing the migration indexes or relative migration times used to annotating analytes
2. A copy of the parameters used to process the data
3. Plots of each extracted ion electropherogram which can be used to check for proper peak annotation and integration
4. A table containing the migration times of peaks or expected peak positions
5. A table containing the peak areas of each peak

### Usage

PeakSeeker has only two requirements to get you up and running:
  * Convert all your data files to open-source mzML files and save them all in a folder titled "mzML Files"
  * Provide a targeted mass list and corresponding parameters using the provided "Mass List and Parameters.xlsx" template

### Copyright


