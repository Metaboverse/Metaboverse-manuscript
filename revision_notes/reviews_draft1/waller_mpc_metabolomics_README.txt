Description of experiment

On 14-15 October 2014, T. Cameron Waller cultivated 3 biological replicate clones for each of 4 different strains in synthetic yeast medium with raffinose as a carbon source. Cameron used a block randomization strategy to randomize the sample designations of each culture and their sequence in the procedure, as recorded in the laboratory notebook and in file "2014-10-14_Randomization.xlsx". On 15 October 2014 at time point 0 hours, Cameron introduced uniform carbon-13 glucose to the cultures to a final concentration of 0.5 % weight/volume. At specific time points thereafter (0, 0.25, 0.5, 1, 2, and 5 hours after label introduction), Cameron harvested cells from these cultures and froze them, recording details on the harvest volumes in file "2014-10-15_CellCollection.xlsx". On 18-21 November 2014, J. Alan Maschek on the team of James E. Cox processed the samples (derivitization, gas chromatography, mass spectrometry) for metabolomics (original MassLynx sample list "AM2-105 Rutter Cam Nat Mut Long.spl"). J. Alan Maschek performed initial analyses of flux isotopic labeling of metabolites in the samples, and Cameron continued this analysis in order to increase the number of analytes available.

Cameron studied the chromatogram ion traces to identify more analytes from their gas chromatograph (GC) retention times and their mass-to-charge ratios in mass spectrometry (MS). He compiled this information in files "Analyte.xlsx" and "Target.xlsx". He used these parameters to define methods in MassLynx software to detect ion trace peaks for each analyte and quantify the area under each peak. He defined these methods for two different methods of peak detection: Standard and Apex. Cameron ran the quantification methods in QuanLynx, from which he saved the quantification files.

On 7 April 2015, Cameron performed a thorough growth assay of the same yeast strains in the same media.

Description of files

original path: original path to file on Cameron's computer
original file: original file name on Cameron's computer
novel file: new file name in package for transfer to Jordan A. Berg

original path: "/.../data/remote/archive/2016/research_rutter/ResearchNotebook/"
original file: "research_notebook.pdf"
novel file:
description: This file consitutes the laboratory notebook of T. Cameron Waller during his work on the team of Jared P. Rutter. Notes from 14-15 October 2014 are on pages 886-889. The notebook describes yeast strains and references several additional files with information on culture turbidities.

original path: "/.../data/remote/archive/2016/research_rutter/Results/Spectrophotometer_AmershamUltrospec/"
original file: "2015-04-07_CultureTurbidity_Complete.xlsx"
novel file:
description: This table records the culture turbidity measurements from the growth assay of the same yeast strains as were used in the flux metabolomics experiment.

original path: "/.../data/remote/archive/2016/research_rutter/Results/Calculations/"
original file: "2014-10-14_Randomization.xlsx"
novel file:
description: This table details the block randomization of biological replicate clones of yeast strains for flux isotopic metablic labeling procedure.

original path: "/.../data/remote/archive/2016/research_rutter/Results/Calculations/"
original file: "2014-10-14_Sequence_With_Strain.xlsx"
novel file:
description: This table details the block randomization of biological replicate clones of yeast strains for flux isotopic metablic labeling procedure.

original path: "/.../data/remote/archive/2016/research_rutter/Results/Calculations/"
original file: "2014-10-15_CellCollection.xlsx"
novel file:
description: This table records culture densities at times of harvest and volumes of cultures harvested for metabolomics.

original path: "/.../data/remote/archive/2016/research_rutter/Projects/2016_Waller_Metabolomics/Method/Target/"
original file: "Analyte.xlsx"
novel file:
description: This table documents retention times and mass-to-charge ratios of analytes that I identified in the data.

original path: "/.../data/remote/archive/2016/research_rutter/Projects/2016_Waller_Metabolomics/Method/Target/"
original file: "Analyte_Identification.xlsx"
novel file:
description: This table documents characteristic mass-to-charge peaks for identification of analytes.

original path: "/.../data/remote/archive/2016/research_rutter/Projects/2016_Waller_Metabolomics/Method/Target/"
original file: "Target.xlsx"
novel file:
description: This table documents retention times and mass-to-charge ratios of analytes that I identified in the data. I used this table as a program-readable script to control a program that I wrote for AutoIt, a GUI automation program, to build method files in MassLynx. The program would read the table and then use the values in columns to set parameters in the MassLynx method file. I still have the AutoIt script, but it would need to be re-optimized for the screen on which it was used. It might be possible to use these same parameters to define method files for MassLynx in a more convenient way, such as through some tabular interface for MassLynx or another program capable of peak detection and quantification from raw GC-MS data. Alternatively, it would also be possible to use the QuanLynx reports that I still have.


original path: "/.../data/remote/archive/2016/research_rutter/Projects/2016_Metabolite_Analysis/Quantity_Method/"
original file: "Script_Python_QuanLynx_Parser.py"
novel file:
description: This script parses export files from QuanLynx. Or it at least used to work.

original path: "/.../data/remote/archive/2016/research_rutter/Projects/2016_Waller_Metabolomics/Results/Quantity_2015-07-25/"
original file: ...
novel file:
description: This directory contains the most recent reports that I have from QuanLynx. I used parameters (analyte retention times, mass-to-charge ratios, mass error windows) from files "Analyte.xlsx" and "Target.xlsx" to program a method in MassLynx for both Standard and Apex methods of peak detection. Then I used those methods in QuanLynx to detect and quantify analyte peaks in the raw chromatography mass spectrometry data from the Waters GC-MS. I think that this latest quantification of the data did not include all the samples. I think it only included the samples for the knock-out and wild type rescue strains. I have an old Python script that I used to parse the QuanLynx export files. The QuanLynx export also includes the information from the MassLynx sample lists, so it indicates the original sample number. It is then necessary to reference the block randomization to determine to which yeast strain each sample matches.

----------
----------
----------
original path: "/.../data/remote/archive/2016/research_rutter/Results/Mass_Spectrometry/2014-10-15/"
directory size: 2.3 Gigabytes

original path: "/.../data/remote/archive/2016/research_rutter/Results/Mass_Spectrometry/2014-10-15/Sample/"
original file: ...
novel file:
description: This directory contains the original sample lists from the GC-MS runs in MassLynx format. These sample lists are necessary to match the raw data to the original samples.

original path: "/.../data/remote/archive/2016/research_rutter/Results/Mass_Spectrometry/2014-10-15/Data/"
original file: ...
novel file:
description: This directory contains the original MassLynx data files from the GC-MS runs, of which there were apparently 168 (different injection ratios of 10:1 or 100:1 and short or long chromatography runs).

original path: "/.../data/remote/archive/2016/research_rutter/Results/Mass_Spectrometry/2014-10-15/Method/Quantity_2015-07-25/"
original file: ...
novel file:
description: This directory contains the MassLynx method files for QuanLynx quantification. I think that these methods correspond to the QuanLynx export files that I still have from 25 July 2015.

