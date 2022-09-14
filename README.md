# R Metafer cellprofiler
Scripts for quantification of DSB repair foci in fixed samples acquired with Metafer WideField microscope.
This repository contains Cell profiler pipelines for analysis and a R script to analyze and visualize the data.

Issues, questions or suggestions can be reported here: https://github.com/maartenpaul/R_Metafer_cellprofiler/issues

## Data format
For automatic processing of image files in Cell Profiler, files should be properly named.
The CellProfiler pipeline uses the following file name format, which can be adjusted in Cell Profiler:
It is recommended to already stick with this file name format at the Metafer microscope so you do not need to rename file names.

"Cell line-treatment-replicate[0-9].ImageNumber.Channel.TIF"

## Running CellProfiler
The pipeline can be loaded imported in CellProfiler via "File->Import"

Subsequently you can drag the images you want to analyze into the Images window.

At the next step "Metadata" you can check if the file names are properly processed

You can use "Start Test Mode" to test if the pipeline is working good enough for your images.

Depending on the output you can change the thresholds or other parameters to properly segment the nuclei and to detect the foci. It is recommended to check the output for a few images from different conditions if the results are satisfactory.

When running the entire pipeline images will we exported that 


## Batch analysis
When analyzing large data sets it happens frequently that Cell Profiler runs out of memory. I therefor found a way to split up the data and run in parallel. This work for Linux. A solution for Windows needs to be sorted out.

To run cellprofiler on the Linux analysis pc : open terminal

Run:

```bash
cellprofiler
```

Load the project and select the files to analyze, so all the cannels

- Check if the file names are properly processed in the [Metadata] step
- Make sure to use groups to export separate csv file per slide and run the different slides in parallel to avoid memory issues
- Make sure to include as last step in analysis - Batch processing - this will export a batch.hdf5 file
- Change in ExportToSpreadsheet [Filename prefix] to correct data

for parallel analysis:

Type in command line

```bash
cd ~/batch/
cellprofiler --get-batch-commands Batch_data.h5
```

example output:

```bash
CellProfiler -c -r -p Batch_data.h5 -f 1 -l 109
CellProfiler -c -r -p Batch_data.h5 -f 110 -l 176
CellProfiler -c -r -p Batch_data.h5 -f 177 -l 239
CellProfiler -c -r -p Batch_data.h5 -f 240 -l 298
...
```

Copy output to Libreoffice Calc/Excel

Copy rows with the numbers and include in the following order in the code:

Run this code from command line

```bash
parallel --xapply -j 6 cellprofiler -c -r -p batch/Batch_data.h5 -f {1} -l {2} ::: 1 110 177 240 ::: 109 176 239 298
```

## Processing and plotting data in R
