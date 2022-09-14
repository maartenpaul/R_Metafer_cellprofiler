# R Metafer cellprofiler
Scripts for quantification of DSB repair foci in fixed samples acquired with Metafer WideField microscope
This repository contains Cell profiler pipelines for analysis and a R script to analyze and visualize the data

## Data format
For automatic batch processing in Cell Profiler, files should be properly named.
The current file name format:

"Cell line-treatment-replicate[0-9].ImageNumber.Channel.TIF"

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
