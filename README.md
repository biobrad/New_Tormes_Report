# New_Tormes_Report
Work in progress of recreating Tormes pipeline report using python

This is a work in progress, not complete, is functional though.

### Usage

Install dependencies either in a conda environment or virtual environment  

```
conda create -n tormes_report -c conda-forge datapane=0.14.0 plotly biopython
```

After running the tormes pipeline, unpack tormes report files:  

```
tar -xvgf tormes_report.tgz
```

copy tormes_report_wip.py into the tormes_report folder

then run it:  

```
python tormes_report_wip.py
```
