# New_Tormes_Report

Version 1.0

### Note: Currently this reporting option doesn't work with less than 3 sequences or sequence runs where you have requested to skip the pangenome / phylogeny steps of the pipeline

### TO DO
 * allow for single/double genome reporting or where user has opted to skip pangenomics and phylogeny steps.
 * add virulence and genera reporting options

### Usage

Install dependencies either in a conda environment or virtual environment  

```
conda create -n tormes_report -c conda-forge datapane=0.14.0 plotly biopython
```

After running the tormes pipeline, unpack tormes report files in your tormes results folder:  

```
tar -xvzf tormes_report.tgz
```

Activate the tormes_report environment
```
conda activate tormes_report
```

copy tormes_report-1.0.py into the tormes_report folder

then run it:  

```
python tormes_report-1.0.py
```

### Output

tormes_report_datapane.html
