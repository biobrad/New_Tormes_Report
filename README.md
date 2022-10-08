# New_Tormes_Report

Version 1.5
### Improvements
* largely code based, no cosmetic improvements:
* script now automatically decompresses and re-compresses report_files.tgz
* script will not fail if pangenome and phylogenetic details are not present (single, double genome analysis)

### TO DO
 * add virulence and genera reporting options

### Usage

Install dependencies either in a conda environment or virtual environment  

```
conda create -n tormes_report -c conda-forge datapane=0.14.0 plotly biopython
```

After running the tormes pipeline, copy tormes_report-1.5.py into your tormes output folder.

Activate your tormes_report environment

```
conda activate tormes_report
```

then run the python script:  

```
python tormes_report-1.5.py
```

### Output

tormes_report_datapane.html
