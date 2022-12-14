# New_Tormes_Report

Alternative reporting option for Tormes Genome Pipeline results:  
https://github.com/nmquijada/tormes

Version 1.75
### Improvements
* code based, use of list instead of dictionary for collecting visualisation objects

### TO DO
 * add virulence and more genera reporting options

### Usage

Install dependencies in a conda environment  

```
conda create -n tormes_report -c conda-forge datapane=0.14.0 plotly biopython scipy
```

After running the tormes pipeline, copy tormes_report-1.75.py into your tormes output folder.

Activate your tormes_report environment

```
conda activate tormes_report
```

then run the python script:  

```
python tormes_report-1.75.py
```

### Output

tormes_report_datapane.html

![image](https://user-images.githubusercontent.com/55652506/197070644-118ed30d-9023-4801-bbef-8ab5e73fd9c6.png)
![image](https://user-images.githubusercontent.com/55652506/197070720-1cbb522f-fd7d-4fba-ae76-2a86ff8ef3e6.png)
