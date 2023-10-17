# SNAPMe
### <ins>S</ins>urveying <ins>N</ins>utrient <ins>A</ins>ssessment with <ins>P</ins>hotographs of <ins>Me</ins>als

This repository contains the analytic code and data used for reproducing results and visualizations from the submitted manuscript [citation to follow]

#### Use and Requirements

Each Python script provides addition information on the input and output files.

Python scripts use the following package versions:
- python==3.8.8
- pandas==1.4.4
- nltk==3.7

R scripts use the following package versions:
- R==4.1.0
- ggplot2==3.3.5
- dplyr==1.0.8
- boot==1.3-28
- RVAideMemoire==0.90-81-2

#### Supplemental files
Supplemental files provide additional information:
- File 1: ingredientized ASA24 food descritpions and which database was used for this process (NA if ingredientization was not done)
- File 2: predictions for each food photo using the Inverse Cooking algorithm and F1 scores
- File 3: predictions for each food photo using the Im2Recipe algorithm and F1 scores
