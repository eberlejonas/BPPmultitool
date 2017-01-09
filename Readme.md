# BPPmultitool
BPPmultitool is a collection of R scripts that is useful to analyze the output of BPPmulti which can be found here:
http://github.com/eberlejonas/BPPmulti

Currently supported is:
- automatic reading of BPPmulti settings
- plotting guide trees
- parsing ESS (effective sampling size) values
- writing of summarizing result tables
- extraction of single analysis results
- plots of single analysis results
- summarizing plots of species delimitation results from multiple guide trees and repeated runs

## Requirements
An R installation with the packages *ape* and *stringr* installed. Use `install.packages("ape")` and `install.packages("stringr")` in the R console if you don't have them installed already.

## Usage
1. Download all BPPmultitool files to a folder on your computer.
2. Open template.R in a texteditor
3. Search for the command `setwd("/home/jonas/my_BPPmulti_run")` (line 21) and put the path to the working directory inside the quotes (the one you also specified in your BPPmulti analysis). On Windows you can for isntance write `setwd("C:/Users/jonas/my_BPPmulti_run")`.
4. Execute all desired commands in template.R on the R console. `source("subs.R")`, `setwd("...")`, `settings <- getSettings()`, and `BPPresults <- parseBPPmulti(settings)` are obligatory.

Please also refer to the comments in template.R.
