# CyTOF Processor

The CyTOF Processor R Shiny application was created to provide a supplementary interface for the gateR package. 


## Installation

### The application uses the following packages:

| Package    | Version |
|------------|---------|
| dplyr      | 1.0.2   |
| flowCore   | 2.0.1   |
| gateR      | 0.1.9   |
| ggplot2    | 3.3.2   |
| logR       | 1.2.3   |
| shiny      | 1.5.0   |
| shinyFiles | 0.9.0   |

In order to launch the application, use the following command in an R console.
```
runGitHub( "CyTOF", "SubiNair", ref = "main")
```

## CyTOF Data Prep

Input Flow Cytometry data needs to be formatted to the [gateR](https://cran.r-project.org/web/packages/gateR/gateR.pdf) package's specifications. The gateR shiny app only requires  .fcs uploads and will reformat the data for gating. The following format is accepted by gateR. In the following example, the ID column for each of the cells is first, followed by the conditions, and the markers (each marker is associated with a channel which what the columns are named after):

|  ID              | Condition 1   |  Condition 2                    |Channel 1             |Channel ...            
|----------------|-------------------------------|-----------------------------|---|--|
|Cell 1|0           |control          |*x* | ...
|Cell 2          |10    |control         |*y* | ...
|Cell 3          |0|case|*z* | ...

*x*, *y*, and *z* in this example represent the numeric cell specific data from the experiment. The shiny application additionally appends a column for the file each cell came from. This is for tracking and is not used in the analysis. The condition data is extracted from the separate metadata.csv file that is required.

The `fcsprocessor.R` file is designed to take raw .fcs data and transform it to meet the requirements of the package.

The `fcsprocessor.R` script is designed to take raw .fcs data and transform it to meet the requirements of the package. The application uses the same script. 

The application uses the [flowCore](https://bioconductor.org/packages/release/bioc/vignettes/flowCore/inst/doc/HowTo-flowCore.pdf) package to read in the .fcs data and store it in a `flowFrame` object.  

 
