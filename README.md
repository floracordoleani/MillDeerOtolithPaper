# MillDeerOtolithPaper

This repository provides the data and R codes used to perform the spring-run Chinook otolith analyses and to produce the figures presented in the manuscript "Threatened salmon rely on a rare life history strategy in a modified and warming landscape".

STEP I: Set up.

Follow these steps to download software, model input files, additional code, and libraries required to replicate our study.

    Download R and RStudio.

    Get model input files and additional code at https://github.com/floracordoleani/MillDeerOtolithPaper.

    Either clone the repository or download as MillDeerOtolithPaper.zip; when unzipped locally, this directory will serve as your R project directory.

    Confirm that files are stored with the following structure within the MillDeerOtolithPaper directory.

    code
        MasterCode.R - this is the main analysis script.
        ClusterAnalysis.R - this script contains all the model functions and is sourced from the model script.

    Data.in
        MillDeerIncrements.csv - 
        MillDeerNatalExit.csv - 
        MillDeerOtoliths.csv -
        OR_FL_FINALforR.csv -
        RSTChinMillDeer.csv - 

    Create a new R project in RStudio.

STEP II: Run the clustering analysis.

The clustering analysis code will produce the data and figure listed below:

    Data.out
    MillDeerClusters.csv - 

    Figures
    FigureS1.png - 
    
STEP III: Run the otolith analyses.

The otolith analysis code will produce the figures listed below:

    Figures
    Figure1.png - 
    Figure2.png - 
    Figure3.png - 
    Figure4.png - 
    FigureS2.png - 
    FigureS3.png - 

