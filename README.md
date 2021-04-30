# MillDeerOtolithPaper

This repository provides the data and R codes used to perform the spring-run Chinook otolith analyses and to produce the figures presented in the manuscript "Threatened salmon rely on a rare life history strategy in a modified and warming landscape".

Authors: Flora Cordoleani, Corey C. Phillis, Anna Sturrock, Alyssa M. FitzGerald, George Whitman, Anthony Malkassian, Peter K. Weber, and Rachel C. Johnson.

STEP I: Set up.

Follow these steps to download software, model input files, additional code, and libraries required to replicate our study.

    Download R and RStudio.

    Get model input files and codes at https://github.com/floracordoleani/MillDeerOtolithPaper.

    Either clone the repository or download as MillDeerOtolithPaper.zip; when unzipped locally, this directory will serve as your R project directory.

    Confirm that files are stored with the following structure within the MillDeerOtolithPaper directory.

    code
        MasterCode.R - this is the main otolith analysis script.
        ClusterAnalysis.R - this is the otolith strontium profile clustering analysis script.

    Data.in
        MillDeerOtoliths.csv - Otolith strontium isotope data 
        MillDeerIncrements.csv - Otolith increment data from the microchemistry analysis
        MillDeerNatalExit.csv - Fish's otolith radius at natal exit (estimated visually) data
        OR_FL_FINALforR.csv - Central Valley juvenile size and otolith radius data for spring-run fish size back-calculation model
        RSTChinMillDeer.csv - Combined Mill and Deer Creek juvenile rotary screw trap data provided by the California Department of Fish and Wildlife

    Create a new R project in RStudio.

STEP II: Run the clustering analysis.

The clustering analysis code will produce the data and figure listed below:

    Data.out
    MillDeerClusters.csv - Juvenile life history type classification table

    Figures
    FigureS1.png - Strontium profile groups as identified by the clustering analysis 
    
STEP III: Run the otolith analyses.

The otolith analysis code will produce the figures listed below:

    Figures
    Figure1.png - Spring-run Chinook salmon life history diversity
    Figure2.png - Migrant size distributions at natal and freshwater exit
    Figure3.png - Migrant sizes and life history diversity across years
    Figure4.png - Early-life (first 15 days after emergence) salmon growth across life history types
    FigureS2.png - juvenile trapping data, fish size reconstruction, migrant size and timing distributions at natal and freshwater exit
    FigureS3.png - Early-life (first 30 days after emergence) salmon growth across life history types

