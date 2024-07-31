# Modelling the immunogenicity, efficacy, and impact of dengue vaccines

This repository provides the code and data to reproduce the main results of the PhD thesis Modelling the immunogenicity, efficacy, and impact of
dengue vaccines. Each folder entitled Chapter 2, Chapter 3, etc., contains the code to reproduce the results of that chapter. Scripts should be run in their numbered order. They take as input anything within the data folders, and rely on functions provided in the R folders. All results and figures are saved to the folder outputs. Code for each chapter is set up to provide a demo of the results that are run. To reproduce the results presented in the thesis, the Stan models of Chapters 2 and 3 were run for 10,000 iterations, the impact model was run for 10,000 simulations per scenario in Chapter 4, and the machine learning classifiers were fitted to 100 bootstrap sets of the test and train data in Chapter 5.  
