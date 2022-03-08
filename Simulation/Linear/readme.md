# Readme file for linear simulation

DM_sim_linear_pipe_batch_parallelised_DAJ.R sets up the simulation scenarios  
DM_sim_linear_1obperday_final_parallelised.R is the code for each iteration  
dma_linear_DAJ_source.R is the edited dma package that allows batch updating  
submit.qsub is the bash script used to execute the simulation  

The simulation was run on the computational shared facility (http://ri.itservices.manchester.ac.uk/csf3) at Manchester University.  
All files were placed in the same folder and then the submit.qsub file was submitted. This file is BASH script and executes the R files.  
It first runs the DM_sim_linear_pipe_batch_parallelised_DAJ.R file. This sets up a martix with all of the parameter choices and runs the DM_sim_linear_1obperday_final_parallelised.R code. This generates a function (called "DM.sim.linear.1") that runs a single iteration of the simulation given a set of paremeters.  
The batch file then reads the matrix and executes the DM.sim.linear.1 with each row (unique set of parameter choices).  

In the submit.qsub file please ensure the number at the end of line 8 is the number of rows in the matrix. In the simulation I have 400 parameter choices when we observe one observation at each time point and 1000 iterations. Hence, 400000 rows are in my matrix.


