# Nek5000 Resolvent Tools : Resolvent modes, SSRM

##Setup

Copy in the folder relevant .rea/.map files for the flow of interest. In the .rea file 

1. set a constant time step (**DT/param(12)**<0)
2. set **npert** to -1
3. set a restart option to load the flow for which the analysis is to be performed. The restart file should have step and time equals to zero.

Setup the *params.input* file, which has 5 lines, namely

1. **runType**:  <0 to integrate an initial condition in time until it vanishes (for relaxation time estimation), or >0 for iterations (automatically set)
2. **relaxTime** : relaxation time,
3. **maxFreq** : highest frequency of interest,
4. **nFreq** : number of frequencies (frequencies will be linearly spaced between 0 and the max frequency). Note that the maximum number of frequencies that can be used is set by the parameter **nmaxFreq** in *FFT.f*. Increase it if needed.
5. **norm** : norm at which inputs are normalised before each iteration (a value of 1 is the usual choice)

To define the relaxation time set **runType** to -1 (direct run) or -2 (adjoint run) in *params.input*, and run Nek5000. The perturbation field will be randomly initialised, and the perturbation norm stored in *runNorm.txt*. Check the time needed for the norm to reach a desired tolerance and modify *params.input* accordingly.  

To use  **B** and **C** matrices create 'B.fld' and 'C.fld' files. The velocity fields will be used as diagonal entries into B and C. If the files are not present, B and/or are assumed to be the identity matrix.

Finally run iterations to compute resolvent modes as
>  `./runArnoldi [reafile] [nproc] [firstIter] [lastIter] [nOctSections(optional)]`


* **[reafile]** : name of the file to be run
* **[nproc]** : number of processes
* **[firstIter]/[lastIter]** : iteration range to be run. E.g. set it to "0 4" to run the first 5 iterations. Supplementary  iterations can be obtained on a second run setting those to "5 10".
* **[nOctSections]** : number of parallel octave sections used to compute the inputs for the next iteration. The default is **[nproc]**. For large problems and large interaction count, memory may limit the number of frequencies that can be simultaneously processed.  Q smaller **[nOctSections]** reduces the memory requirements, at the price of a longer iteraction setup time.

Each iteration is saved on folders "IterArXX", with subfolders for the direct and adjoint run. Inputs and outputs of each run are saved and used latter for post-processing. Gains and modes are computed and saved on the last iteration, and as stored as

*  **./Gains_SS_XXX.dat** : containing to gains estimated using the first XXX iterations,
*  **./Resolvent/resresp_XX_[real/imag]0.fYYYYY** : XX-ith leading real/imaginary part of the response modes for the YYYYY-th frequency computed,
*  **./Resolvent/resforce_XX_[real/imag]0.fYYYYY** : likewise for forces.

## Comments
Nek5000 .rea parameter IOSTEP is copied into IOPERT during run time, after   IOSTEP and IOTIME are set to zero, IOPERT is used for saving perturbation snapshots. This avoids writing the baseflow to disk.

Other than for debug purposes, set IOSTEP and IOTIME to zero to minimise disk usage. 


##Examples


Examples for the computation of modes for 2/3D channel flows are. To reproduce validation generate .map files for *channel_[2,3].rea*, modify *SIZE* to setup Nek5000 for a 2/3D run, compile and run *runArnoldi* as
> `makenek channel`   
> `./runArnoldi channel_[2/3]D [number of cores] 0 20`

Alternatively, run the **runExamples** file.

##Files
The main files, other than the inputs previously mentioned, are

* **run** : A bash script that calls nek5000, to compute direct/adjoint runs and to apply filters, and octave, to compute filter coefficients.
* **Nek_TRM.ust** : Nek5000 user file to compute runs.
* **Nek_GetFilter.m** : Octave/Matlab script that reads previous runs frequency norms and designs a filter to normalise its spectra. 
* **Nek_ClearFiles.m** : Octave/Matlab script that reads snapshot norms and lists snapshots with norms bellow the tolerance for deletion.
* **readnek.m** and **writenek.m** : routines to load and write files for nek5000 in matlab/octave. These were obtained from [another project](https://github.com/nfabbiane/nekmatlab).

# References

1. Martini, Eduardo, Daniel Rodríguez, Aaron Towne, and André V. G. Cavalieri. “Efficient Computation of Global Resolvent Modes,”  Journal of Fluid Mechanics 919 (2021): A3. [https://doi.org/10.1017/jfm.2021.364](https://doi.org/10.1017/jfm.2021.364)



