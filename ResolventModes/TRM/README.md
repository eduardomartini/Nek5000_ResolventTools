# Nek5000 Resolvent Tools : Resolvent modes, SSRM

##Setup

Copy in the folder relevant .rea/.map files for the flow of interest. In the .rea file 

1. set a constant time step (**DT/param(12)**<0)
2. set **npert** to -1
3. set a restart option to load the flow for which the analysis is to be performed. The restart file should have step and time equals to zero.

Setup the *params.input* file, which has 6 lines, namely

1. **runType** :  -1 to apply FIR filter, positive even for direct runs and positive odds for adjoint runs
2. **time** : step used. Overwrites timestep in .rea file. Used also for filter construction.  
3. **IOSTEP** : also overwrites the one one in .rea file). The sampling period is IOSTEP*DT
4. **fcutoff** : Cut-off for filter construction.
5. **dfTar** : Target filter resolution
6. **tol**  : Norm tolerance. Runs are time marched until the relative norm (curr norm/max norm) becomes smaller than the tolerance.

Finally, input the list of desired frequencies in *freqList.input*. Note that the maximum number of frequencies that can be used is set by the parameter **nmaxFreq** in *FFT.f*. Increase it if needed.

To run iterations
> ` ./runArnoldi [reafile] [nproc] [firstIter] [lastIter] [keepfiles:y/n(optional)]`

* **[reafile]** : name of the file to be run
* **[nproc]** : number of processes
* **[firstIter]/[lastIter]** : iteration range to be run. E.g. set it to "0 4" to run the first 5 iterations. Supplementary  iterations can be obtained on a second run setting those to "5 10".
* **[keepfiles]** : "y" to save iterations time histories, "n" to only keep the last one. The last time history is needed on the next iteration. Previous can be used for post-processing extra frequencies.

Each iteration is saved on folders "IterXX". Even XX corresponds to a direct run, and odd XX to adjoin runs: a full iteration corresponds to two folders. Fourier transforms each iteration for the time series are stored in their corresponding folders as *[c/s]01[reaFile]0.fXXXXX*, corresponding to the cosine and sine components of the XXXXX-th frequency.

After the last iteration is completed the *run* script calls *computeGains.m* and saves power-iteration (*gains_powerIter.twt*) and Krylov-space (*gains_krylov_XXX.txt*) gains estimations, and norms for runs before (*runNormss.dat*) and after (*run_freqNorms_filtered.dat*) application of the temporal filter.

Check that the filter parameters used are sufficient for proper regularisation of the amplitudes. plotGains plots frequency norms for all runs, the ratio between the max and min amplitudes should ideally be lower than ~100.

Leading response/force modes estimated with the power-iteration  are the Fourier components saved the last direct/adjoint run.

## Comments
Nek5000 .rea parameters IOSTEP and IOTIME are set to zero during runtime, and instead IOPERT, initialised from the *params.input* file is used for saving perturbation snapshots. 

Note also that the iteration count in this module is slightly different, with a iteration as described in [1] corresponding to two folders, e.g. the direct and adjoint runs for the first iteration are stored in Iter00 and Iter01, respectively.

##Examples


Examples for the computation of modes for 2/3D channel flows are. To reproduce validation generate .map files for *channel_[2,3].rea*, modify *SIZE* to setup Nek5000 for a 2/3D run, compile and run *run* as
>  `makenek Nek_SSRM`
>  `./run channel_[2/3]D [number of cores] 0 20 [y/n]`



##Files
The main files, other than the inputs previously mentioned, are

* **runArnoldi** : A bash script that calls nek5000, to compute direct/adjoint runs, and octave, to compute the Arnoldi decomposition and inputs for the next run.
* **Nek_SSRM.ust** : Nek5000 user file to compute runs.
* **Nek_Arnoldi.m** : Octave/Matlab script that reads previous results,  computes an Arnoldi factorisation, saves the inputs for the next run, and if requested, writes estimation of the resolvent gains and modes to disk. 
* **readnek.m** and **writenek.m** : routines to load and write files for nek5000 in matlab/octave. These were obtained from [another project](https://github.com/nfabbiane/nekmatlab).



# References

1. Martini, Eduardo, Daniel Rodríguez, Aaron Towne, and André V. G. Cavalieri. “Efficient Computation of Global Resolvent Modes,” 2020. [https://arxiv.org/abs/2008.10904](https://arxiv.org/abs/2008.10904) 



