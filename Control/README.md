# Nek5000 Resolvent Tools : Control 
**Nek5000 ResolventTools :  Control** provides an implementation of the Wiener-Hopf regulator described in [1]. The estimation and control Kernels are obtained on a matrix-free approach from the post processing of readings from time marching of the adjoint and direct linearised equations.

The code makes uses the **Shapes** module for sensor/actuator/target spatial supports.  Please see the README file on /Libs/Shapes for a list of available shapes. Each shape is parametrised on the file *shapes.input*, one shape per line. Any of such shapes can be used to provide the desired spatial support. The code assumes that sensor and actuators used are sequential shapes, e.g. sensors corresponding to shapes 29 to 31 are used as actuators, so please organize the *shapes.input* file accordingly.



##Setup

Copy in the folder relevant .rea/.map files for the flow of interest. For the linearised runs, set in the .rea file 

1. a constant time step (**DT/param(12)**<0)
2.  **npert** to -1
3.  a restart option to load the flow for which the analysis is to be performed. The restart file should have step and time equals to zero.


Setup the *params.input* file, which has 5 lines, namely

1. **runType**:  +-1 for direct/adjoint run using a sensor as initial condition ; +-2 for a direct/adjoint run using a previous solution as a full rank external forcing ; +-2 for a direct/adjoint run using a previous solution as a low rank external forcing ; 4/5 to run the uncontrolled/controlled linear system disturbed by white-noise forcing
2. **shapeRange1\_start shapeRange1\_end** : range of shapes 1
3. **shapeRange2\_start shapeRange2\_end** : range of shapes 2
4. **tolerance** : tolerance used to terminate dir/adj runs. Runs will be terminated when the norm is smaller than the largest previous norm times the tolerance used.

The shape ranges have different uses: when using low rank forces on a run, range 1 specifies the shapes used as inputs; for control runs shape ranges 1 and 2 specify the sensors and actuators that will be used for control.

The **run** script automates most of the setup necessary to perform dir/adj/control runs. Run it without arguments to check its syntax. The **runExamples** script run all the necessary data for construction of the kernels presented in [1].

After running the necessary dir/adj runs, use the script *ComputeKernel.m*. to obtain the optimal control kernels. The script also provides plots, and if low rank targets are used, estimations of the control performance. It assumes that the adj/dir runs are stored in the following format `./Runs_[ishape]_[runType]_full`, where *runType*. can be *dir*, *dir-adj*, *dir-adj-dir*, *adj*, *adj-dir* or *adj-dir-adj*.

## Construction estimation and control Kernels
The folder *PostProcessing* contains Matlab scripts to control and estimation Kernels. The main script in *ComputeKernel.m*, which reads sensor readings from the dir-* /adj-* runs and performs all computations to obtain control and estimation Kernels. In all the scripts used the data is assumed to be uniformly spaced in time (constant sampling rate), and all the matrices constructed as M(i,j,it), e.g. the index related to the time dependence being the third one.

The scripts are described bellow : 

| File | Description   |
| :---:   | :-: |
| ComputeKernel.m | Octave/Matlab script to obtain estimation and control kernels. Control Kernels are saved to disk for later use in control runs. Variables "iy,ia,iz" indicate which shapes will be used as sensors/actuators/targets. |
| readData1.m(readData2.m) | Loads readings from the linearised runs that will be used for construction of the control laws. The file reads data for full-rank forces and low-(full-) rank targets.  | 
| interpTensors.m | Interpolate tensors/matrices from one time discretisation into another. | 
| multiprod.m | Library for fast multiplication of several matrices. | 
| WienerHopfDec.m | Performs Wiener-Hopf matrix multiplicative factorization. | 
| getHGs.m | Uses data from readData?.m to construct the relevant Hl,Hr,Gl,Gr,HrCzGr terms. | 
| getEstimationKernels.m | Uses data from getHGs to construct estimation kernels (causal, non-causal, and truncated non-causal). | 
| getControlKernels.m | Uses data from getHGs to construct control kernels (causal, non-causal, and truncated non-causal). | 
| chilbert.m | Generalizes Matlab's *hilbert* function to compute Hilbert transforms of complex data. | 


###Main files description
| File | Description   |
| :---:   | :-: |
| Nek\_EstCon.usr   | Nek5000 user file to interface with the module. Can be used to run the necessary dir/adj systems, and  for control of the linear and non-linear problems. |
| ComputeKernel.m | Octave/Matlab script to obtain estimation and control kernels. Control Kernels are saved to disk for later use in control runs. |
| run | Script that facilitates the running of adj/dir runs. Run without commands to view the syntax to be used|
| runExamples | Available on the examples branch, run all the necessary data for construction of the kernels presented in [1]|
| B.fld (not yet implemented) | A .fld file, following Nek5000 standards, containing the diagonal entries of a high rank **B** matrix, used only for high rank forces. If not present, **B=I** is assumed. |
| C.fld (not yet implemented) | A .fld file, following Nek5000 standards, containing the diagonal entries of a high rank **C** matrix, used only for high targets forces. If not present, **C=I** is assumed. |
 

##Comments
Please not that the code has been written as to be an efficient manner to obtain control and estimation kernels. As currently implemented, the use of the control kernel for flow control introduces significant overhead, and is not parallelised yet, meaning that it easily be the run's bottleneck if a large number of processors is used. 

A more efficient implementation for the use of the control Kernel will be developed and later added to the repository. 

To apply white noise disturbances the code uses a random number generator which is initialised on the first iteration based on the current processor ID. In order to compare different control strategies use the same number of processors to obtain the same random time series on each of the runs.

# Roadmap
The following features will be developed soon

* Parallelisation of the control module, and use of a higher order interpolation scheme to reduce the number of points needed in the control law.
* Inclusion of B.fld and C.fld files to easily implement diagonal B and C matrices.
* Generalize for the use of low-rank forces.

# References
1. E. Martini, J. Jung ,  A. V. G. Cavalieri, P. Jordan, L. Lesshafft & A. Towne, Optimal causal Resolvent-based control using the Wiener-Hopf formalism. Journal of Fluid Mechanics, (2021) [submitted]


