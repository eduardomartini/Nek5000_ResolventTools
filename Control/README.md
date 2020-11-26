# Nek5000 Resolvent Tools : Control 
**Nek5000 ResolventTools :  Control** provides an implementation of the Wiener-Hopf regulator described in [1]. The estimation and control Kernels are obtained on a matrix-free approach from the post processing of readings from time marching of the adjoint and direct linearized equations.

The code makes uses of the **Shapes** module for sensor/actuator/target spatial suports.  Please see the README file on /Libs/Shapes for a list of available shapes. Each shape is parametrised on the file *shapes.input*, one shape per line. Any of such shapes can be used to provide the desired spatial support. The code assumes that sensor and actuators used are sequential shapes, e.g. sensors corresponding to shapes 29 to 31 are used as actuators, so please organize the shape files accordingly.



##Setup

Copy in the folder relevant .rea/.map files for the flow of interest. For the linearized runs, set in the .rea file 

1. a constant time step (**DT/param(12)**<0)
2.  **npert** to -1
3.  a restart option to load the flow for which the analysis is to be performed. The restart file should have step and time equals to zero.


Setup the *params.input* file, which has 5 lines, namely

1. **runType**:  +-1 for direct/adjoint run using a sensor as initial condition ; +-2 for a direct/adjoint run using a previous solution as a full rank external forcing ; +-2 for a direct/adjoint run using a previous solution as a low rank external forcing ; 4/5 to run the a uncontrolled/controled linear system disturbe by white-noise forcing
2. **shapeRange1\_start shapeRange1\_end** : range of shapes 1
3. **shapeRange2\_start shapeRange2\_end** : range of shapes é
4. **tolerance** : tolerance used to terminate dir/adj runs. Runs will be terminated one the norm is smaller than the largest previous norm times the tolerance used.

The shape ranges have different uses: when using low rank forces on a run, range 1 specifies the shapes used as inputs; for control runs shape ranges 1 and 2 specify the sensors and actuators that will be used for control.

The **run** script automates most of the setup necessary to perform dir/adj/control runs. Run it without arguments to check its syntax. The **runExamples** script run all the necessary data for construction of the kernels presented in [1].

After running the necessary dir/adj runs, use the script *ComputeKernel.m*. to obtain the optimal control kernels. The script also provides plots, and if low rank targets are used, estimations of the control performance. It assumes that the adj/dir runs are stored in the following format `./Runs_[ishape]_[runType]`, where *runType*. can be *dir*, *dir-adj*, *dir-adj-dir*, *adj*, *adj-dir* or *adj-dir-adj*.



###Main files description
| File | Description   |
| :---:   | :-: |
| Nek\_EstCon.usr   | File to be compiled for dir/adj runs, and for control of the linearized problem. |
| Nek\_EstCon\_NL.usr | File to be compiled control of the non-linear problem. | 
| ComputeKernel.m | Octave/Matlab script to obtain estimation and control kernels. Control Kernels are saved to disk for later use in control runs. |
| run | Script that facilitates the running of adj/dir runs. Run without commands to view the syntax to be used|
| runExamples | Available on the examples branch, run all the necessary data for construction of the kernels presented in [1]|
| B.fld | A .fld file, following Nek5000 standards, containing the diagonal entries of a high rank **B** matrix, used only for high rank forces. If not present, **B=I** is assumed. |
| C.fld | A .fld file, following Nek5000 standards, containing the diagonal entries of a high rank **C** matrix, used only for high targets forces. If not present, **C=I** is assumed. |
 

##Comments
Please not that the code has been written as to be an efficient manner to obtain control and estimation kernels. As currently implemented, the use of the control kernel for flow control introduces significant overhead, and is not parallel, meaning that it easily be the run's bottleneck if a large number of processors is used. 

A more efficient implementation for the use of the control Kernel will be developed and latter added to the repository. 

To apply white noise disturbances the code uses a random number generator which is initialised on the first iteration based on the current processor ID. In order to compare different control strategies use the same number of processors to obtain the same random time series on each of the runs.


# References

1. Martini, Eduardo, Daniel Rodríguez, Aaron Towne, and André V. G. Cavalieri. ...

