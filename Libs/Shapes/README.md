# Shapes 
**Nek5000 ResolventTools :  Shapes** is used to obtain parametric spatial supports, which can be used, for example, in the definition of sensors and actuators. 

The code has as input a single file *shapes.input*, which contains one line for each shape to be used. For each shape a line containing
> `[G\g][G\g][G\g] x0 y0 z0 Sx Sy Sz sx sy sz`

The first three letters indicate a Gaussian support for the shape, which can be normalised to have maximum value of 1 (G) are area equal to 1 (g), with each letter corresponding to one of the three spatial coordinated. The following three parameters (x0,y0,z0) indicate the center of the Gaussian, followed by the half-width on each of these directons (Sx,Sy,Sz) and the flow variable on which it acts (sx,sy,sz).  For example
> `GgG 1 2 3 .1 .25 .35 1 .5 0`

can be used as a Gaussian-distributed sensor, located at the coordinate (1,2,3), with widths given by (.1,.25,.35), that reads velocity components on the x direction, with weight 1, and on the y direction, with weight .5 . 

The module has routines to recover the shapes, obtain readings (scalar products with a field) and to write the shapes to disk.

Other shapes, as well as spatial support provided in .fld files, still need to be implemented.