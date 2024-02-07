## Simulation based Robot Modeling, Path Planning and Control
 Series of simulations and testing based on the dynamics, control, inverse kinematics, path and trajectory optimization of a robot. are performed in MATLAB Application of RRT algorithm for various scenarios, evaluation of the inverse kinematics and implementation of control stability is performed.

 ### Steps to view the code:
* Have MATLAB R2019a pr R2019b + installed. Alternatively, one can also use MATLAB Online at: https://www.mathworks.com/products/matlab-online.html. Octave is a secondary alternative.
* Note: These files are MATLAB Live scripts, .mlx files instead of the standard .m files.
* Download the folder and extract all the .mlx files in the same directory (This step is important, as MATLAB identifies all the related functions and mentions of the functions by having the files in the same directory).
 
* **Best Alternate Approach**: View the project with results in my portfolio: https://www.nishantkb.info/robotics.html

### Relevant Information:

*  Important Support functions are in a separate folder provided for reference. They include functions for Forward Kinematics (*fk*), Space Jacobian (*JacS*), Skew-symmetric matrix converter for both 3x3 and 4x4 matrices (*bracket3* & *bracket4*), Adjoint calculation (*adjointM*) (different from the adjoint function provided in MATLAB), Coriolis matrix evaluation (*coriolis*) and R2AxisAngle evaluation (*r2axisangle*).
*  We have covered redundant robots (more than required DoF) as well and that also helps us figure out the singularity cases in a robot arm/structure.
