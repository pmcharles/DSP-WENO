# DSP-WENO Reconstruction and entropy stable flow solver
### Author: Phil Charles and Deep Ray, University of Maryland, College Park
### Webpage: pmcharles.github.io, deepray.github.io
### Emails: charlesp@umd.edu, deepray@umd.edu
### Date: 12 May 2025

**DSP-WENO** is a sign-preserving weighted essentially non-oscillatory (WENO) reconstruction method which incorporates an artificial neural network. This repository contains the code for the entropy stable TeCNO solvers for multi-dimensional scalar and systems of conservation laws. Additionally, the trainer is provided to train your own **DSP-WENO** if desired. Details on TeCNO solvers can be found in the following paper: [TeCNO](https://epubs.siam.org/doi/abs/10.1137/110836961?casa_token=YutiwMi6rcwAAAAA:jPQqHct6_kLnfO3yJudb6MMauXY1ENRPgrMqU0B1R0QDl7Lq7SRhYCzcE7y8kxmMs7qxTAuTLN4).

## Table of contents 

* [Compiling the code](#compiling-the-code)
* [Writing a param.in file](#writing-a-paramin-file)
    * [Scalar case](#scalar-case)
    * [Systems case](#systems-case)
* [Extracting 1D Line Data](#extracting-1d-line-data)
* [Training DSP-WENO](#training-dsp-weno)

## Compiling the code

The following instructions are to compile the systems code located in the enflo directory.

Set the following in your .bashrc file (or .profile file in unix or.zshrc if using zsh)

```
export ENFLO_HOME = <path to enflo>
```

e.g.,

```
export ENFLO_HOME = /home/deep/enflo
```

Set the compiler in ```$ENFLO_HOME/src/makefile.in``` and type ```make``` inside src directory to compile the code.

You should also add following lines to your .bashrc file

```
PATH=$PATH:$ENFLO_HOME/src 
PATH=$PATH:$ENFLO_HOME/utils 
PATH=$PATH:$ENFLO_HOME/extern/delaundo/std 
export PATH
```

Assuming the code compiles properly, test out the solver by going to ```$ENFLO_HOME/examples/1D/ModSod``` and running the following command

```
enflo -i param.in
```

The procedure is similar for the scalar conservation laws solver.

## Writing a param.in file

There are separate solvers for the scalar and systems cases, so the param.in files will necessarily be different. A solver can handle problems of one, two, or three dimensions. Example param.in files are available in the corresponding examples folders within each solver, and in this section, we detail what the general components of a file are.

### Scalar case

Here is an example of a param.in file for the scalar solver. It corresponds to a complicated Burgers' problem:

~~~matlab
grid
{
   x_min   -pi
   x_max   pi
   Nx_cell { 100,200,400,600,800,1000 }
   y_min   0
   y_max   1
   Ny_cell { 2,2,2,2,2,2 }
   z_min   0
   z_max   1
   Nz_cell { 2,2,2,2,2,2 }
}

numeric
{
   time_scheme  ssprk3
   cfl          0.5
   max_iter     50000
   final_time   0.5
   min_residue  1.0e-6
   reconstruct  sp_weno_dl
   net_file 	../../networks/dsp_weno.txt
}

material
{
   model       linadv
   vel         { 1 0 0 }
   flux        lin_tecno4
}

constants
{
}

initial_condition
{
   pow(sin(x),4)
}

exact_soln
{
   available   yes
   pow(sin(x-t),4)
}


boundary_condition
{
   xleft
   {
      type periodic
   }
   xright
   {
      type periodic
   }
   yleft
   {
      type periodic
   }
   yright
   {
      type periodic
   }
   zleft
   {
      type periodic
   }
   zright
   {
      type periodic
   }
}

output
{
   OOC_study        yes
   problem_dim      1
   time_stamps      5
   output_path .
}
~~~

The grid parameter specifies the dimensions of the problem. Note that three dimensions must always be specified. For one-dimensional or two-dimensional problems, simply set the number of cells of the unused dimensions to two. The cell sizes are lists of sizes to use, so you must ensure that each list is of the same size. If you are not doing an order of convergence study, simply use a single mesh size for each dimension.

The numeric parameter specifies the numerical details of the solver, such as the time-marching scheme, CFL, final time, and the reoncstruction method used in the TeCNO framework. Options include sp_weno, sp_wenoc, sp_weno_dl (**DSP-WENO**), as well as eno k (where k is the order of ENO to use). The net_file parameter specifies the location of the trained network's file.

> **_NOTE:_** If you are using the DSP-WENO reconstruction, ensure that the net_file parameter points to an existing weight/bias parameter file for the trained network. If such a file does not exist, the solver will not run.

The material parameter specifies the PDE model to use (linadv or burgers) as well as the flux. Here, lin_tecno4 means a fourth-order accurate TeCNO for linear advection will be used. If solving a linear advection problem, you must include the vel parameter. If solving a Burgers' equation problem, omit the vel parameter.

The constants parameter defines constants that may be used in the rest of the file when describing the problem.

The initial_condition parameter specifies the initial condition of the problem.

The exact_soln parameter firstly specifies if the exact solution is available and if so, specifies the exact solution.

The boundary_condition parameter specifies the type of boundary conditions to use (periodic, dirichlet, neumann, or wall). For one-dimensional or two-dimensional problems, simply set the conditions of the unused boundaries to periodic.

The output parameter specifies other details about the problem and output format. If OOC_study is yes, perform an order of convergence mesh refinement study, using the cell sizes specified in the grid parameter. The dimension of the problem is specified by problem_dim. The time_stamps parameter specifies how many times at which to save the solution. The output_path parameter specifies the location of the output; if not specified, creates an output directory in the working directory.

### Systems case

Here is an example of a param.in file for the systems solver. It corresponds to the Shu-Osher problem:

~~~matlab
grid
{
   x_min   -5
   x_max   5
   Nx_cell { 400 }
   y_min   0
   y_max   1
   Ny_cell { 2 }
   z_min   0
   z_max   1
   Nz_cell { 2  }
}

numeric
{
   time_scheme  ssprk3
   cfl          0.4
   max_iter     50000
   final_time   1.8
   min_residue  1.0e-6
   reconstruct  sp_weno_dl
   net_file     ../../../networks/dsp_weno.txt
}

material
{
   gamma       1.4
   gas_const   1.0
   viscosity
   {
      model    constant
      mu_ref   0.0
   }
   prandtl     1.0
   model       euler
   flux        kepes_tecno4_roe
   EC1_fix     no
}

constants
{
}

initial_condition
{
   density           3.857143*(x<-4) + (1+0.2*sin(5*x))*(x>=-4)
   xvelocity         2.629369*(x<-4)
   yvelocity         0.0
   zvelocity         0.0
   pressure          10.33333*(x<-4) + 1*(x>=-4)
}

exact_soln
{
   available   no
}


boundary_condition
{
   xleft
   {
      type neumann
   }
   xright
   {
      type neumann
   }
   yleft
   {
      type periodic
   }
   yright
   {
      type periodic
   }
   zleft
   {
      type periodic
   }
   zright
   {
      type periodic
   }
}

output
{
   OOC_study        no
   problem_dim      1
   time_stamps      5
   output_path .
}
~~~

Generally, this file follows the format of a scalar parameter file with some differences. Notably, the material parameter requires the specification of more information about the problem (e.g. gamma, gas_const, prandtl, etc.). Additionally, the initial_condition requires the specification of initial conditions for all variables. For Euler equations, for example, we have density, pressure, and all components of the velocity. Otherwise, the format is identical to the scalar case.

## Extracting 1D line data

When solving a 1D system of conservation laws, you may want to directly see the density, velocity, and pressure profiles (of course, you can do the same thing to look at a 1D scalar problem solution as well). This can easily be done with the line_extract.py script in the utils folder by runnnig the following command from the solution file directory:

```
python line_extract.py x1 x2 Nx
```

where x1 and x2 are the 1D boundary points and Nx is the number of cells (this command may vary depending on where line_extract.py is located, of course).

For 2D problems, we recommend examining the generated vtk files with a software such as Visit.

## Training DSP-WENO

In this repo, a trained DSP-WENO is present in the networks folders of the scalar and systems solver. If you want to train your own DSP-WENO, the python script and data are provided as well. Note that due to the randomness of training, not every network will result in the improved shock-capturing performance, but the trained network we have provided is one such example of a good network.

In order to train a network, first modify the config.py file to tune any hyperparameters. The data file to use is DSP-WENO_training_data.dat. Then, simply run DSP-WENO_network_trainer and navigate to the DSP-WENO_trained_models directory to find a directory containing the trained model and other information about its training run. To use with the entropy stable solver, use the generated dsp_weno.txt file.