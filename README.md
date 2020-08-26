# Double-Slit-Experiment
A MATLAB script for modeling and animating the double slit experiment using the Finite Element Method and Euler averages 

## Problem Statement
We use FEM to numerically solve for u in the wave equation given as:

![](./wave_equation.jpg)
 
And are given the initial and boundary conditions: 

![](./boundary_equation.jpg)

Where Γ_g, Γ_h, Ω, I  are the Neumann boundary, Dirichlet boundary, physical domain, and temporal domain, respectively. The Dirichlet function g(t) is given as 0.5*sin(6πt), and models the wavelength of light coming into the double slits. The Neumann function p is given as 0 along the Neumann boundary, which is the boundaries of our 2-D box (this represents the approximation that no light will travel through the box).  

# Code Functionality

```Double_Slit_Experiment.m``` is a MATLAB implementation of the solution that proivdes functionality to explore Euler averaging averaging methods and mass lumping. In addition to solving the problem, the code outputs an animated video of the solution. 
![](./boundary_equation.jpg)
# Calculations
The  ``` Double_Slit_Project_Writeup.pdf ``` contains the detailed calculations and examples of Euler averaging. 
