# Mitochondrial dynamics

A set of code to simulate the effect of mitochondrial dynamics on the bioenergetics of a cardiomyocyte.

## Requirements

- MATLAB (R2018a or higher)
- Parallel Computing Toolbox (version 6.12)
- MATLAB Distributed Computing Server (version 6.12)

## Usage
To run the model and save the output in a variable called `output` use
```Matlab
output = model;
```
Population values can can be accessed using dot notation. See <https://au.mathworks.com/help/matlab/matlab_prog/access-data-in-a-structure-array.html> for more information on using dot notation.

### Running code on non Linux systems
This code uses [MEX](https://www.mathworks.com/help/matlab/ref/mex.html) files for performance. To run this code on a non Linux system, replace `myo_odefun_mex` on line 34 of `main_myo.m` with `myo_odefun`. 
