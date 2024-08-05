# Chi-Square_Function

This repository contains the Python implementation of the Chi-Square function designed to numerically validate the 4-zeroes texture model. This function utilizes four parameters: $A_{u}$, $A_{d}$, $\phi_{1}$, and $\phi_{2}$ to compute the Chi-Square criterion.

## Installation

To use the Chi-Square function, you need to have Python installed on your machine (the code was tested on Python 3.11). In order to use or test this function you will need the following packages:

```bash
numpy
cmath
```

## Function Overview
The Chi-Square function is defined as follows:
    
```bash

    def Chi_Square(Au, Ad, phi1, phi2):
    '''
    Calculate the Chi-Square statistic for the 4-zeroes texture model.

    Parameters:
    - Au: Free parameter Au
    - Ad: Free parameter Ad
    - phi1: Free parameter phi1
    - phi2: Free parameter phi2

    Returns:
    - chi_square: Calculated Chi-Square criterion
    '''
```
For the case study I, where $\eta_{u} = +1$ and $\eta_{d} = +1$, the the search space is defined by the following ranges for the parameters $A_{u}$, $A_{d}$, $\phi_{1}$, and $\phi_{2}$:

- $m_{c}< A_{u}< m_{t}$
- $m_{s}< A_{d}< m_{b}$
- $0<\phi_{1}<2\pi$
- $0<\phi_{2}<2\pi$

 The experimental values for $m_{c}$, $m_{t}$, $m_{s}$, $m_{b}$, along with their respective uncercanities, are already defined within the file.

 ## Usage

 To use the Chi-Square function, simply call it with the required parameters. Here's a basic example:

```bash
criterion = Chi_Square(139297.25173570923,2393.291447697453,1.6830151602082093,6.1851624875199525)
print(criterion)
#output
0.7904494037073188
```
# License

This project is licensed under the GPL-3.0 License. 
