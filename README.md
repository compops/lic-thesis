lic-thesis
====

Sequential Monte Carlo for inference in nonlinear state space models

This code was downloaded from < https://github.com/compops/lic-thesis > or from < http://liu.johandahlin.com/ > and contains the code used to produce the results in the Licentiate's thesis

* J. Dahlin, *Sequential Monte Carlo for inference in nonlinear state space models*. Linköping Studies in Science and Technology. Thesis. No. 1652, April 2014.

The thesis is available from < http://liu.johandahlin.com/ >. Put all the code into the same folder with the data directory as an intact subdirectory and execute first the Python code and then the R code to generate the plots. Some results might differ slightly from the thesis due to random differences in the data or the algorithms.

Requirements
--------------
The main programs are written in Python 2.7 and makes use of NumPy 1.7.1, SciPy 0.12.0, Matplotlib 1.2.1, Pandas 0.13.1 and GPy 0.4.9. Please have these packages installed, on Ubuntu they can be installed using "sudo pip install --upgrade *package-name* ". The plotting is done in R 3.0.1 and does not require any special packages. If any are missing from our system, they can be installed by executing the command "install.packages("packagename")" in the R console.

Included files (examplesForThesis-python and examplesForThesis-R)
--------------
**ch2-example-likelihoodtheory**
Recreates the plot in "Score and information matrix for the LGSS model" in Example 2.6 in Section 2.3.

**ch3-example-earthquakefilteringsmoothing**
Recreates the plot in "State inference in the earthquake count model" in Example 3.2 in Section 3.3.2.

**ch3-example-forgettingproperties**
Recreates the plots in "Mixing property in the LGSS model" in Example 3.5 in Section 3.4.1.

**ch3-example-garchpathdegenercy**
Recreates the plots in "Path degeneracy in the GARCH(1,1) model" in Example 3.3 in Section 3.3.2.

**ch3-example-hwsvfilteringsmoothing**
Recreates the plot in "State inference in the Hull-White SV model in Example 3.6 in Section 3.4.1.

**ch3-example-hwsvllscoreinfo**
Recreates the plot in "Score and information matrix in the Hull-White SV model" in Example 3.7 in Section 3.4.2.

**ch3-example-importancesampling-hwsv**
Recreates the plot in "IS for Bayesian parameter inference in the HWSV model" in Example 3.1 in Section 3.2.

**ch3-example-llestimationCLT**
Recreates the plot in "Bias and varaiance of the log-likelihood estimate" in Example 3.4 in Section 3.3.4.

**ch4-example-earthinference**
Recreates the plot in "GPO for ML inference in the earthquake count model" in Example 4.9 in Section 4.4.3.

**ch4-example-garchinference**
Recreates the plot in "PMH0 for parameter inference in the GARCH(1,1) model" in Example 4.4 in Section 4.3. This code takes some time to run (in the order of hours).

**ch4-example-gp-acqfunc**
Recreates the plot in "GPO using different acqusition rules" in Example 4.7 in Section 4.4.2.

**ch4-example-gpo-garchinference**
Recreates the plot in "GPO for ML inference in the GARCH(1,1) model" in Example 4.8 in Section 4.4.3.

**ch4-example-gp-priorrealisations**
Recreates the plot in "GP kernels" in Example 4.5 in Section 4.4.1.

**ch4-example-gp-regression**
Recreates the plot in "GP regression" in Example 4.6 in Section 4.4.1.

**ch4-example-lgssinference-mixing**
Recreates the plot in "Parameter inference in the LGSS model" in Example 4.1 in Section 4.2.

**ch4-example-lgssinference**
Recreates the plot in "Parameter inference in the LGSS model" in Example 4.1 in Section 4.2.

**ch4-example-lgss-inputdesign**
Recreates the plot in "Input design in the LGSS model using GPO" in Example 4.10 in Section 4.4.3.

Supporting files (helpers-python)
--------------
**pmh.py**
Defines the general class for the particle MH algorithm and helper functions for this.

**smc.py**
Defines the general class for sequential Monte Carlo algorithm.

**kf.py**
Defines the general class for Kalman filtering and smoothing algorithm.

**gpohelpers.py**
Defines the acqusition rules for the GPO algorithm.

**classes.py**
Defines the different system models and generates the data.

**helpers.py**
Defines different helpers for the other functions.
