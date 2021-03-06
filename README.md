# Twitter Spread R-package

R-package with the code used for the paper: _Modelling rapid online cultural transmission: Evaluating neutral models on Twitter data with Approximate Bayesian Computation_ in _Palgrave Communication_ special issue on Cultural Evolution. Vignettes are availalbe online:
- [spreadrt.html](https://simoncarrignon.github.io/spreadrt/vignettes/spreadrt.html)
- [abc_spreadrt.html](https://simoncarrignon.github.io/spreadrt/vignettes/abc_spreadrt.html)

## Install the package

To install the 'spreadrt' package directly from this git repository you can use different option:
1. Clone the repository 
```bash
git clone https://github.com/simoncarrignon/spreadrt.git
```
and install from source:
```R
install.packages("spreadrt",repos=NULL,type="source")
library(spreadrt)
```
2. Or you can use `install_github` function from package `devtools`:
```R
library(devtools)
devtools::install_github("simoncarrignon/spreadrt")
```

## Start using the package
Once you installed the package, the best way to start learning and using the package is to use the vignettes available in `vignettes/`. You can read directly the `.Rmd` files and run directly the chunk of code from here or you can compile the code using the library
```R
library(devtools)
```
and the function:
```R
build_vignettes()
```

__warning__: Building the vignettes takes time (more than 30 minutes) as it has to run a great number of simulations and generate numerous graphs.


Once it's done, you can open the resulting html files that should have been created in `inst/doc/` at the follwing URLs: [spreadrt.html](inst/doc/spreadrt.html) and  [abc_spreadrt.html](inst/doc/abc_spreadrt.html).

A simpler and more straightforward way is to go to the compiled version of the vignettes in the pages branch of the github: 
- [spreadrt.html](https://simoncarrignon.github.io/spreadrt/vignettes/spreadrt.html)
- [abc_spreadrt.html](https://simoncarrignon.github.io/spreadrt/vignettes/abc_spreadrt.html)

## Modify the package
If you plan to modify the package and want to check your modifications you need to clone the repository:
```bash
git clone https://github.com/simoncarrignon/spreadrt.git

Then start `R` from the git folder and use the packages `devtools` and `roxygen2`
```R
library(devtools)
library(roxygen2)
```

then if you modify some code within the package you use: 
* `load_all()` to update the code  
* `documentation()` to update the documentation


## Approximate Bayesian Computation and `exec` subfolder

Some of the scripts used to run the ABC have been stored in `exec/abcdir` as well as some script used to generate the plot of the paper in `exec/palcomm`. Though those script are very hardware and problem specific they could be of some use (for people using SLURM and MPI mainly). 

__TODO__: Part of script and comment in `exec` should be integrated in functions and vignette in the main package.
