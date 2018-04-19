# ViSA: Visualization of Sensitivity analysis

**ViSA** is a R package that aims to apply the sensitivity analysis workflow to quantify and visualize the parameter sensitivity from [**multi-compartment models**](https://en.wikipedia.org/wiki/Multi-compartment_model). It includes some functions that can use to check the convergence as well as graphical checking the output from R [**sensitivity**](https://cran.r-project.org/web/packages/sensitivity/index.html) package. 

This package can be installed via the devtools package using:  
```
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("nanhung/ViSA")
```
