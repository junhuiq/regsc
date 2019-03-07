# What is regsc?
**regsc** is an R package for estimating linear regressions with structural changes using GFL (group-fused Lasso), which is proposed in Qian and Su (2016). We consider the following regression model with possibly multiple structural changes in <br>
y<sub>t</sub> = x<sub>t</sub> ' &beta;<sub>t</sub>  + z<sub>t</sub> ' &gamma; + u<sub>t</sub>,<br>
where x<sub>t</sub> is a p-by-1 vector of variables that may have time-varying effects on y, z<sub>t</sub> is a q-by-1 vector of variables that have constant effects on y. We assume that although the number of structural changes in &beta;<sub>t</sub> is unknown, it is much smaller than the sample size.

# How to install regsc?
## From GitHub
``` r
install.packages("devtools")
library(devtools)
install_github("junhuiq/regsc")
```
# An Example
After installing **regsc**, one can try estimating a monetary policy rule of the Federal Reserve:
``` r
library(regsc)
res = regsc(fedrate~inflation+gap,policyrule)
summary(res)
```

# References
Junhui Qian, Liangjun Su, 2016, Shrinkage estimation of regression models with multiple structural changes, Econometric Theory, 32 (6), 1376-1433. 
