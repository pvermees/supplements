# Summary

This repository contains all the `R` code needed to reproduce the
figures of "Carbonate <sup>206</sup>Pb/<sup>238</sup>U problems and
potential <sup>207</sup>Pb/<sup>235</sup>U fixes" by Vermeesch, P.,
McLean, N., Vaks, A., Golan, T., Breitenbach, S.F.M and Parrish,
R. Geochronology (in review).

## Prerequisites

1. [`R`](https://r-project.org): a free statistical programming language

2. [`IsoplotR`](https://cran.r-project.org/package=IsoplotR): an `R`-package

## Installation instructions

1. Install `R` on your system

2. Start `R` and install `IsoplotR` from CRAN by entering the following command at the console:

```
install.packages('IsoplotR')
```

## Running the code

1. Start `R`

2. Navigate to the folder containing `main.R`. You can either do this
by clicking the `set working directory` button in your `R` IDE, or by
running`setwd(/path/to/main.R)` at the console.

3. Run the following command at the console:

```
source('main.R')
```

4. Inspect the `output` folder for the results.