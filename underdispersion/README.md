# Summary

This repository contains all the `R` code needed to reproduce the
figures of "Too good to be true: underdispersion in geochronology".

## Contents

1. `main.R`: an `R` script that reproduces the figures in the paper.

2. [`geochron.org`](geochron.org): a folder with `.xlsx` files
containing all the fission track data contained in the
[geochron.org](https://geochron.org) database on 7 October 2025.

3. `parse.R`: an `R` script that extracts all the p-values from the
contents of the `geochron.org` folder and stores them in
`p-values_FT.csv`.

4. `output`: the output of `main.R`, i.e. the four figure panels in
the paper.

## Prerequisites

1. [`R`](https://r-project.org): a free statistical programming
language.

2. `readr`: a dependency of `parse.R`, to be installed from CRAN by
entering the following command at the `R` console:

```
install.packages('readr')
```

## Running the code

1. Start `R`

2. Navigate to the folder containing `main.R`. This can either be done
by clicking the `set working directory` button in the `R` IDE, or by
running `setwd(/path/to/main.R)` at the console.

3. Run the following command at the console:

```
source('main.R')
```

4. Inspect the `output` folder for the results.
