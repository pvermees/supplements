# Summary

This repository contains all the data and code needed to reproduce the
figures of the GChron manuscript entitled
"Too good to be true: underdispersion in geochronology" (in review).

## Repository structure

This repository contains five folders:

1. [`code`](code): `R` and `Julia` code to generate tables and figures

2. [`data`](data): `.csv`-files with data used by `code`

3. [`DLNR`](DLNR): `.csv`-files with Agilent 7900 data for a detrital
zircon U-Pb age spectrum used in the paper

4. [`figures`](figures): the output of the `R` scripts in `code`

5. [`paper`](paper): `LaTeX` source code for the GChron manuscript

## Prerequisites for `code`

1. [`R`](https://r-project.org): a free statistical programming language

2. [`Julia`](https://julialang.org): a free scientific programming language

3. [`IsoplotR`](https://cran.r-project.org/package=IsoplotR): an `R`-package

4. [`KJ`](https://github.com/pvermees/KJ.jl): a `Julia`-package for LA-ICP-MS data reduction

5. [`LaTeX`](https://latex.org): a markup language

## Installation instructions for `code`

1. Install `R` and `Julia` on your system

2. Start `R` and install `IsoplotR` from CRAN by entering the following command at the console:

```
install.packages('IsoplotR')
```

3. Install the `KJ`, `CSV`, `Optim` and `Distributions` packages from the `Julia` REPL:

```
using Pkg; Pkg.install("KJ","CSV","Optim","Distributions")
```

## Instructions to reproduce the figures

1. To generate the p-value distributions of the detrital zircon U-Pb
data, start `Julia` and run the following script at the REPL:

```
include("Droellner.jl")
```

This will create `DLNR.csv` and `DLNRcherries.csv` in [`data`](data)

2. To generate the figures, run the following script in `R`:

```
source("synthetic.R")
source("real.R")
```

This will create 11 `.pdf` files in [`figures`](figures)

3. Render the LaTeX code in [`paper`](paper) in a terminal window:

```
pdflatex paper.tex
bibtex paper.tex
pdflatex paper.tex
pdflatex paper.tex
```

This will create `paper.pdf`, which contains the GChron manuscript.
