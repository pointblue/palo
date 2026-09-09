
<!-- README.md is generated from README.Rmd. Please edit that file -->

# palo

<!-- badges: start -->

<!-- badges: end -->

This R package contains useful functions for analyzing data collected at
the Palomarin Field Station and closely-associated monitoring programs.

## Installation

You can install the development version of palo from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("pointblue/palo")
```

Fitting the BBS-style hierarchical models requires also installing
[JAGS](https://sourceforge.net/projects/mcmc-jags/files/).

## Example

This is a basic example which shows you how to solve a common problem of
compiling and summarizing multiple point count data files downloaded
from the [AKN](https://avianknowledge.net).

``` r
library(palo)

# combine one or more PC data files:
data_compiled = compile_PC_dat(dir = 'C:/directory_containing_PC_data', 
                               pattern = '.csv')

# convert distance bins to numeric mindist and maxdist fields and limit the maxdist
data_filtered = data_compiled |> 
  add_distance_key() |> 
  filter(maxdist <= 50)

# subset to one species from one project, and summarize the total count per station and visit (regardless of detection distance)
SOSP_data = summarize_PC_dat(df = data_filtered,
                             species = 'SOSP',
                             project = 'PINN')
```

<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->

<!-- ```{r cars} -->

<!-- summary(cars) -->

<!-- ``` -->

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/v1/examples>. -->

<!-- You can also embed plots, for example: -->

<!-- ```{r pressure, echo = FALSE} -->

<!-- plot(pressure) -->

<!-- ``` -->

<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
