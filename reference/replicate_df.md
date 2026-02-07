# Replicate a prediction data frame over time

Useful for replicating prediction grids across time slices used in model
fitting.

## Usage

``` r
replicate_df(dat, time_name, time_values)
```

## Arguments

- dat:

  Data frame.

- time_name:

  Name of time column in output.

- time_values:

  Time values to replicate `dat` over.

## Value

A data frame replicated over `time_values` with a new column based on
`time_name`.

## Examples

``` r
df <- data.frame(variable = c("a", "b"))
replicate_df(df, time_name = "year", time_values = 1:3)
#>   variable year
#> 1        a    1
#> 2        b    1
#> 3        a    2
#> 4        b    2
#> 5        a    3
#> 6        b    3

head(qcs_grid)
#>     X    Y    depth depth_scaled depth_scaled2
#> 1 456 5636 347.0834    1.5608122    2.43613479
#> 2 458 5636 223.3348    0.5697699    0.32463771
#> 3 460 5636 203.7408    0.3633693    0.13203724
#> 4 462 5636 183.2987    0.1257046    0.01580166
#> 5 464 5636 182.9998    0.1220368    0.01489297
#> 6 466 5636 186.3892    0.1632882    0.02666303
nd <- replicate_df(qcs_grid, "year", unique(pcod$year))
head(nd)
#>     X    Y    depth depth_scaled depth_scaled2 year
#> 1 456 5636 347.0834    1.5608122    2.43613479 2003
#> 2 458 5636 223.3348    0.5697699    0.32463771 2003
#> 3 460 5636 203.7408    0.3633693    0.13203724 2003
#> 4 462 5636 183.2987    0.1257046    0.01580166 2003
#> 5 464 5636 182.9998    0.1220368    0.01489297 2003
#> 6 466 5636 186.3892    0.1632882    0.02666303 2003
table(nd$year)
#> 
#> 2003 2004 2005 2007 2009 2011 2013 2015 2017 
#> 7314 7314 7314 7314 7314 7314 7314 7314 7314 
```
