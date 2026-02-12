# Get TMB parameter list

Get TMB parameter list

## Usage

``` r
get_pars(object)
```

## Arguments

- object:

  Fit from
  [`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md)

## Value

A named list of parameter values

## Examples

``` r
fit <- sdmTMB(present ~ 1, data = pcod_2011, family = binomial(), spatial = "off")
pars <- get_pars(fit)
names(pars)
#>  [1] "ln_H_input"          "b_j"                 "b_j2"               
#>  [4] "bs"                  "ln_tau_O"            "ln_tau_Z"           
#>  [7] "ln_tau_E"            "ln_kappa"            "thetaf"             
#> [10] "ln_student_df"       "gengamma_Q"          "logit_p_extreme"    
#> [13] "log_ratio_mix"       "ln_phi"              "ln_tau_V"           
#> [16] "rho_time_unscaled"   "ar1_phi"             "re_cov_pars"        
#> [19] "re_b_pars"           "b_rw_t"              "omega_s"            
#> [22] "zeta_s"              "epsilon_st"          "b_threshold"        
#> [25] "b_epsilon"           "ln_epsilon_re_sigma" "epsilon_re"         
#> [28] "b_smooth"            "ln_smooth_sigma"    
```
