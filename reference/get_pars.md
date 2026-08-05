# Get TMB parameter list

Get TMB parameter list

## Usage

``` r
get_pars(object)
```

## Arguments

- object:

  Fit from
  [`sdmTMB()`](https://sdmTMB.github.io/sdmTMB/reference/sdmTMB.md).

## Value

A named list of parameter values.

## Examples

``` r
fit <- sdmTMB(present ~ 1, data = pcod_2011, family = binomial(), spatial = "off")
pars <- get_pars(fit)
names(pars)
#>  [1] "ln_H_input"          "b_j"                 "b_j2"               
#>  [4] "bs"                  "ln_tau_O"            "ln_tau_Z"           
#>  [7] "ln_tau_E"            "ln_kappa"            "log_kappaS_nl"      
#> [10] "kappaT_nl_raw"       "thetaf"              "ln_student_df"      
#> [13] "gengamma_Q"          "logit_p_extreme"     "log_ratio_mix"      
#> [16] "ln_phi"              "ln_tau_V"            "rho_time_unscaled"  
#> [19] "ar1_phi"             "logit_rho_sar"       "re_cov_pars"        
#> [22] "re_b_pars"           "b_rw_t"              "omega_s"            
#> [25] "zeta_s"              "epsilon_st"          "b_threshold"        
#> [28] "b_epsilon"           "ln_epsilon_re_sigma" "epsilon_re"         
#> [31] "b_smooth"            "ln_smooth_sigma"    
```
