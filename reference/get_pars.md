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
#>  [1] "ln_H_input"        "b_j"               "b_j2"             
#>  [4] "b_disp_k"          "bs"                "ln_tau_O"         
#>  [7] "ln_tau_Z"          "ln_tau_E"          "ln_kappa"         
#> [10] "log_kappaS_nl"     "log_kappaT_nl"     "thetaf"           
#> [13] "ln_student_df"     "gengamma_Q"        "psi"              
#> [16] "logit_p_extreme"   "log_ratio_mix"     "ln_phi"           
#> [19] "ln_tau_V"          "rho_time_unscaled" "ar1_phi"          
#> [22] "logit_rho_sar"     "re_cov_pars"       "re_b_pars"        
#> [25] "b_rw_t"            "omega_s"           "zeta_s"           
#> [28] "epsilon_st"        "b_threshold"       "b_epsilon"        
#> [31] "b_smooth"          "ln_smooth_sigma"  
```
