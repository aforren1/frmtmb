### Q1a. Wall time (s), 4 chains x 1000/1000 run one after another

| id | key | n | t_frm | t_codegen | t_stanli_setup | t_stanli_sample | t_rstan_compile | t_rstan_sample | speedup_sample | speedup_total |
|---|---|---|---|---|---|---|---|---|---|---|
| p003 | gaussian x \|ID\| | 200 | 0.38 | 0.16 | 0.07 | 0.69 | 58.09 | 5.44 | 7.9 | 83.6 |
| p041 | rescor x mm() | 200 | 0.36 | 0.19 | 0.09 | 4.19 | 63.92 | 9.51 | 2.3 | 17.2 |
| p060 | von_mises x \|ID\| | 200 | 0.66 | 0.21 | 0.08 | 4.49 | 56.2 | 33.02 | 7.4 | 19.5 |
| p062 | gaussian x smooth | 200 | 0.38 | 0.11 | 0.07 | 1.35 | 46.78 | 7.28 | 5.4 | 38.1 |
| p075 | cumulative x hypothesis_profile | 200 | 0.29 | 0.09 | 0.05 | 0.53 | 45.21 | 3.93 | 7.4 | 84.7 |
| p082 | gaussian x hsgp | 200 | 0.45 | 0.14 | 0.1 | 1.96 | 57.45 | 10.04 | 5.1 | 32.8 |
| p096 | negbinomial x fitted | 200 | 0.31 | 0.06 | 0.08 | 0.33 | 49.52 | 2.43 | 7.4 | 126.7 |
| p121 | poisson x double_bar | 200 | 0.35 | 0.11 | 0.06 | 1.08 | 55.8 | 15.32 | 14.2 | 62.4 |
| p152 | weibull x cens() | 200 | 0.33 | 0.08 | 0.06 | 0.59 | 52.91 | 1.84 | 3.1 | 84.2 |
| big_crossed | big_crossed | 6000 | 1.48 | 0.16 | 0.8 | 18.26 | 57.25 | 96.73 | 5.3 | 8.1 |
| big_pois | big_pois | 10000 | 1.01 | 0.13 | 0.82 | 106.19 | 62.81 | 941.22 | 8.9 | 9.4 |
| big_slope | big_slope | 4000 | 1.56 | 0.26 | 0.58 | 17.77 | 73.26 | 114.83 | 6.5 | 10.3 |

### Q1b. stanli against rstan: diagnostics and posterior means

z = (mean_stanli - mean_rstan) / sqrt(mcse_stanli^2 + mcse_rstan^2) over every variable both report (r_ included, lp__ included); crit = the two-sided 5% Bonferroni bound for that many z.

| id | stanli_div | rstan_div | stanli_maxrhat | rstan_maxrhat | stanli_miness | rstan_miness | sr_nvar | sr_maxz | sr_crit | sr_frac3 | sr_lp_z | stanli_error | rstan_error |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| p003 | 0 | 1 | 1.008 | 1.004 | 695.502 | 719.389 | 75 | 2.315 | 3.403 | 0 | 0.421 |  |  |
| p041 | 0 | 0 | 1.008 | 1.003 | 968.934 | 990.619 | 33 | 1.655 | 3.172 | 0 | -0.744 |  |  |
| p060 | 0 | 0 | 1.004 | 1.006 | 985.193 | 870.768 | 75 | 2.21 | 3.403 | 0 | 1.398 |  |  |
| p062 | 1 | 2 | 1.004 | 1.003 | 1238.807 | 1129.558 | 24 | 1.21 | 3.078 | 0 | 0.966 |  |  |
| p075 | 0 | 0 | 1.002 | 1.001 | 2100.029 | 1829.312 | 10 | 0.825 | 2.807 | 0 | -0.761 |  |  |
| p082 | 19 | 37 | 1.001 | 1.006 | 1260.199 | 934.513 | 17 | 2.789 | 2.974 | 0 | 0.64 |  |  |
| p096 | 0 | 0 | 1.001 | 1.002 | 2191.5 | 2086.615 | 7 | 1.443 | 2.69 | 0 | -0.296 |  |  |
| p121 | 0 | 1 | 1.008 | 1.007 | 767.813 | 914.674 | 47 | 1.905 | 3.273 | 0 | -0.162 |  |  |
| p152 | 0 | 0 | 1.001 | 1.002 | 1676.518 | 2182.451 | 7 | 1.765 | 2.69 | 0 | -1.765 |  |  |
| big_crossed | 0 | 0 | 1.007 | 1.006 | 772.739 | 810.651 | 687 | 2.95 | 3.967 | 0 | 0.182 |  |  |
| big_pois | 0 | 0 | 1.012 | 1.009 | 902.797 | 912.441 | 4006 | 4.135 | 4.369 | 0.002 | 0.532 |  |  |
| big_slope | 0 | 0 | 1.011 | 1.007 | 609.758 | 564.904 | 3013 | 3.456 | 4.306 | 0.001 | 1.541 |  |  |

### Q1c. frm() maximum likelihood against the stanli posterior (brms default priors)

z = (MLE - posterior median) / posterior sd. `mle_maxz_hyp` over b_, sd_, cor_, sigma and other population parameters; `mle_medz_re` = median |z| over r_ levels.

| id | mle_nhyp | mle_maxz_hyp | mle_worst | mle_medz_re | frm_error | mle_error |
|---|---|---|---|---|---|---|
| p003 | 6 | 0.574 | sd_g__sigma_Intercept | 0.047 |  |  |
| p041 | 8 | 0.399 | sd_mmg1g2__y1_Intercept | 0.031 |  |  |
| p060 | 6 | 1.49 | cor_g__Intercept__kappa_Intercept | 0.124 |  |  |
| p062 | 5 | 0.193 | sigma | 0.028 |  |  |
| p075 | 5 | 0.102 | b_Intercept[3] |  |  |  |
| p082 | 4 | 68.986 | sdgp_gpx |  |  |  |
| p096 | 4 | 0.122 | shape |  |  |  |
| p121 | 4 | 0.48 | sd_g__Intercept | 0.029 |  |  |
| p152 | 4 | 0.249 | shape |  |  |  |
| big_crossed | 4 | 0.281 | sd_b__Intercept | 0.026 |  |  |
| big_pois | 3 | 0.254 | sd_g__Intercept | 0.058 |  |  |
| big_slope | 5 | 0.655 | sd_g__Intercept | 0.052 |  |  |

### Q2a. Survey instrument steps (s) and the no-Jacobian shim

Every program was in dev/stan-cache, so `t_rstan_compile` is rstan reusing the program it had just read, not a compile; the cold compile times are Q1a's. `d_*` = max relative difference, stanli shim against rstan, at frm()'s estimate and 5 points around it.

| id | key | npar_u | cache_hit | t_rstan_warm | t_rstan_chains0 | t_rstan_lp | t_stanli_setup | t_stanli_setup_jac | t_stanli_lp | d_noj | d_j | d_grad | d_unconstrain |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| p027 | importance x equidistant | 14 | TRUE | 1.11 | 0.09 | 0 | 0.09 | 0.07 | 0 | 2.057e-16 | 2.053e-16 | 1.776e-15 | 0 |
| p075 | cumulative x hypothesis_profile | 5 | TRUE | 0.66 | 0.07 | 0 | 0.11 | 0.07 | 0 | 0 | 0 | 2.401e-15 | 0 |
| p113 | us x gr_cov | 26 | TRUE | 0.71 | 0.1 | 0.02 | 0.09 | 0.1 | 0 | 0 | 0 | 1.11e-15 | 0 |
| p124 | binomial x importance | 13 | TRUE | 0.92 | 0.11 | 0 | 0.08 | 0.09 | 0 | 5.489e-16 | 5.48e-16 | 1.132e-14 | 0 |

### Q2b. L4 sampling (2 chains x 2000, flat priors) and stanli's rstan-shaped calls

| id | t_rstan_l4 | t_stanli_l4 | l4_nvar | l4_maxz | l4_crit | l4_rhat_rstan | l4_rhat_stanli | asm_chains0 | asm_noj | fp_error |
|---|---|---|---|---|---|---|---|---|---|---|
| p027 | 9.8 | 1.36 | 31 | 1.786 | 3.154 | 1.005 | 1.004 | stanli_stanfit | adjust_transform = FALSE is unsupported: stanli's compiled graph includes the parameter-transform Jacobian | pars not in the model: Lncor, Cortime, ordered_Intercept, fixed_Intercept, theta, Llncor, Lcortime, chol_cor, bQ, zb, zbsp, zbs, zar, zma, hs_local, R2D2_phi, scales, merged_Intercept, zcar, nszcar, z |
| p075 | 2.5 | 0.37 | 9 | 1.805 | 2.773 | 1.002 | 1.005 | stanli_stanfit | adjust_transform = FALSE is unsupported: stanli's compiled graph includes the parameter-transform Jacobian | pars not in the model: Lncor, Cortime, ordered_Intercept, fixed_Intercept, theta, Llncor, Lcortime, chol_cor, bQ, zb, zbsp, zbs, zar, zma, hs_local, R2D2_phi, scales, merged_Intercept, zcar, nszcar, z |
| p113 | 3.33 | 0.82 | 73 | 1.472 | 3.396 | 1.004 | 1.009 | stanli_stanfit | adjust_transform = FALSE is unsupported: stanli's compiled graph includes the parameter-transform Jacobian | pars not in the model: Lncor, Cortime, ordered_Intercept, fixed_Intercept, theta, Llncor, Lcortime, chol_cor, bQ, zb, zbsp, zbs, zar, zma, hs_local, R2D2_phi, scales, merged_Intercept, zcar, nszcar, z |
| p124 | 16.98 | 1.97 | 25 | 1.298 | 3.09 | 1.013 | 1.004 | stanli_stanfit | adjust_transform = FALSE is unsupported: stanli's compiled graph includes the parameter-transform Jacobian | pars not in the model: Lncor, Cortime, ordered_Intercept, fixed_Intercept, theta, Llncor, Lcortime, chol_cor, bQ, zb, zbsp, zbs, zar, zma, hs_local, R2D2_phi, scales, merged_Intercept, zcar, nszcar, z |

### Q2c. The method tier's Fixed_param brmsfit on as_stanmodel() (keep-list with include = TRUE)

| id | key | n_excl | n_excl_cut | t_fp_rstan | t_fp_stanli | epred_rel_diff | vars_same | fixef_rel_diff | error |
|---|---|---|---|---|---|---|---|---|---|
| p075 | cumulative x hypothesis_profile | 21 | 0 | 0.22 | 0.2 | 0 | TRUE | 0 |  |
| p124 | binomial x importance | 25 | 1 | 0.03 | 0.15 | 0 | TRUE | 0 |  |
| p027 | importance x equidistant | 25 | 1 | 0.05 | 0.16 | 0 | TRUE | 0 |  |
| p113 | us x gr_cov | 25 | 4 | 0.04 | 0.18 | 0 | TRUE | 0 |  |

### Q3a. Efficiency: min bulk ESS over every variable per second of wall time (whole call)

| id | n | n_par | n_random | tmb_total | tmb_sampling | tmb_miness | tmb_ess_s | sn_metric | sn_total | sn_sampling | sn_miness | sn_ess_s | ratio |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| big_crossed | 6000 | 344 | 340 | 790.5 | 785.256 | 795.4 | 1.006 | sparse | 388.77 | 376.003 | 2896.5 | 7.45 | 7.41 |
| big_pois | 10000 | 2003 | 2000 | 851.61 | 837.471 | 1511.7 | 1.775 | sparse | 224.97 | 202.04 | 2788 | 12.393 | 6.98 |
| gp | 60 | 64 | 60 | 368.17 | 366.818 | 105.7 | 0.287 | dense | 52.47 | 50.131 | 7.1 | 0.135 | 0.47 |
| slope | 400 | 86 | 80 | 6.74 | 4.711 | 508.9 | 75.504 | dense | 6.31 | 4.615 | 1118 | 177.179 | 2.35 |
| small_pois | 300 | 33 | 30 | 9.08 | 7.364 | 1094.5 | 120.54 | dense | 5.83 | 4.371 | 6274.6 | 1076.261 | 8.93 |
| ald | 200 | 25 | 20 | 7.81 | 5.49 | 674.791 | 86.401 |  | 1.25 |  |  |  |  |
| small_pois-laplace | 300 | 33 | 30 | 57.3 | 55.943 | 3261.6 | 56.921 | diag | 39.14 | 38.379 | 3151.7 | 80.524 | 1.41 |
| weak | 300 | 153 | 150 | 45.26 | 42.601 | 648.3 | 14.324 | sparse | 68.56 | 65.377 | 22.7 | 0.331 | 0.02 |

### Q3b. Diagnostics and agreement (tmbstan route vs SNUTS)

z as in Q1b over every variable both arms report; sd ratio = posterior sd tmbstan / SNUTS.

| id | tmb_div | tmb_maxtd | tmb_maxrhat | tmb_leap | sn_div | sn_maxtd | sn_maxrhat | sn_leap | nvar | maxz | maxz_hyp | crit | worst | sdratio_range | error |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| big_crossed | 0 | 4 | 1.007 | 60000 | 0 | 4 | 1.007 | 60000 | 346 | 2.346 | 1.696 | 3.8 | r_a[130,Intercept] | 0.927..1.083 |  |
| big_pois | 0 | 5 | 1.008 | 124016 | 0 | 4 | 1.007 | 60000 | 2004 | 3.616 | 1.152 | 4.215 | r_g[809,Intercept] | 0.92..1.093 |  |
| gp | 293 | 10 | 1.025 | 381984 | 1158 | 7 | 1.55 | 65804 | 64 | 2.067 | 2.067 | 3.359 | sdgp_gpx | 0.754..2.472 |  |
| slope | 0 | 5 | 1.007 | 61792 | 0 | 4 | 1.006 | 35912 | 89 | 1.963 | 1.142 | 3.449 | r_g[2,Intercept] | 0.952..1.184 |  |
| small_pois | 0 | 5 | 1.008 | 60584 | 0 | 3 | 1.004 | 28000 | 34 | 3.134 | 1.721 | 3.18 | r_g[5,Intercept] | 0.96..1.057 |  |
| ald | 0 | 6 | 1.01 | 126784 |  |  |  |  |  |  |  |  |  |  | error in evaluating the argument 'x' in selecting a method for function 't': error in evaluating the argument 'x' in selecting a method for function 't': leading principal minor of order 20 is not positive |
| small_pois-laplace | 0 | 3 | 1.002 | 19576 | 0 | 3 | 1.002 | 14928 | 4 | 0.791 | 0.791 | 2.498 | sd_g__Intercept | 0.992..1.028 |  |
| weak | 1 | 5 | 1.011 | 97762 | 53 | 7 | 1.134 | 159837 | 154 | 2.764 | 1.614 | 3.595 | r_g[128,Intercept] | 0.885..1.688 |  |

