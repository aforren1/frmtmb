| package | section (as filed) | source | record | bullet |
|---|---|---|---|---|
| frmtmb | Breaking changes | smallfix | `dev/smallfix-findings.md` | A `dirichlet()` prior on a parameter that is not a simplex i |
| frmtmb | Bug fixes | smallfix | `dev/smallfix-findings.md` | A saved fit that needs extension packages this session canno |
| frmtmb | Bug fixes | smallfix | `dev/smallfix-findings.md` | `fixef()` of a fit with a part from a package this session c |
| frmtmb | Bug fixes | smallfix | `dev/smallfix-findings.md` | `predict(newdata = )` names the reason for the model at hand |
| frmtmb | Bug fixes | smallfix | `dev/smallfix-findings.md` | `print()` of a curve (`frm_curve()`, `frm_curve_deriv()`) or |
| frmtmb | Other changes | smallfix | `dev/smallfix-findings.md` | frmtmb's R CMD check (four platforms) and its coverage job f |
| frmtmb | Other changes | smallfix | `dev/smallfix-findings.md` | The pkgdown job, the one that installs every extension, also |
| frmtmb | Breaking changes | bhaz | `dev/bhaz-findings.md` | brms's `bhaz()` addition term is frmtmb.survival's. |
| frmtmb | New features | bhaz | `dev/bhaz-findings.md` | `frmtmb_register_aterm(option =)` registers an addition term |
| frmtmb | Breaking changes | brmsnames | `dev/brmsnames-1024-findings.md` | A term that calls `exp()` has brms's coefficient names. To k |
| frmtmb | Breaking changes | brmsnames | `dev/brmsnames-1024-findings.md` | `default_prior()` lists each coefficient by brms's name, and |
| frmtmb | Breaking changes | brmsnames | `dev/brmsnames-1024-findings.md` | The rows of `default_prior()`, `validate_prior()` and |
| frmtmb | Breaking changes | brmsnames | `dev/brmsnames-1024-findings.md` | A nonlinear parameter of a multivariate model has brms's pre |
| frmtmb | Breaking changes | brmsnames | `dev/brmsnames-1024-findings.md` | `fixef()`, `vcov()` and `summary()` of a fit list the rows i |
| frmtmb | Bug fixes | brmsnames | `dev/brmsnames-1024-findings.md` | `drop1()` and `step()` see a term that calls `exp()`. With t |
| frmtmb | Bug fixes | brmsnames | `dev/brmsnames-1024-findings.md` | `re_formula` can name a group-level `exp()` term. |
| frmtmb | Bug fixes | brmsnames | `dev/brmsnames-1024-findings.md` | A smooth's unpenalized column written in both spellings is o |
| frmtmb | Bug fixes | brmsnames | `dev/brmsnames-1024-findings.md` | A multivariate model whose second response's nonlinear formu |
| frmtmb | New features | sarfcor | `dev/sarfcor-findings.md` | `fcor(M)`, brms's fixed residual correlation term, is suppor |
| frmtmb | New features | sarfcor | `dev/sarfcor-findings.md` | brms 2.23.0 checks that `M` is positive definite with |
| frmtmb | New features | sarfcor | `dev/sarfcor-findings.md` | brms's `sar()` is now supported in frmtmb.spatial. Without t |
| frmtmb | For extension authors | sarfcor | `dev/sarfcor-findings.md` | New seam `frmtmb_register_autocor()`. A package registers a  |
| frmtmb | For extension authors | sarfcor | `dev/sarfcor-findings.md` | The map gets one `sigma` when the term sets `sigma_varies =  |
| frmtmb | For extension authors | sarfcor | `dev/sarfcor-findings.md` | A term's `bounds` slot states the interval of a parameter th |
| frmtmb | For extension authors | sarfcor | `dev/sarfcor-findings.md` | A slot of `frmtmb_register_autocor()` or |
| frmtmb | For extension authors | sarfcor | `dev/sarfcor-findings.md` | New export `frmtmb_covstruct_info()`: what the registry know |
| frmtmb | For extension authors | sarfcor | `dev/sarfcor-findings.md` | `frmtmb_register_covstruct(full_rank =)` states whether an e |
| frmtmb | For extension authors | sarfcor | `dev/sarfcor-findings.md` | autoscale leaves the random slopes of every structure but `u |
| frmtmb | Bug fixes | sarfcor | `dev/sarfcor-findings.md` | The weights check, which refuses a weighted likelihood with  |
| frmtmb | Bug fixes | sarfcor | `dev/sarfcor-findings.md` | The Student-t latents of `gr(dist = "student")` were not ful |
| frmtmb.latent | Other changes | smallfix | `dev/smallfix-findings.md` | The `mixture_mvn()` tests compare each result as a ratio to  |
| frmtmb.latent | Other changes | smallfix | `dev/smallfix-findings.md` | `predict(newdata = )` of a `mixture_mvn()` fit says why it i |
| frmtmb.sample | Breaking changes | brmsnames | `dev/brmsnames-1024-findings.md` | The draws of a model with an `exp()` term carry brms's names |
| frmtmb.sample | Breaking changes | brmsnames | `dev/brmsnames-1024-findings.md` | `fixef()` of the draws lists the rows in brms's order, which |
| frmtmb.sample | New features | sarfcor | `dev/sarfcor-findings.md` | The draws of a model with frmtmb.spatial's `sar()` carry brm |
| frmtmb.sample | New features | sarfcor | `dev/sarfcor-findings.md` | The default prior on `lagsar` and `errorsar` is brms's: flat |
| frmtmb.sample | Bug fixes | sarfcor | `dev/sarfcor-findings.md` | The derived `sd_` and `cor_` columns asked core's |
| frmtmb.survival | New features | bhaz | `dev/bhaz-findings.md` | brms's `bhaz()` addition term. `time \| cens(c) + bhaz(df = 6 |
| frmtmb.survival | New features | bhaz | `dev/bhaz-findings.md` | `cox(knots =, Boundary.knots =)` set the internal and the bo |
| frmtmb.survival | New features | bhaz | `dev/bhaz-findings.md` | An argument that `bhaz()` and `cox()` both set must have one |
| frmtmb.survival | New features | bhaz | `dev/bhaz-findings.md` | `bhaz(gr = g)`, brms's stratified Cox model with one baselin |
| frmtmb.survival | New features | bhaz | `dev/bhaz-findings.md` | A shortened argument name that brms reads is read the same w |
| frmtmb.survival | New features | bhaz | `dev/bhaz-findings.md` | Boundary knots that do not cover every observed time stop. b |
| frmtmb.survival | Other changes | bhaz | `dev/bhaz-findings.md` | `cox(bhaz = list(df = 6))`, the spelling of brms before 2.22 |
| frmtmb.survival | Other changes | bhaz | `dev/bhaz-findings.md` | A `df` below `degree + intercept` is still refused when the  |
| frmtmb.spatial | New features | sarfcor | `dev/sarfcor-findings.md` | `sar(M, type = )`, brms's spatial simultaneous autoregressiv |
| frmtmb.spatial | New features | sarfcor | `dev/sarfcor-findings.md` | `fitted()`, `predict()` and `simulate()` take `newdata` with |
| frmtmb.spatial | New features | sarfcor | `dev/sarfcor-findings.md` | `rho` is estimated on `(1 / min eigenvalue, 1 / max eigenval |
| frmtmb.spatial | New features | sarfcor | `dev/sarfcor-findings.md` | A matrix whose smallest or largest eigenvalue repeats for al |
| frmtmb.spatial | New features | sarfcor | `dev/sarfcor-findings.md` | A row-standardized `M` of an irregular graph, spdep's defaul |
| frmtmb.spatial | New features | sarfcor | `dev/sarfcor-findings.md` | spdep's `nb` and `listw` objects are read as brms reads them |
