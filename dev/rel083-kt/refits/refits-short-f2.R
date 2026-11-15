# Short check of the refits fix round, F2: resolve_cluster(). Seed 1,
# 120 rows (probe2 data).
source("C:/Users/adf44/source/r/frmtmb-wt-rel083/dev/rel083-kt/refits/refits-review-lib.R")
cat("build has fix:",
    any(grepl("fit_refit_data", deparse(frmtmb:::resolve_cluster))), "\n")
G <- 20
mkd <- function() {
  set.seed(1)
  dd <- data.frame(g = factor(rep(seq_len(G), each = 6)),
                   cl = factor(rep(seq_len(G / 2), each = 12)),
                   x = rnorm(G * 6))
  dd$y <- 1 + 0.5 * dd$x + rnorm(G * 6) + rnorm(G / 2)[dd$cl]
  dd
}
se_of <- function(v) round(sqrt(diag(v)), 7)
dd <- mkd()
fit <- frm(y ~ x, data = dd)
ref <- se_of(vcov_cluster(fit, dd$cl))
go("S0 explicit dd$cl", ref)

# a cluster variable absent from the data, present in the formula env
cl2 <- factor(rep(1:6, each = 20))
d_no <- dd[, c("y", "x", "g")]
f_no <- frm(y ~ x, data = d_no)
go("S1 ~ cl2, cl2 only in the workspace = explicit cl2",
   list(se_of(vcov_cluster(f_no, ~ cl2)), se_of(vcov_cluster(f_no, cl2))))
h <- function(fit) { cl3 <- factor(rep(1:4, each = 30)); vcov_cluster(fit, ~ cl3) }
go("S2 ~ cl3 local to a function = explicit",
   list(se_of(h(f_no)), se_of(vcov_cluster(f_no, factor(rep(1:4, each = 30))))))
# the same with NA rows dropped: the workspace vector has length n_data
d_na <- d_no; d_na$x[c(5, 60)] <- NA
f_na <- frm(y ~ x, data = d_na)
go("S3 ~ cl2 (length 120) with 2 NA rows dropped", se_of(vcov_cluster(f_na, ~ cl2)))
go("S3 explicit cl2[-c(5, 60)]", se_of(vcov_cluster(f_na, cl2[-c(5, 60)])))
cl118 <- cl2[-c(5, 60)]
go("S4 ~ cl118 (length 118) = explicit", se_of(vcov_cluster(f_na, ~ cl118)))

# a cluster variable that is also a model term
fg <- frm(y ~ x + cl, data = dd)
go("S5 ~ cl, cl a fixed term = explicit dd$cl",
   list(se_of(vcov_cluster(fg, ~ cl)), se_of(vcov_cluster(fg, dd$cl))))
fr <- frm(y ~ x + (1 | cl), data = dd)
go("S6 ~ cl, cl a grouping factor = explicit dd$cl",
   list(se_of(vcov_cluster(fr, ~ cl)), se_of(vcov_cluster(fr, dd$cl))))
cl <- factor(rep(1:2, length.out = 120))
go("S7 S5 with a global cl of length n", se_of(vcov_cluster(fg, ~ cl)))
go("S8 S6 with a global cl of length n", se_of(vcov_cluster(fr, ~ cl)))
rm(cl)
# a model term with NA rows: the fit's own column, aligned
dna <- dd; dna$x[c(3, 50)] <- NA
fgn <- frm(y ~ x + cl, data = dna)
go("S9 ~ cl model term, NA rows = explicit",
   list(se_of(vcov_cluster(fgn, ~ cl)), se_of(vcov_cluster(fgn, dna$cl[-c(3, 50)]))))
# an expression over a data column and a workspace column
k <- 2L
go("S10 ~ interaction(cl, g %% k): data cl + data g + global k",
   { v <- vcov_cluster(fit, ~ interaction(cl, as.integer(g) %% k)); se_of(v) })
go("S10 explicit",
   se_of(vcov_cluster(fit, interaction(dd$cl, as.integer(dd$g) %% k))))
# a cluster variable in the data and the workspace, data have NA rows
dna2 <- mkd(); dna2$x[c(3, 50)] <- NA
fna2 <- frm(y ~ x, data = dna2)
cl <- factor(rep(1:2, length.out = 120))
go("S11 global cl beside data cl, NA rows",
   list(se_of(vcov_cluster(fna2, ~ cl)), se_of(vcov_cluster(fna2, dna2$cl[-c(3, 50)]))))
rm(cl)
# character name
go("S12 cluster = 'g' (character, model frame lacks g)",
   se_of(vcov_cluster(fit, "g")))
go("S12 cluster = ~ g", se_of(vcov_cluster(fit, ~ g)))
go("S12 explicit dd$g", se_of(vcov_cluster(fit, dd$g)))
