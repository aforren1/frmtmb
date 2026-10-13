# Reviewer of lane optrap: rows of dev/optrap-rev-lmprec.R identical
# between the builds (both pinned to P cores), by spread.
a <- read.delim("dev/optrap-rev-log/lmprec-lane.tsv", colClasses = "character")
b <- read.delim("dev/optrap-rev-log/lmprec-r7.tsv", colClasses = "character")
same <- a$beta_se == b$beta_se & a$sigma_se == b$sigma_se & a$llgap == b$llgap
print(tapply(same, a$sd, sum))
print(tapply(a$unit, a$sd, function(u) paste(unique(u), collapse = ",")))
