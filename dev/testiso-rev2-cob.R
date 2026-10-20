# Reviewer scratch: what filename= and <source> covr's to_cobertura writes
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
cat(packageVersion("covr") |> as.character(), "\n")
x <- deparse(covr::to_cobertura)
i <- grep("filename|source|attr|display_name|relative|root", x)
cat(x[i], sep = "\n")
