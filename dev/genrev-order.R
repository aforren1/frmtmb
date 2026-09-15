b <- deparse(base::loadNamespace)
i <- grep("onLoad|registerS3methods|namespaceExport|sealNamespace|lockEnvironment|S3methods|onAttach", b)
for (k in i) cat(sprintf("%5d: %s\n", k, b[k]))
