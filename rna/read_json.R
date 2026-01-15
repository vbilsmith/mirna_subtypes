library("rjson")
json_file <- "rna/metadata.repository.2025-06-10.json"
json_data <- fromJSON(paste(readLines(json_file), collapse=""))
