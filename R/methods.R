
print.lod <- function(x, ...){
	cat("Data summary\n")
	print(x$dataSum)
	cat("\nStandards summary\n")
	print(x$standardsSum)
	cat("\nAssay summary\n")
	print(x$assaySum)
}

summary.lod <- function(x, ...){
	cat("Assay summary\n")
	print(x$assaySum)
}
