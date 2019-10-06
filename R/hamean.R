#[export]
colhameans <- function(x, parallel = FALSE)  {
	dim(x)[1]/ Rfast::colsums(1/x, parallel = parallel)
}

#[export]
rowhameans <- function(x)  {
	dim(x)[1]/ Rfast::rowsums(1/x)
}