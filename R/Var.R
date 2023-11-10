#[export]
rowVars <- function(x, std = FALSE, na.rm = FALSE,parallel = FALSE,cores = 0) {
	.Call(Rfast_row_vars,x,std,na.rm,parallel,cores)
}

#[export]
colVars <- function(x, std = FALSE, na.rm = FALSE, parallel = FALSE,cores = 0) {
	UseMethod("colVars")
}

#[export s3]
colVars.matrix <- function(x, std = FALSE, na.rm = FALSE, parallel = FALSE,cores = 0) {
	.Call(Rfast_col_vars,x,std,na.rm,parallel,cores)
}

#[export s3]
colVars.data.frame <- function(x, std = FALSE, na.rm = FALSE, parallel = FALSE,cores = 0) {
	.Call(Rfast_col_vars,x,std,na.rm,parallel,cores)
}

#[export]
Var <- function(x,std = FALSE,na.rm = FALSE) {
	.Call(Rfast_var,x,std,na.rm)
}