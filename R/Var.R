#[export]
rowVars <- function(x, std = FALSE, na.rm = FALSE,parallel = FALSE) {
	.Call(Rfast_row_vars,x,std,na.rm,parallel)
}

#[export]
colVars <- function(x, std = FALSE, na.rm = FALSE, parallel = FALSE) {
	UseMethod("colVars")
}

#[export s3]
colVars.matrix <- function(x, std = FALSE, na.rm = FALSE, parallel = FALSE) {
	.Call(Rfast_col_vars,x,std,na.rm,parallel)
}

#[export s3]
colVars.data.frame <- function(x, std = FALSE, na.rm = FALSE, parallel = FALSE) {
	.Call(Rfast_col_vars,x,std,na.rm,parallel)
}

#[export]
Var <- function(x,std = FALSE,na.rm = FALSE) {
	.Call(Rfast_var,x,std,na.rm)
}