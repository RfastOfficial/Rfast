#[export]
rowVars <- function(x, suma = NULL, std = FALSE, na.rm = FALSE,parallel = FALSE) {
	.Call(Rfast_row_vars,x,std,na.rm,parallel)
}

#[export]
colVars <- function(x, suma = NULL, std = FALSE, na.rm = FALSE, parallel = FALSE) {
	.Call(Rfast_col_vars,x,std,na.rm,parallel)
}

#[export]
Var <- function(x,std = FALSE,na.rm = FALSE) {
	.Call(Rfast_var,x,std,na.rm)
}