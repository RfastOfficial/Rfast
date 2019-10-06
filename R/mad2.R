#[export]
rowMads <- function(x,method = "median", na.rm = FALSE,parallel = FALSE) {
	.Call(Rfast_row_mads,x,method,na.rm,parallel)
}

#[export]
colMads <- function(x,method = "median", na.rm = FALSE,parallel = FALSE) {
	.Call(Rfast_col_mads,x,method,na.rm,parallel)
}

#[export]
mad2 <- function(x,method = "median",na.rm = FALSE) {
	#.Deprecated("Rfast::Mad")
	.Call(Rfast_mad2,x,method,na.rm)
}


#[export]
Mad <- function(x,method = "median",na.rm = FALSE) {
	.Call(Rfast_mad2,x,method,na.rm)
}