#[export]
med <- function(x,na.rm=FALSE) {
#	.Deprecated("Rfast::Median")
  .Call(Rfast_med,x,na.rm)
}

#[export]
Median <- function(x,na.rm=FALSE) {
  .Call(Rfast_med,x,na.rm)
}

#[export]
colMedians <- function(x,na.rm=FALSE,parallel = FALSE) {
	.Call(Rfast_col_meds,x,na.rm,parallel)
}

#[export]
rowMedians <- function(x,na.rm=FALSE,parallel = FALSE) {
	.Call(Rfast_row_meds,x,na.rm,parallel)
}