#[export]
rep_col <- function(x,n) {
	.Call(Rfast_rep_col,x,n)
}

#[export]
rep_row <- function(x,n) {
	.Call(Rfast_rep_row,x,n)
}