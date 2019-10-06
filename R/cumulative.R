#[export]
colCumProds <- function(x) {
	.Call(Rfast_col_cum_prods,x)
}

#[export]
colCumMins <- function(x) {
	.Call(Rfast_col_cum_mins,x)
}

#[export]
colCumSums <- function(x) {
	.Call(Rfast_col_cum_sums,x)
}

#[export]
colCumMaxs <- function(x) {
	.Call(Rfast_col_cum_maxs,x)
}