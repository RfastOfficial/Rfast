#[export]
Rank <- function(x,method = "average",descending = FALSE) {
  .Call(Rfast_rank,x,method,descending)
}

#[export]
rowRanks <- function(x,method = "average",descending = FALSE,stable = FALSE, parallel = FALSE,cores = 0) {
	if(parallel){
  		.Call(Rfast_row_ranks_p,x,method,descending,stable,cores)
	}else{
		.Call(Rfast_row_ranks,x,method,descending,stable)
	}
}

#[export]
colRanks <- function(x,method = "average",descending = FALSE,stable = FALSE,parallel = FALSE,cores = 0) {
	.Call(Rfast_col_ranks,x,method,descending,stable,parallel,cores)
}