#[export]
Rank <- function(x,method = "average",descending = FALSE) {
  .Call(Rfast_rank,x,method,descending)
}

#[export]
rowRanks <- function(x,method = "average",descending = FALSE,stable = FALSE, parallel = FALSE) {
	if(parallel){
  		.Call(Rfast_row_ranks_p,x,method,descending,stable)
	}else{
		.Call(Rfast_row_ranks,x,method,descending,stable)
	}
}

#[export]
colRanks <- function(x,method = "average",descending = FALSE,stable = FALSE,parallel = FALSE) {
	if(parallel){
  		.Call(Rfast_col_ranks_p,x,method,descending,stable)
	}else{
		.Call(Rfast_col_ranks,x,method,descending,stable)
	}
}