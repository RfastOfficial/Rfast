#[export]
colmeans <- function(x,parallel = FALSE) {
	if(parallel){
		.Call(Rfast_col_mean_p,x)
	}else{
  		as.vector(.Call(Rfast_col_means,x))
  	}
}

#[export]
rowmeans <- function(x) {
  	as.vector(.Call(Rfast_row_means,x))
}