#[export]
colsums <- function(x,indices = NULL,parallel = FALSE,na.rm = FALSE,cores = 0) {
	if(parallel){
		.Call(Rfast_col_sums_p,x,cores)
	}else{
  		.Call(Rfast_col_sums,x,indices,na.rm)
	}
}

#[export]
rowsums <- function(x,indices = NULL,parallel = FALSE,na.rm = FALSE,cores = 0) {
	if(parallel){
  		.Call(Rfast_row_sums_p,x,cores)
  	}else{
  		.Call(Rfast_row_sums,x,indices,na.rm)
  	}
}