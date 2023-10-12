
#[export]
rowAll <- function(x,parallel = FALSE,cores = 0) {
	if(parallel){
		.Call(Rfast_row_all_p,x,cores)
	}else{
		.Call(Rfast_row_all,x)
	}
}

#[export]
colAll <- function(x,parallel = FALSE,cores = 0) {
	if(parallel){
		.Call(Rfast_col_all_p,x,cores)
	}else{
		.Call(Rfast_col_all,x)
	}
}