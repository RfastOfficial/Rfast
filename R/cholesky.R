#[export]
cholesky <- function(x,parallel = FALSE) {
	if(parallel){
		.Call(Rfast_cholesky_par,x)
	}else{
  		.Call(Rfast_cholesky,x)
	}
}