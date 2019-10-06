
#[export]
lower_tri <- function(x,suma = FALSE,diag = FALSE) {
	if(is.vector(x) && length(x)==2){
		.Call(Rfast_lower_tri_b, x[1L],x[2L],diag)
	}else if(suma){
		.Call(Rfast_sum_lower_tri, x,diag)
	}else{
		.Call(Rfast_lower_tri, x,diag)
	}
}

#[export]
upper_tri <- function(x,suma = FALSE,diag = FALSE) {
  	if(is.vector(x) && length(x)==2){
		.Call(Rfast_upper_tri_b,x[1L],x[2L],diag)
	}else if(suma){
		.Call(Rfast_sum_upper_tri,x,diag)	
	}else{
		.Call(Rfast_upper_tri,x,diag)
	}
}

#[export]
upper_tri.assign <- function(x,v,diag = FALSE) {
	.Call(Rfast_upper_tri_assign,x,v,diag)
}

#[export]
lower_tri.assign <- function(x,v,diag = FALSE) {
	.Call(Rfast_lower_tri_assign,x,v,diag)
}