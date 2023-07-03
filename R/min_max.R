#[export]
rowMinsMaxs <- function(x) {
	x <- .Call(Rfast_row_min_max,x)
	rownames(x) <- c("min","max")
	x
}

#[export]
rowMins <- function(x,value=FALSE) {
	if(value){
		.Call(Rfast_row_min,x)
	}else{
		.Call(Rfast_row_min_indices,x)
	}
}

#[export]
rowMaxs <- function(x,value=FALSE) {
	if(value){
		.Call(Rfast_row_max,x)
	}else{
		.Call(Rfast_row_max_indices,x)
	}
}

#[export]
colMins <- function(x,value=FALSE,parallel = FALSE) {
	if(value){
		.Call(Rfast_col_min,x,parallel)
	}else{
    	.Call(Rfast_col_min_indices,x)
	}
}

#[export]
colMaxs <- function(x,value=FALSE,parallel = FALSE) {
	if(value){
		.Call(Rfast_col_max,x,parallel)
	}else{
    	.Call(Rfast_col_max_indices,x)
	}
}

#[export]
colMinsMaxs <- function(x,parallel = FALSE) {
	x <- .Call(Rfast_col_min_max,x,parallel)
	rownames(x) <- c("min","max")
	x
}

#[export]
min_max<-function(x,index=FALSE,percent=FALSE){
	if(percent){
		x <- .Call(Rfast_min_max_perc,x)
		names(x) <- c("min","max","negative%","positive%")
		return(x)
	}
	x <- .Call(Rfast_min_max,x,index)
	names(x) <- c("min","max")
	x
}