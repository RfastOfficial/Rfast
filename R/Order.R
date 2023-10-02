#[export]
Order<-function(x,stable=FALSE,descending=FALSE,partial=NULL,parallel=FALSE){
	if(is.character(x)){
		x <- as.numeric(x)
	}
	if(is.null(partial)){
		.Call(Rfast_Order,x,stable,descending,parallel)
	}else{
		.Call(Rfast_partial_sort_index,x,partial,descending,parallel)
	}
}

#[export]
rowOrder <- function(x,stable=FALSE,descending=FALSE,parallel=FALSE,cores = 0) {
  if(parallel){
  	.Call(Rfast_row_order_p,x,stable,descending,cores)
  }else{
  	.Call(Rfast_row_order,x,stable,descending)
  }
}

#[export]
colOrder <- function(x,stable=FALSE,descending=FALSE,parallel=FALSE,cores = 0) {
  if(parallel){
  	.Call(Rfast_col_order_p,x,stable,descending,cores)
  }else{
  	.Call(Rfast_col_order,x,stable,descending)
  }
}