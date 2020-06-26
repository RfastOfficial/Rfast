#[export]
rownth <- function(x,elems,num.of.nths=1,descending=FALSE,na.rm = FALSE,index.return = FALSE,parallel = FALSE) {
	if(parallel){
		.Call(Rfast_row_nth_p,x,elems,num.of.nths,descending,na.rm,index.return)
	}else{
  		.Call(Rfast_row_nth,x,elems,num.of.nths,descending,na.rm,index.return)
	}
}

#[export]
colnth <- function(x,elems,num.of.nths=1,descending=FALSE,na.rm = FALSE,index.return = FALSE,parallel = FALSE) {
	if(parallel){
  		.Call(Rfast_col_nth_p,x,elems,num.of.nths,descending,na.rm,index.return)
	}else{
		.Call(Rfast_col_nth,x,elems,num.of.nths,descending,na.rm,index.return)
	}
}

#[export]
nth <- function(x,k,num.of.nths=1,descending=FALSE,index.return=FALSE,na.rm = FALSE) {
	if(is.integer(x)){
		.Call(Rfast_nth_int,x,k)
	}else{
  		.Call(Rfast_nth,x,k,num.of.nths,descending,na.rm,index.return)
  	}
}