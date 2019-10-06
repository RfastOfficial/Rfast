#[export]
columns <- function(x,indices) {
	.Call(Rfast_columns,x,indices)
}

#[export]
rows <- function(x,indices) {
	.Call(Rfast_rows,x,indices)
}

#[export]
XopY.sum<-function(x,y=NULL,oper="*"){
	if(is.null(y)){
		.Call(Rfast_sum_XopX,x,oper)
	}
	else{
		.Call(Rfast_sum_XopY,x,y,oper)
	}
}

#[export]
submatrix <- function(x,rowStart=1,rowEnd=1,colStart=1,colEnd=1) {
  .Call(Rfast_submatrix,x,rowStart,rowEnd,colStart,colEnd)
}

#[export]
transpose <- function(x) {
	.Call(Rfast_transpose,x)
}

#[export]
mat.mult <- function(x,y) {
	.Call(Rfast_mat_mult_p,x,y)
}

#[export]
mat.mat <- function(x, y) {
	.Call(Rfast_mat_mat,x,y)
}