#[export]
Table.sign <- function(x,names = TRUE,useNA = FALSE) {
  .Call(Rfast_table_sign,x,useNA,names)
}

#[export]
Table <- function(x,y=NULL,names = TRUE,useNA = FALSE,rm.zeros = FALSE) {
	if(names){
		if(is.null(y)){
			.Call(Rfast_table_with_names,x)
		}else{
			x <- .Call(Rfast_table2_with_names,x,y,rm.zeros)
			f <- x$f
			rownames(f) <- x$x
			colnames(f) <- x$y
			f
		}
	}else if(is.null(y)){
		.Call(Rfast_table_c,x,useNA)
	}else{
		.Call(Rfast_table2_c,x,y,rm.zeros)
	}
}