#[export]
data.frame.to_matrix <- function(x,col.names = NULL,row.names = NULL) {
	if(is.null(col.names) && is.null(row.names)){
		x <- .Call("Rfast_frame_to_matrix",x)
	}else{
		y <- .Call("Rfast_frame_to_matrix",x)
		if(!is.null(col.names) && col.names == TRUE){
			colnames(y) <- colnames(x)
		}
		if(!is.null(row.names) && row.names == TRUE){
			rownames(y) <- rownames(x)
		}
	}
	x
}