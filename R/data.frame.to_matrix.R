#[export]
data.frame.to_matrix <- function(x,col.names = NULL,row.names = NULL) {
	if(is.null(col.names) && is.null(row.names)){
		x <- .Call(Rfast_frame_to_matrix,x)
	}else{
		if(col.names == TRUE){
			col.names <- colnames(x)
		}
		if(row.names == TRUE){
			row.names <- rownames(x)
		}
		x <- .Call(Rfast_frame_to_matrix,x)
		colnames(x) <- col.names
		rownames(x) <- row.names
	}
	x
}