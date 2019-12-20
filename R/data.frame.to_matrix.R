#[export]
data.frame.to_matrix <- function(x,col.names = NULL,row.names = NULL) {
	if(is.null(col.names) && is.null(row.names)){
		x <- .Call(Rfast_frame_to_matrix,x)
	}else{
		if(length(col.names) == 1 && col.names == TRUE){
			col.names <- colnames(x)
		}
		if(length(row.names) == 1 && row.names == TRUE){
			row.names <- rownames(x)
		}
		x <- .Call(Rfast_frame_to_matrix,x)
		colnames(x) <- col.names
		rownames(x) <- row.names
	}
	x
}
