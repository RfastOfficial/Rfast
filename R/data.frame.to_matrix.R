#[export]
data.frame.to_matrix <- function(x,col.names = NULL,row.names = NULL) {
	y <- .Call("Rfast_frame_to_matrix",x)
	if(is.character(col.names) && length(col.names)==ncol(y)){
		colnames(y) <- col.names
	}else if(col.names == TRUE){
		colnames(y) <- colnames(x)
	}
	if(is.character(row.names) && length(row.names)==nrow(y)){
		rownames(y) <- row.names
	}else if(col.names == TRUE){
		rownames(y) <- rownames(x)
	}
	y
}
