#[export]
dista <- function(xnew,x,type = "euclidean",k=0,index=FALSE,trans = TRUE,square = FALSE,p = 0,parallel = FALSE) {
	if(k>0){
		if(index){
			x <- .Call(Rfast_dista_index,t(xnew),t(x),k,type)
		}else{
			x <- .Call(Rfast_dista_values,t(xnew),t(x),k,square,type)
		}
	}else{
		x <- .Call(Rfast_dista,t(xnew),t(x),square,type,p,parallel)
	}
	if(trans)	x <- t(x)
	x
}