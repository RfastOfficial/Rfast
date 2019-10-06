#[export]
binary_search <- function(x, v, index=FALSE) {
	if(index){
  		return (.Call(Rfast_lowerbound,x,v))
	}
	.Call(Rfast_binarysearch,x,v)
}