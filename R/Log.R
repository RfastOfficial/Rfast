#[export]
Log <- function(x, na.rm = FALSE) {
	if(na.rm){
		x <- x[-is.na(x)]
	}
  	.Call(Rfast_Log,x)
}