

#[export]
Dist <- function(x,method = "euclidean", square = FALSE,p=0,vector = FALSE) {
	if(method != "harvesine")
  		x <- t(x)
	if(method == "hellinger"){
		x <- sqrt(x)
	}
	if(vector){
		.Call(Rfast_dist_vec,x,method,square,p)
	}else{
		.Call(Rfast_dist,x,method,square,p)
	}
}