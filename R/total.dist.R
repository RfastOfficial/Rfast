

#[export]
total.dist <- function(x,method = "euclidean", square = FALSE,p=0) {
  if(method != "harvesine")
	x <- t(x)
  if(method == "hellinger"){
	x <- sqrt(x)
  }
  .Call(Rfast_total_dists,x,method,square,p)
}