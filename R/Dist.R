

#[export]
Dist <- function(x,method = "euclidean", square = FALSE,p=0,vector = FALSE) {
	if (method == "canberra1" || method == "canberra2") {
        .Deprecated("The replacement type is \"canberra\"", "Rfast")
    }
	if(method != "haversine")
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

#[export]
edist <- function(x, y = NULL){
	if(is.list(x)){
		## x is a list
	  p <- length(x)
	  dii <- numeric(p)
	  dis <- matrix(0, p, p)
	  ni <- numeric(p)
	  for (i in 1:p) { 
	    dii[i] <- Rfast::total.dist(x[[ i ]])
	    ni[i] <- dim(x[[ i ]])[1]  ## poses grammes exei
	  }
	  for ( i in 1:(p - 1) ) {
	    for ( j in (i + 1):p ) {
	      mij <- Rfast::total.dista(x[[ i ]], x[[ j ]])
	      dis[i, j] <- (2 * mij - ni[ j ] * dii[i] / ni[ i ] - ni[i] * 
	                   dii[j]/ni[j] ) / (ni[i] + ni[j])
	      dis[j, i] <- dis[i, j]
	    }
	  }
	  dis
	}else{
    	.Call(Rfast_edist, t(x), t(y))
    }
}

#[export]
total.dist <- function(x,method = "euclidean", square = FALSE,p=0) {
  if (method == "canberra1" || method == "canberra2") {
        .Deprecated("The replacement method is \"canberra\"", "Rfast")
    }
  if(method != "haversine")
	x <- t(x)
  if(method == "hellinger"){
	x <- sqrt(x)
  }
  .Call(Rfast_total_dists,x,method,square,p)
}

#[export]
vecdist <- function(x) {
  	.Call(Rfast_vecdist,x)
}