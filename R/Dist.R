

#[export]
Dist <- function(x,method = "euclidean", square = FALSE,p=0, result = "matrix" ,vector = FALSE, parallel = FALSE) {
	if (vector) {
        .Deprecated("Use options \"result\" instead", "Rfast")
    }
	if (method == "canberra1" || method == "canberra2") {
        .Deprecated("The replacement type is \"canberra\"", "Rfast")
    }
	if(method != "haversine")
  		x <- t(x)
	if(result == "vector"){
		.Call(Rfast_dist_vec,x,method,square,p, parallel)
	}else if(result == "matrix"){
		.Call(Rfast_dist,x,method,square,p, parallel)
	}else if(result == "sum"){
  		.Call(Rfast_total_dists,x,method,square,p, parallel)
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
	   for ( i in 1:p ) {
            dii[i] <- 2 * Rfast::Dist( x[[ i ]], result = "sum" )
            ni[i] <- dim( x[[ i ]] )[1]
        }
        for ( i in 1:(p - 1) ) {
            for ( j in (i + 1):p ) {
                mij <- Rfast::dista( x[[ i ]], x[[ j ]], result = "sum" )
                dis[i, j] <- ( 2 * mij - ni[j] * dii[i] / ni[i] - 
                  ni[i] * dii[j] / ni[j] ) / (ni[i] + ni[j])
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
  .Deprecated("Use Dist(x, result = \"sum\")", "Rfast")
  if (method == "canberra1" || method == "canberra2") {
        .Deprecated("The replacement method is \"canberra\"", "Rfast")
    }
  if(method != "haversine")
	x <- t(x)
  .Call(Rfast_total_dists,x,method,square,p,FALSE)
}

#[export]
vecdist <- function(x) {
  	.Call(Rfast_vecdist,x)
}

#[export]
coeff <- function(x,method,vector = FALSE) {
	if(vector){
		.Call(Rfast_coeff_vec,t(x),method)
	}else{
		.Call(Rfast_coeff,t(x),method)
	}
}