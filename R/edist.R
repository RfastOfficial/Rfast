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
    	.Call("Rfast_edist", PACKAGE = "Rfast", t(x), t(y))
    }
}
 
