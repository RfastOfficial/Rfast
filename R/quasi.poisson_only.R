#[export]
quasi.poisson_only <- function(x,y,tol = 1e-09, maxiters = 100) {
	.Call(Rfast_quasi_poisson_only,x,y,sum(y*log(y),na.rm=TRUE),tol,maxiters)
}