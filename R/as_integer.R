#[export]
as_integer <- function(x,result.sort = TRUE,init = 1, parallel = FALSE) {
  	.Call(Rfast_as_integer,x,result.sort,init,parallel)
}