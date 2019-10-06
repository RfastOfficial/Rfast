
#[export]
dvar <- function(x) {
  .Call(Rfast_dvar,t(x))
}