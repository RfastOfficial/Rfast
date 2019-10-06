#[export]
Lchoose <- function(x,k) {
  .Call(Rfast_Lchoose,x,k)
}

#[export]
Choose <- function(x,k) {
  .Call(Rfast_Choose,x,k)
}

#[export]
Trigamma <- function(x) {
  .Call(Rfast_Trigamma,x)
}

#[export]
Digamma <- function(x) {
  .Call(Rfast_Digamma,x)
}

#[export]
Lgamma <- function(x) {
  .Call(Rfast_Lgamma,x)
}