#[export]
colPmax <- function(x,y) {
  .Call(Rfast_col_pmax,x,y)
}

#[export]
colPmin <- function(x,y) {
  .Call(Rfast_col_pmin,x,y)
}

#[export]
Pmin_Pmax <- function(x,y,na.rm = FALSE) {
  .Call(Rfast_pmin_pmax,x,y,na.rm)
}

#[export]
Pmax <- function(x,y,na.rm = FALSE) {
  .Call(Rfast_pmax,x,y,na.rm)
}

#[export]
Pmin <- function(x,y,na.rm = FALSE) {
  .Call(Rfast_pmin,x,y,na.rm)
}