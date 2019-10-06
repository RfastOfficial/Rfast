#[export]
colAny <- function(x) {
  .Call(Rfast_col_any,x)
}

#[export]
rowAny <- function(x) {
  .Call(Rfast_row_any,x)
}