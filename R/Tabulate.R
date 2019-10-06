
#[export]
rowTabulate <- function(x, max_number = max(x)) {
  .Call(Rfast_row_tabulate,x,max_number)
}

#[export]
colTabulate <- function(x, max_number = max(x)) {
  .Call(Rfast_col_tabulate,x,max_number)
}