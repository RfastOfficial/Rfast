#[export]
coldiffs <- function(x) {
  .Call(Rfast_col_diffs,x)
}