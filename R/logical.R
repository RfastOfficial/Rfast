#[export]
rowTrueFalse <- function(x) {
  x <- .Call(Rfast_row_true_false,x)
  rownames(x) <- c("FALSE","TRUE")
  x
}

#[export]
colTrueFalse <- function(x) {
  x <- .Call(Rfast_col_true_false,x)
  rownames(x) <- c("FALSE","TRUE")
  x
}

#[export]
colTrue <- function(x) {
  .Call(Rfast_col_true,x)
}

#[export]
rowTrue <- function(x) {
  .Call(Rfast_row_true,x)
}

#[export]
rowFalse <- function(x) {
  .Call(Rfast_row_false,x)
}

#[export]
colFalse <- function(x) {
  .Call(Rfast_col_false,x)
}