#[export]
dirknn <- function(xnew, x, y, k, type = "C", parallel = FALSE) {
  .Call(Rfast_dir_knn,t(xnew),t(x),y,k,type,parallel)
}
