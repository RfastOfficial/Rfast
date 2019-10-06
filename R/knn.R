
#[export]
knn <- function(xnew, y, x, k, dist.type = "euclidean", type = "C", method = "average", freq.option = 0, mem.eff = FALSE) {
  .Call(Rfast_k_nn,xnew, y, x, k, dist.type, type, method, freq.option, mem.eff)
}
