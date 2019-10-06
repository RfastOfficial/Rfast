#[export]
eigen.sym <- function(A, k, vectors = TRUE) {
  .Call(Rfast_eigs_sym_c,A,k,vectors)
}