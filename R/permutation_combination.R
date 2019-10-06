#[export]
bincomb <- function(n) {
  .Call(Rfast_bincomb,n)
}

#[export]
rowShuffle <- function(x) {
	.Call(Rfast_row_shuffle,x)
}

#[export]
colShuffle <- function(x) {
	.Call(Rfast_col_shuffle,x)
}

#[export]
permutation.prev <- function(x,nperm=gamma(length(x)+1)) {
  .Call(Rfast_permutation_prev,x,nperm)
}

#[export]
permutation.next <- function(x,nperm=gamma(length(x)+1)) {
  .Call(Rfast_permutation_next,x,nperm)
}

#[export]
permutation <- function(x,nperm=gamma(length(x)+1)) {
  .Call(Rfast_permutation,x,nperm)
}

#[export]
comb_n <- function(n,k,simplify=TRUE) {
	if(k<0){
	  	stop("K must be a positive number.")
	}
	if(length(n)==1) {
		n<-if(n<0) n:-1 else 1:n
	}
	.Call(Rfast_comb_n,n,k,simplify)
}