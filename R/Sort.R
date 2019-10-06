#[export]
Sort <- function(x,descending=FALSE,partial=NULL,stable=FALSE,na.last=NULL) {
	if(stable){
		.Call(Rfast_stable_sort,x,descending)
	}
	else if(!is.null(partial)){
	 	.Call(Rfast_partial_sort,x,partial,descending)
	}else if(is.character(x)){
		.Call(Rfast_sort_string,x,descending)
	}else{
		if(identical(na.last,FALSE)){
			.Call(Rfast_Sort_na_first,x,descending)
		}else{
			.Call(Rfast_Sort,x,descending,na.last)
		}
	}
}

#[export]
Sort.int <- function(x) {
  .Call(Rfast_sort_int,x)
}

#[export]
rowSort <- function(x,descending=FALSE,stable=FALSE,parallel=FALSE) {
	.Call(Rfast_sort_mat,x,descending,TRUE,stable,parallel)
}

#[export]
colSort <- function(x,descending=FALSE,stable=FALSE,parallel=FALSE) {
	.Call(Rfast_sort_mat,x,descending,FALSE,stable,parallel)
}

##[export]
#sort_mat <- function(x,by.row=FALSE,descending=FALSE,stable=FALSE,parallel=FALSE) {
#	.Defunct(if(by.row) "Rfast::rowSort" else "Rfast::colSort","Rfast")
#}

#[export]
sort_cor_vectors <- function(x, base, stable = FALSE, descending = FALSE) {
  x[Rfast::Order(base,stable,descending)]
}