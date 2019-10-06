#[export]
sort_unique.length <- function(x) {
  .Call(Rfast_len_sort_unique_int,x)
}

#[export]
sort_unique <- function(x) {
  if(is.double(x)){
  	.Call(Rfast_sort_unique_double,x)
  }else{
  	.Call(Rfast_sort_unique_int,x)
  }
}