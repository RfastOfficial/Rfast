
#[export]
hash.list <- function(key,x) {
  .Call(Rfast_Hash_list,key,x)
}

#[export]
Hash.key.multi<-function(x,...,sep = " "){
    x<-.Call(Rfast_Hash_key_multi,x,paste(...,collapse = sep),sep)
    if(x == "") NULL else x
}

#[export]
hash.find <- function(x,key) {
  .Call(Rfast_hash_find,x,key)
}

#[export]
hash2list <- function(x,sorting=FALSE) {
  .Call(Rfast_hash2list,x,sorting)
}

#[export]
Hash<-function(keys=NULL,values=NULL){
	keys <- if(is.character(keys)) keys else as.character(keys)
	i <- 1
	len <- length(keys)
	x <- new.env(size=len)
	x$.length <- len
	for(key in keys){
		x[[key]] <- values[i]
		i <- i + 1
	}
	class(x)<-"Hash"
	x
}