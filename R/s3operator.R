##[dont read]

delete <- new.env(size=1)
class(delete)<-"delete"
lockEnvironment(delete,TRUE)

#[export s3]
"[<-.Hash"<-function(x,...,sep = " ",value){
    key<-paste0(c(...),collapse = sep)
    if(is.null(x[[key]])){
        if(class(value) == "delete"){
            return(x)
        }else{
            x$.length <- x$.length + 1 
        }
    }else if(class(value) == "delete"){
        x$.length <- x$.length - 1 
    }
    x[[key]]<-value
    x
}
#[export s3]
"[.Hash"<-function(x,...,sep = " "){
    key<-paste0(c(...),collapse = sep)
    val <- x[[key]]
    if(class(val) == "delete"){
        NULL
    }else{
        val
    }
}
#[export s3]
length.Hash<-function(x){
    x$.length
}


#[export s3]
print.Hash<-function(x,...){
	len <- length(x)
    cat("<Hash> contains ",len," keys-values pairs\n")
    is.Ctype <- function(t){
    	return(is.character(t) || is.numeric(t) || is.double(t) || is.integer(t))
    }
    if(len > 0) {
    	nam <- names(x)
        nam <- nam[nam != ".length"]
        size_of_each_name <- sapply(nam,nchar)
        max_size_of_keys<-max(size_of_each_name)
        for(i in 1:length(nam)){
            key <- nam[i]
            val <-  x[[key]]
            if(class(val) == "delete"){
            	next
            }
            space <- paste0(rep(" ",max_size_of_keys-size_of_each_name[i]),collapse = "")
            vals <- if(is.function(val)){
            	paste0(class(val),"(",length(formals(val)),")",collapse = "")
            }else if(!is.null(dim(val))){
            	paste0(class(val),"(",nrow(val),",",ncol(val),")",collapse = "")
            }else if(is.environment(val) || is.list(val) || is.object(val)){
            	paste0(class(val),"(",length(val),")",collapse = "")
            }else if(is.Ctype(val)){
            	if(length(val)>1){
            		paste0(class(val),"(",length(val),")",collapse = "")
            	}
            	else{
            		val
            	}
            }else if(is.null(val)){
            	"NULL"
            }else{
            	"Error: uknown type of value."
            }
            cat(key,space,": ",vals,"\n")
        }
    }
}


##########################################################################################################################
#[export s3]
print.environment<-function(x,all.names=FALSE,...){
	.print_helper_for_environment<-function(x,all.names,count_depth){
	    all.vars<-ls(x,all.names = all.names)
	    for(var in all.vars){
	    	spaces <- rep("  ",count_depth)
	        cat(spaces)
	        cat(var,"= ")
	        val<-x[[var]]
	        if(is.environment(val)){
	            len<-length(val)
	            cat(sprintf("environment(%s){",len))
	            if(len>0){
	                cat("\n")
	                .print_helper_for_environment(val,all.names,count_depth+1)
	            }
	            if(count_depth > 0)
	            	cat(rep(" ",count_depth),"}\n")
	            else
	            	cat("}\n")
	        }else{
	            cl<-class(val)
	            cat(cl)
	            len<-if(is.matrix(val) || is.data.frame(val)){ dm<-dim(val);sprintf("%s,%s",dm[1],dm[2]) }
	            else if(is.function(val)) length(formals(val)) else length(val)
	            cat(sprintf("(%s)\n",len))
	        }
	    }
	}
    .print_helper_for_environment(x,all.names,0)
}

##########################################################################################################################

#[export]
"Elem<-"<-function(x,value){
    UseMethod("Elem<-")
}
#[export s3]
"Elem<-.iterator"<-function(x,value){
    if(x$.method=="ceil"){
        x[[".variable"]][x[[".value"]]] <- value
    }else if(x$.method=="col"){
        x[[".variable"]][,x[[".value"]]] <- value
    }else if(x$.method=="row"){
        x[[".variable"]][x[[".value"]],] <- value
    }else{
        stop("Error...who knows...")
    }
    x
}
#[export]
Elem<-function(x){
    UseMethod("Elem")
}
#[export s3]
Elem.iterator<-function(x){
    if(x$.method=="ceil"){
        x$.variable[x$.value]
    }else if(x$.method=="col"){
        x$.variable[,x$.value]
    }else if(x$.method=="row"){
        x$.variable[x$.value,]
    }else{
        stop("Error...who knows...")
    }
}
#[export s3]
print.iterator<-function(x,...){
    cat("variable:\n     ")
    if(is.null(dim(x$.variable))){
        cat("length: ",length(x$.variable),"\n")
    }else{
        cat("dimensions: ",dim(x$.variable),"\n")
    }
    cat("     class : ",class(x$.variable),"\n")
    cat("method: ",x$.method,"\n")
    cat("by    : ",x$.by,"\n")
    cat("value : ",x$.value,"\n")
    cat("type  : ",x$.type,"\n")       
}
#[export s3]
"==.iterator"<-function(x,y){
    if(class(y)=="iterator"){
        identical(y$.variable,x$.variable) && identical(y$.by,x$.by) && 
        identical(y$.type,x$.type) && identical(y$.method,x$.method) && identical(y$.value,x$.value)
    }else{
        stop("Error in function '==', second argument is not class 'iterator'.")
    }
}
#[export s3]
"!=.iterator"<-function(x,y){
    !(x==y)
}


##########################################################################################################################
#[export s3]
"[.ufactor"<-function(x,i){
	x$levels[x$values[i]]
}
#[export s3]
print.ufactor<-function(x,...){
    cat("Values: \n ")
    options(digits = 15)
    cat(x$levels[x$values])
    cat("\nLevels: \n ",x$levels)
}

##########################################################################################################################



##########################################################################################################################



##########################################################################################################################