#[export]
all_equals<-function(x,y,round_digits = FALSE,without_attr=FALSE,fast_result=FALSE){
	if(!fast_result){
		dmx <- dim(x) ; dmy <- dim(y)
		if(is.null(dmx) && is.null(dmy)){ #if are vectors/factors
			if(length(x)!=length(y)){
				return("mismatch: length")
			}
		}else{ #matrix/dataframe
			if(!identical(dmx,dmy)){
				return("mismatch: dimension")
			}
		}
		if(without_attr){
			attributes(x)<-attributes(y)<-NULL
		}else{
			the_attrs <- identical(attributes(x),attributes(y))
			if(!the_attrs){
				return ("mismatch: attributes")
			}
		}
	}
	if(round_digits!=FALSE){
		x<-Rfast::Round(x,round_digits)
		y<-Rfast::Round(y,round_digits)
	}
	identical(x,y)
}