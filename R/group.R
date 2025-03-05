
#[export]
group.sum <- function(x, ina,ina.max = NULL,ina.min = NULL) {
	.Defunct("group(x,ina,method=\"sum\")","Rfast")
}

#[export]
group<-function(x,ina,method="sum",ina.min=NULL,ina.max = NULL,
	ina.length.unique=NULL,mad.method="median", std = FALSE, sorted = TRUE) {
	# if(method=="med"){
	# 	ina.max<- if(is.null(ina.max)) max(ina) 
	# 	ina.min<- if(is.null(ina.length.unique)) length(unique(ina))
	# }
	.Call(Rfast_group,x,ina,method,mad.method, std, sorted)
}