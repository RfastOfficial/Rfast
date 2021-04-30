
#[export]
group.sum <- function(x, ina,ina.max = NULL,ina.min = NULL) {
	#.Defunct("group(x,ina,method=\"sum\")","Rfast")
	group(x,ina,ina.min=ina.min,ina.max=ina.max)
}


##[export]
#group.min_max <- function(x, ina,ina.max = max(ina)) {
#	.Defunct("group(x,ina,method=\"min.max\")","Rfast")
#}
#
#
##[export]
#group.min <- function(x, ina,ina.max = max(ina)) {
#	.Defunct("group(x,ina,method=\"min\")","Rfast")
#}
#
#
##[export]
#group.med <- function(x,ina,ina.length.unique=length(unique(ina))) {
#	.Defunct("group(x,ina,method=\"med\")","Rfast")
#}
#
#[export]
group.mean <- function(x, ina,ina.max = max(ina)) {
	#.Defunct("group(x,ina,method=\"mean\")","Rfast")
	group(x,ina,ina.max=ina.max)
}
#
##[export]
#group.max <- function(x, ina,ina.min = NULL,ina.max = NULL){
#	.Defunct("group(x,ina,method=\"max\")","Rfast")
#}
#
##[export]
#group.mad <- function(x,ina,method = "median") {
#	.Defunct("group(x,ina,method=\"mad\",mad.method=\"median or mean\")","Rfast")
#}
#
##[export]
#group.any <- function(x, ina,ina.max = max(ina)) {
#	.Defunct("group(x,ina,method=\"any\")","Rfast")
#}
#
##[export]
#group.all <- function(x, ina,ina.max = max(ina)) {
#	.Defunct("group(x,ina,method=\"all\")","Rfast")
#}
#
##[export]
#group.var <- function(x, ina,ina.max = max(ina)) {
#	.Defunct("group(x,ina,method=\"var\")","Rfast")
#}

#[export]
group<-function(x,ina,method="sum",ina.min=NULL,ina.max = NULL,ina.length.unique=NULL,mad.method="median") {
	if(method=="med"){
		ina.max<- if(is.null(ina.max)) max(ina) 
		ina.min<- if(is.null(ina.length.unique)) length(unique(ina))
	}
	.Call(Rfast_group,x,ina,method,ina.min,ina.max,mad.method)
}