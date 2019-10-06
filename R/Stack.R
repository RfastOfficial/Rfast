#[export]
Stack<-function(x,type=NULL){
    y<-new.env()
    class(y)<-"Stack"
    var<- if(class(x)=="Stack") y$.variable else x
    y$.variable<-if(is.null(type)){ y$.total<-length(var);var }else{ y$.total<-x;type(x)}
    y$.top_element<-y$.total + 1
    y$pop = function(){ 
        if(y$.top_element>0){
            y$.top_element<-y$.top_element-1
            y$.variable[y$.top_element]
        }else{
            stop("Stack is empty. You have to fill it first.")
        }
    }
    y$top = function(){
        if(y$.top_element>0){
            y$.variable[y$.top_element-1]
        }else{
            stop("Stack is empty. You have to fill it first.")
        }
    }
    y$push = function(x){
        if(y$.top_element>y$.total){
            y$.variable=c(y$.variable,x)
        	y$.total=y$.total+1
        }else{
            y$.variable[y$.top_element]=x
        }
        y$.top_element=y$.top_element+1
    }
    y$clear = function(){
    	y$.top_element<-1
    }
    lockEnvironment(y)
    y
}