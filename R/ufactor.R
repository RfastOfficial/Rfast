
#[export]
ufactor<-function(x){
    y<-new.env()
    un<- if(is.character(x)) Rfast::Sort(unique(x)) else Rfast::sort_unique(x)
    y$values<-Rfast::Match(x,un)
    y$levels<-un
    class(y)<-"ufactor"
    lockEnvironment(y)
    lockBinding("values",y)
    lockBinding("levels",y)
    y
}

