#[export]
iterator<-function(x,method="ceil",type="vector",by=1){
    if((type=="vector" || type=="list") && method!="ceil"){
        stop("Vector myst have method ceil.")
    }
    y<-new.env()
    class(y)<-"iterator"
    var <- if(class(x)=="iterator") x$.variable else x
    y$.variable<-var
    y$.method <- method
    y$.value <- 1
    y$.by <- by
    y$.type <- type
    y$copy<-function(){
        newy<-iterator(y,y$.method,y$.type,y$.by)
        newy$value <- y$value
        newy
    }
    y$end<-function(){
        endy<-iterator(y)
        if(endy$.method=="ceil"){
            endy$.value<-length(endy$.variable) + 1
        }else if(endy$.type=="matrix" || endy$.type=="data.frame"){
            endy$.value<-dim(endy$.variable)[(if(endy$.method=="col") 2 else 1)] + 1
        }
        endy
    }
    y$nextElem<-function(){
        y$.value <- y$.value + y$.by
    }
    y$prevElem<-function(){
        y$.value <- y$.value - y$.by
    }
    y
}