
#[export]
env.copy<-function(x,all.names=FALSE){
    y<-new.env()
    all.vars<-ls(x,all.names = all.names)
    for(var in all.vars){
        val<-x[[var]]
        y[[var]] <- if(is.environment(val)) env.copy(val,all.names) else x[[var]]
    }
    y
}