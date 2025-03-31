
r.files = "C:/Users/epapadakis/Documents/GitHub/Rfast2/R/"
rd.files = "C:/Users/epapadakis/Documents/GitHub/Rfast2/man/"
package.dir = "C:/Users/epapadakis/Documents/GitHub/Rfast2/"
nsp.dir = "C:/Users/epapadakis/Documents/GitHub/Rfast2/NAMESPACE"

r.files = "C:/Users/epapadakis/Documents/GitHub/Rfast/R/"
rd.files = "C:/Users/epapadakis/Documents/GitHub/Rfast/man/"
package.dir = "C:/Users/epapadakis/Documents/GitHub/Rfast/"
nsp.dir = "C:/Users/epapadakis/Documents/GitHub/Rfast/NAMESPACE"


r.files = "C:/Users/epapadakis/Documents/GitHub/Rnanoflann/R/"
rd.files = "C:/Users/epapadakis/Documents/GitHub/Rnanoflann/man/"
package.dir = "C:/Users/epapadakis/Documents/GitHub/Rnanoflann/"
nsp.dir = "C:/Users/epapadakis/Documents/GitHub/Rnanoflann/NAMESPACE"

Rfast::AddToNamespace(nsp.dir,r.files)
Rfast::checkAliases(rd.files,r.files)
Rfast::checkUsage(rd.files,r.files)
Rfast::checkExamples(rd.files,"Rfast",print.names = TRUE)
Rfast::checkTF(rd.files)
Rfast::sourceR(r.files,print.errors = TRUE)
Rfast::sourceRd(rd.files,print.errors = TRUE)

Rfast2::benchmark(Rfast::checkExamples(rd.files,print.errors = "a.txt"),checkExamples2(rd.files,print.errors = "b.txt"),times = 10)
Rfast2::benchmark(Rfast::sourceRd(rd.files),sourceRd2(rd.files),times = 10)

checkExamples2<-function(path.man,each = 1,print.errors = stderr(),print.names = FALSE){
  examples_files <- .Call("Rfast_read_examples",path.man)
  error_files<-vector("character")
  examples <- examples_files$examples
  file_names<-examples_files$files
  if(!is.null(print.errors)){
    warning_error_function <-function(err){
      write(paste(file_names[i],":","\n",err),print.errors)
      error_files <<- c(error_files,file_names[i])
    }
  }else{
    warning_error_function <-function(err){
      error_files <<- c(error_files,file_names[i])
    }
  }
  if(print.names){
    foreach::foreach(i=1:length(examples)) %dopar% {
      print(file_names[i])
      for(j in 1:each){
        tryCatch(eval(parse(text=examples[i])),error=function(err){
          write(paste(file_names[i],":","\n",err),print.errors)
          error_files <<- c(error_files,file_names[i])
        }, warning=function(err){
          write(paste(file_names[i],":","\n",err),print.errors)
          error_files <<- c(error_files,file_names[i])
        })
      }
    }
  }else{
    foreach::foreach(i=1:length(examples)) %dopar% {
      for(j in 1:each){
        tryCatch(eval(parse(text=examples[i])),error=function(err){
          write(paste(file_names[i],":","\n",err),print.errors)
          error_files <<- c(error_files,file_names[i])
        }, warning=function(err){
          write(paste(file_names[i],":","\n",err),print.errors)
          error_files <<- c(error_files,file_names[i])
        })
      }
    }
  }
  list("Errors"=error_files,"Big Examples"=examples_files$long_lines,"dont read"=examples_files$`dont read`)
}

#[export]
sourceRd2 <- function(path,print.errors=FALSE) {
  file_names <- Rfast::read.directory(path)
  error_files<-vector("character")
  if(print.errors){
    warning_error_function <-function(err){
      write(paste(file_names[i],":","\n",err),stderr())
      error_files <<- c(error_files,file_names[i])
    }
  }else{
    warning_error_function <-function(err){
      error_files <<- c(error_files,file_names[i])
    }
  }
  error<-0
  foreach(i=1:length(file_names)) %dopar% {
    error<-tools::checkRd(sprintf("%s%s",path,file_names[i]))
    if(length(error)!=0){
      warning_error_function(error)
    }
  }
  if(length(error_files)==0){
    return("Everything is ok..!")
  }
  return(error_files)
}



y=as.character(Rfast2::Runif(10^7))
x=sample(y,10^6)
Rfast2::benchmark(x %in% y,myHash(x,y,"fnv1a"),myHash(x,y,"pengy"),myHash(x,y,"superfast"),myHash(x,y,"murmur"),
                  myHash_rcpp(x,y,"fnv1a"),myHash_rcpp(x,y,"pengy"),myHash_rcpp(x,y,"superfast"),myHash_rcpp(x,y,"murmur"),
                  myHash_rcpp_list(x,y,"fnv1a"),myHash_rcpp_list(x,y,"pengy"),myHash_rcpp_list(x,y,"superfast"),myHash_rcpp_list(x,y,"murmur"),times = 10)

Rfast2::benchmark(myHash_table(y,"cpp"),myHash_table(y,"fnv1a"),myHash_table(y,"pengy"),myHash_table(y,"superfast"),myHash_table(y,"murmur"),times = 10)
Rfast2::benchmark(myHash_table(y,"cpp"),myHash_table(y,"fnv1a"),myHash_table(y,"pengy"),myHash_table(y,"superfast"),myHash_table(y,"murmur"),myHash_table_rcpp(y,"fnv1a"),myHash_table_rcpp(y,"pengy"),myHash_table_rcpp(y,"superfast"),myHash_table_rcpp(y,"murmur"),times = 10)





# Base R boosted by Rfast
mostfreqval5 <- function(x,k=1){
  x %>% Rfast::Table() %>% Rfast::nth(k=k,num.of.nths = k,descending = T) %>% names() %>% head(k)
}
