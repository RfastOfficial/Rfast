#[export]
dista <- function(xnew, x, type = "euclidean", k = 0, index = FALSE, trans = TRUE, square = FALSE, p = 0, result = "matrix", parallel = FALSE) {
    if (type == "canberra1" || type == "canberra2") {
        .Deprecated("The replacement type is \"canberra\"", "Rfast")
    }
    if(type != "haversine"){
        xnew <- t(xnew)
        x <- t(x)
    }
    if(result == "matrix"){
		.Call(Rfast_dista, xnew, x, type, square, p, k, index, parallel)
	}else if(result == "sum"){
  		.Call(Rfast_total_dista, xnew, x, type, square, p, k, parallel)
	}
    if (trans) x <- t(x)
    x
}

#[export]
total.dista <- function(xnew, x, type = "euclidean", k = 0, square = FALSE, p = 0, parallel = FALSE) {
    .Deprecated("Use dista(x, result = \"sum\")", "Rfast")
    if(type != "haversine"){
        xnew <- t(xnew)
        x <- t(x)
    }
    .Call(Rfast_total_dista, xnew, x, type, square, p, k, parallel)
}