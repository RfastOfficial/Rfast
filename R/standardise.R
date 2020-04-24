#[export]
standardise <- function (x, center = TRUE, scale = TRUE) {
  if ( center & scale ) {
    y <- t(x) - Rfast::colmeans(x)
    y <- y / sqrt( Rfast::rowsums(y^2) ) * sqrt( (dim(x)[1] - 1) )
	y <- t(y)
  } else if ( center & !scale ) {
    m <- Rfast::colmeans(x)
    y <- eachrow(x, m, oper ="-" )
  } else if ( !center & scale ) {
    s <- Rfast::colVars(x, std = TRUE)
    y <- eachrow(x, s, oper = "/")
  } else {
    y <- x
  }
}
