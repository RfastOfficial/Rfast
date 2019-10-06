%[dont read]

\name{Rfast-package}
\alias{Rfast-package}
\docType{package}
\title{
Really fast R functions
}

\description{
A collection of Rfast functions for data analysis. 
Note 1: The vast majority of the functions accept matrices only, not data.frames. 
Note 2: Do not have matrices or vectors with have missing data (i.e NAs). We do no check about them and C++ internally transforms them into zeros (0), so you may get wrong results.
Note 3: In general, make sure you give the correct input, in order to get the correct output. We do no checks and this is one of the many reasons we are fast.   
}

\details{
\tabular{ll}{
Package: \tab Rfast\cr
Type: \tab Package\cr
Version: \tab 1.9.6 \cr
Date: \tab 2019-09-30\cr
License: \tab GPL-2\cr
}
}

\author{
Manos Papadakis <papadakm95@gmail.com>, Michail Tsagris <mtsagris@yahoo.gr>, Marios Dimitriadis <kmdimitriadis@gmail.com>, Stefanos Fafalios <stefanosfafalios@gmail.com>, 
Ioannis Tsamardinos <tsamard@csd.uoc.gr>, Matteo Fasiolo <matteo.fasiolo@gmail.com>, Giorgos Borboudakis <borbudak@gmail.com>, John Burkardt <jburkardt@fsu.edu> and 
Kleanthi Lakiotaki <kliolak@gmail.com>, Changliang Zou <nk.chlzou@gmail.com>, Christina Chatzipantsiou <chatzipantsiou@gmail.com>
}

\section{Maintainers}{
Manos Papadakis <papadakm95@gmail.com>
}

\note{
Acknowledgments:

We would like to acknowledge Professor Kurt Hornik, Doctor Uwe Ligges (and the rest of R core team) for 
their invaluable help with this R package. 
Erich Studerus for his invaluable comments.
Neal Fultz for his suggestions. 
Vassilis Vasdekis for his invaluable help with the random effects models. 
Marios Dimitriadis work was funded by the Special Account for Research Funds of the University of Crete, Department of Computer Science. 
Phillip Si is greatly acknowledged for his help with the Floyd-Warshal algorithm. 
Keefe Murphy for his invaluable help with NEWS file and for his suggestions. 
Zacharias Papadovassilakis gave us the inspiration for the memory efficient version of the k-NN algorithm. Yannis Pantazis explained us how the orhtogonal matching pursuit works. 
Achim Zeileis for his help with the quasi Poisson regression models. 
Pedro J. Aphalo for finding a bug.
Dimitris Kyriakis for finding a bug. 
Cyrille Conord for finding a bug.
Aaron Robotham for finding a bug.
Calvin Pan from the Department of Human Genetics at UCLA found a bug in the function glm_logistic and he is greatly acknowledged for that.
Adam Muschielok from Rodenstock GmbH found a bug in the function vmf.mle and he is greatly acknowledged for that.
Calvin Pan from the Department of Human Genetics at UCLA found a bug in the function 
glm_logistic and he is greatly acknowledged for that.
Adam Muschielok from Rodenstock GmbH found a bug in the function vmf.mle and he is 
greatly acknowledged for that.


From now on the Rfast can be used in C++ via LinkingTo mechanism.
The main namespace is "Rfast". Inside "Rfast" you will find two more namespaces, "vector" and "matrix".
Namespace "vector" for calling functions using an Rcpp's or RcppArmadillo's vector.
Namespace "matrix" for calling functions using an Rcpp's or RcppArmadillo's matrices.
The signatures of the functions and the arguments are the same that are exported in R.
For namespace "vector" the functions that are available are:
  
	median(x)
	var(x, std = false, na_rm = false)
	mad(x, method = "median", na_rm = false)  
	shuffle(x,engine = Engine(time(0)) // Engine by default is default_random_engine. You can use anyone from C++.

For namespace "matrix" the functions that are available are:
  
	transpose(x)
	matrix_multiplication(x,y)
	colSort(x, descend = false, stable = false, parallel = false)
	rowSort(x, descend = false, stable = false, parallel = false)
	is_symmetric(x)
	colMedian(x, na_rm = false, parallel = false)
	rowMedian(x, na_rm = false, parallel = false)
	colVars(x, std = false, na_rm = false, parallel = false)
	rowVars(x, std = false, na_rm = false, parallel = false)
	colMads(x, method = "median", na_rm = false, parallel = false)
	rowMads(x, method = "median", na_rm = false, parallel = false)
	colShuffle(x,engine = Engine(time(0)) // Engine by default is default_random_engine. You can use anyone from C++.
	rowShuffle(x,engine = Engine(time(0)) // Engine by default is default_random_engine. You can use anyone from C++.
	
How to use it:
	1. Just add in "LinkingTo" in your NAMESPACE file the "Rfast" or in Rstudio "//[[Rcpp::depends(Rfast)]]".
	2. Include in your cpp files the header "Rfast.h" and enjoy!
}

