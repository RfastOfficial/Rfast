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
Version: \tab 2.0.6 \cr
Date: \tab 2022-02-15\cr
License: \tab GPL-2\cr
}
}

\author{
\itemize{
\item Manos Papadakis <papadakm95@gmail.com>
\item Michail Tsagris <mtsagris@uoc.gr>
\item Marios Dimitriadis <kmdimitriadis@gmail.com>
\item Stefanos Fafalios <stefanosfafalios@gmail.com>
\item Ioannis Tsamardinos <tsamard@csd.uoc.gr>
\item Matteo Fasiolo <matteo.fasiolo@gmail.com>
\item Giorgos Borboudakis <borbudak@gmail.com>
\item John Burkardt <jburkardt@fsu.edu>
\item Kleanthi Lakiotaki <kliolak@gmail.com>
\item Changliang Zou <nk.chlzou@gmail.com>
\item Christina Chatzipantsiou <chatzipantsiou@gmail.com>
}
}

\section{Maintainers}{
Manos Papadakis <papadakm95@gmail.com>
}

\note{
Acknowledgments: We would like to acknowledge:

\itemize{
\item Professor Kurt Hornik, Doctor Uwe Ligges (and the rest of R core team) for 
their invaluable help with this R package. 
\item Erich Studerus for his invaluable comments.
\item Neal Fultz for his suggestions. 
\item Vassilis Vasdekis for his invaluable help with the random effects models. 
\item Marios Dimitriadis' work was funded by the Special Account for Research Funds of the University of Crete, Department of Computer Science. 
\item Phillip Si is greatly acknowledged for his help with the Floyd-Warshal algorithm. 
\item Keefe Murphy for his invaluable help with NEWS file and for his suggestions. 
\item Zacharias Papadovassilakis gave us the inspiration for the memory efficient version of the k-NN algorithm. 
\item Yannis Pantazis explained us how the orhtogonal matching pursuit works. 
\item Achim Zeileis for his help with the quasi Poisson regression models. 
\item Pedro J. Aphalo for finding a bug.
\item Dimitris Kyriakis for finding a bug. 
\item Cyrille Conord for finding a bug.
\item Aaron Robotham for finding a bug.
\item Calvin Pan from the Department of Human Genetics at UCLA found a bug in the function glm_logistic and he is greatly acknowledged for that.
\item Adam Muschielok from Rodenstock GmbH found a bug in the function vmf.mle and he is greatly acknowledged for that.
\item Bret Presnell for detecting and correcting a bug in the function rvmf.
\item Gabriel Hoffman for spotting a mistake in the fuction dirimultinom.mle.
\item  Lutz von der Burchard for spotting a bug in the function bic.corfsreg.
}
From now on the Rfast can be used in C++ via LinkingTo mechanism.
\itemize{
\item The main namespace is "Rfast". Inside "Rfast" you will find two more namespaces, "vector" and "matrix".
\item Namespace "vector" for calling functions using an Rcpp's or RcppArmadillo's vector.
\item Namespace "matrix" for calling functions using an Rcpp's or RcppArmadillo's matrices.
\item The signatures of the functions and the arguments are the same that are exported in R.
}
For namespace "vector" the functions that are available are:
\enumerate{
	\item median(x)
	\item var(x, std = false, na_rm = false)
	\item mad(x, method = "median", na_rm = false)  
	\item shuffle(x,engine = Engine(time(0)) // Engine by default is default_random_engine. You can use anyone from C++.
}
For namespace "matrix" the functions that are available are:
\enumerate{
	\item transpose(x)
	\item matrix_multiplication(x,y)
	\item colSort(x, descend = false, stable = false, parallel = false)
	\item rowSort(x, descend = false, stable = false, parallel = false)
	\item is_symmetric(x)
	\item colMedian(x, na_rm = false, parallel = false)
	\item rowMedian(x, na_rm = false, parallel = false)
	\item colVars(x, std = false, na_rm = false, parallel = false)
	\item rowVars(x, std = false, na_rm = false, parallel = false)
	\item colMads(x, method = "median", na_rm = false, parallel = false)
	\item rowMads(x, method = "median", na_rm = false, parallel = false)
	\item colShuffle(x,engine = Engine(time(0)) // Engine by default is default_random_engine. You can use anyone from C++.
	\item rowShuffle(x,engine = Engine(time(0)) // Engine by default is default_random_engine. You can use anyone from C++.
}

How to use it:
\enumerate{
	\item Just add in "LinkingTo" in your NAMESPACE file the "Rfast" or in Rstudio "//[[Rcpp::depends(Rfast)]]".
	\item Include in your cpp files the header "Rfast.h" and enjoy!
}
}

