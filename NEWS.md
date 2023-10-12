<h1 align="center">Rfast</h1>


### **Version 2.0.9 - Hades**
***
> <u> **New** </u>
>
>  |	    Function	    |           What's new!             |
>  | -------------------  | ------------------------------- |
>  |  coeff				| Coefficient matrix.  |
>
>
> <u> **Improved** </u>(_**by speed, correctness or options**_) 
>
>  |	    Function	    |           What's new!             |
>  | -------------------  | ------------------------------- |
>  |  cor.fsreg | Bug fix in the adjusted R-squared computation. |
>  |  ompr | Bug fix in the adjusted R-squared computation. |
>  |  rmdp | Support parallel. |
>  |  colMads				| Supports data.frame  |
>  |  eachrow | Fix bug.  |
>  |  colMedians | Fix bug about data.frame.  |
>  |  rowMedians | Fix bug.  |
>  |  colVars | Fix bug about data.frame.  |
>  |  Dist(...,method = "minimum") | Bug fix.  |
>  |  Dist(...,method = "bhattacharyya") | Bug fix .  |
>  |  Dist(...,method = "itakura_saito") | Bug fix. Supports parallel.  |
>  |  Dist(...,method = "canberra1") | Depricated.  |
>  |  Dist(...,method = "canberra2") | Depricated.  |
>  |  Dist(...,method = "jensen_shannon") | Supports parallel.  |
>  |  Dist(...,method = "kullback_leibler") | Supports parallel.  |
>  |  Dist | New methods ("canberra","chi_square","soergel","sorensen","cosine","wave_hedges","motyka","harmonic_mean","jeffries_matusita","gower","kulczynski")  |
>  |  dista | All the methods from the Dist function are now supported.  |
>  |  total.dist  | All the methods from the Dist function are now supported.  |
>  |  total.dista | All the methods from the Dist function are now supported.  |
>  |  Sort | Supports parallel. Available only in Windows and most unix.  |
>  |  Order | Supports parallel. Available only in Windows and most unix.  |
>  |  Rank | Supports parallel. Available only in Windows and most unix. Add extra option for stable algorithm.  |
>  |  checkUsage | checking if usage lines are wider than 90 characters.  |
>
>
> <u> **LinkingTo** </u>(_**by speed, correctness or options**_) 
>
>  |	    Function/Structure	    |           What's new!             |
>  | -------------------  | ------------------------------- |
>  |  colMads				| Supports data.frame  |
>


</br>
</br>

### **Version 2.0.8 - Poseidon**
***
> <u> **Improved** </u>(_**by speed, correctness or options**_) 
>
>  |	    Function	    |           What's new!             |
>  | -------------------  | ------------------------------- |
>  |  colMads				| Supports data.frame  |
>  |  colShuffle				| Supports data.frame  |
>  |  colMedian				| Supports data.frame  |
>  |  colRanks				| Supports data.frame and fix bug about method "min".  |
>  |  colSort				| Supports data.frame  |
>  |  colMaxs				| Supports data.frame  |
>  |  colMins				| Supports data.frame  |
>  |  colMinsMaxs				| Supports data.frame and add extra option for parallelism.  |
>  |  colrange				| Supports data.frame and add extra option for parallelism.  |
>  |  checkExamples				| Add extra option for running examples inside environment. Treat special functions with new export attribute. |
>  |  sort_mat              | Depricated. |
>
>
> <u> **LinkingTo** </u>(_**by speed, correctness or options**_) 
>
>  |	    Function/Structure	    |           What's new!             |
>  | -------------------  | ------------------------------- |
>  |  colMads				| Supports data.frame  |
>  |  colShuffle				| Supports data.frame  |
>  |  colMedian				| Supports data.frame  |
>  |  colSort				| Supports data.frame  |
>  |  colMaxs				| Supports data.frame  |
>  |  colMins				| Supports data.frame  |
>  |  colMinsMaxs				| Supports data.frame and add extra option for parallelism.  |
>  |  NA				| A structure for manage R's NA values. Supports {int, double, bool, string}.  |
>  |  Type				| A structure for manage R's types. Supports {REAL,INT,CHAR,FACTOR,LIST,DATAFRAME}.  |
>
> ### **Comments**
***
> From now on all the function will be inside namespace Rfast. Sub-namespaces "vector" and "matrix" have been removed.
>
>


</br>
</br>

### **Version 2.0.7 - Zeus**
***
> <u> **Improved** </u>(_**by speed, correctness or options**_) 
>
>  |	    Function	    |           What's new!             |
>  | -------------------  | ------------------------------- |
>  |  colVars				| Supports data.frame  |
>  |  colmeans				| Supports data.frame  |
>  |  invdir.mle			| A small optimization |


</br>
</br>

### **Version 2.0.6**
***
> <u> **Improved** </u>(_**by speed, correctness or options**_) 
>
>  |	    Function	    |           What's new!             |
>  | -------------------  | ------------------------------- |
>  |  Dist(method = "bhattacharyya") | Minor optimization. |


</br>
</br>

### **Version 2.0.5**
***
> <u> **Improved** </u>(_**by speed, correctness or options**_) 
>
>  |	    Function	    |           What's new!             |
>  | -------------------  | ------------------------------- |
>  |  Dist, total.dist, vecdist       | New distance algorithm named "haversine". |
>  |  gaussian.nb           | 3 times faster than the old version. Add new option for parallel which decreases execution 18 times. |
>  |  checkUsage            | Fix a bug about "=function". Optimized empty return values and error messages. |
>  |  checkAliases, checkTF | Optimized empty return values. |
>  |  edist                 | Fix a bug. |
>  |  dirimultinom.mle      | Fix a bug. |


</br>
</br>

### **Version 2.0.4**
***
> <u> **Improved** </u>(_**by speed, correctness or options**_) 
>
>  |	    Function	    |           What's new!             |
>  | ------------------- | ------------------------------- |
>  |         nth           | Fix a bug using na_rm. In github is #52. |
>  |  colsums,rowsums      | Add new option for removing NAs and optimize option indices. |

</br>
</br>

### **Version 2.0.0**
***
> <u> **Improved** </u>(_**by speed, correctness or options**_) 
>
>  |	    Function	    |           What's new!             |
>  | ------------------- | ------------------------------- |
>  |         Var           | Time improvement for removing NAs.|
>  |  checkUsage           | From now on check for S3methods and fix bug.        |
>  |  Tcrossprod           | Same as R's tcrossprod.        |
>  |  Crossprod            | Same as R's crossprod.        |
>  |  matrix_multiplication | Add new arguments for perfoming cross or tcross product. |

</br>
</br>

###  **Version 1.9.9**
***
> <u> **Improved** </u>(_**by speed, correctness or options**_) 
>
>  |	    Function	    |           What's new!             |
>  | ------------------- | ------------------------------- |
>  | Var                   | Time improvement for removing NAs.|
>  | colVars               | Time improvement and depricate argument "suma".|
>  | rowVars               | Time improvement and depricate argument "suma".|
>  | rowMedians            | Add option for removing NAs if any.|
>  | colMads               | Add option for removing NAs if any, add option for choosing method, fix bug and time improvement.|
>  | rowMads               | Add option for removing NAs if any, add option for choosing method, fix bug and time improvement.|
>  | mad2                  | Deprecated. Use "Mad" instead.|
>  | med                   | Deprecated. Use "Median" instead.|
>  | colShuffle            | Fix a bug.|
>  | rowShuffle            | Fix a bug.|
>  | rmdp                  | Fix a bug.|
>  |  checkUsage           | From now on check for S3methods.        |
>  |  AddToNamespace       | From now on export S3methods.           |
>  |  checkAliases         | From now on check the S3methods.        |
>  |  checkNamespace       | Don't use it for now.             |
>
>
> ### **Comments**
***
> From now on the Rfast can be used in C++ via "LinkingTo" mechanism.
> The main namespace is "Rfast". Inside "Rfast" you will find two more namespaces, "vector" and "matrix".
> Namespace "vector" for calling functions using an Rcpp's or RcppArmadillo's vector.
> Namespace "matrix" for calling functions using an Rcpp's or RcppArmadillo's matrices.
> The signatures of the functions and the arguments are the same that are exported in R.
>
> For namespace "vector" the functions that are available are
>   
> - median(x)
> - var(x, std = false, na_rm = false)
> - mad(x, method = "median", na_rm = false)
> - shuffle(x,engine = Engine(time(0)) // Engine by default is default_random_engine. You can use anyone from C++.
>   
> For namespace "matrix" the functions that are available are
>   
> * transpose(x)
> * matrix_multiplication(x,y)
> * colSort(x, descend = false, stable = false, parallel = false)
> * rowSort(x, descend = false, stable = false, parallel = false)
> * is_symmetric(x)
> * colMedian(x, na_rm = false, parallel = false)
> * rowMedian(x, na_rm = false, parallel = false)
> * colVars(x, std = false, na_rm = false, parallel = false)
> * rowVars(x, std = false, na_rm = false, parallel = false)
> * colMads(x, method = "median", na_rm = false, parallel = false)
> * rowMads(x, method = "median", na_rm = false, parallel = false)
> * colShuffle(x,engine = Engine(time(0)) // Engine by default is default_random_engine. You can use anyone from C++.
> * rowShuffle(x,engine = Engine(time(0)) // Engine by default is default_random_engine. You can use anyone from C++.
>
> How to use it
>
> 1. Just add in "LinkingTo" in your NAMESPACE file the "Rfast" or in Rstudio add in the file "//[[Rcppdepends(Rfast)]]".
>
>  
>  Include in your cpp files the header "Rfast.h" and enjoy!

</br>
</br>

###  **Version 1.9.8**
***
> <u> **Improved** </u>(_**by speed, correctness or options**_) 
>
>  |	    Function	    |           What's new!             |
>  | ------------------- | ------------------------------- |
>  | spml.mle              |         Fix of a bug.             |
>  | dirknn                |         Time improvement.         |
>  | glm_logistic          |  Made the code more robust.       |
>  | vmf.mle               |  Made the code more robust.       |
>  | cor.fbed              |         Fix of a bug.             |
>  | Dist                  |         Fix of a bug.             |

</br>
</br>

###  **Version 1.9.4**
***
>
> <u> **Improved** </u>(_**by speed, correctness or options**_) 
>
>  |	    Function	    |           What's new!             |
>  | ------------------- | ------------------------------- |
>  | allbetas              | Time improvement |
>  | cor.fbed              | Time improvement |
>  | cor.fsreg             | Time improvement |
>  | omp                   | Time improvement |
>  | ompr                  | Time improvement |
>  | score.betaregs        | Time improvement |
>  | score.gammaregs       | Time improvement |
>  | score.glms            | Time improvement |
>  | score.invgaussregs    | Time improvement |
>  | score.negbinregs      | Time improvement |
>  | score.ztpregs         | Time improvement |
>  | group.any         | Deprecated and replaced by group |
>  | group.all         | Deprecated and replaced by group |
>  | group.min         | Deprecated and replaced by group |
>  | group.max         | Deprecated and replaced by group |
>  | group.min_max     | Deprecated and replaced by group |
>  | group.mean        | Deprecated and replaced by group |
>  | group.med         | Deprecated and replaced by group |
>  | group.mad         | Deprecated and replaced by group |
>  | group.var         | Deprecated and replaced by group |
>  | group.sum         | Deprecated and replaced by group |
>  | groupcolVars      | Deprecated and replaced by "Rfast2colGroup(...,method="var")" |
>  | sort_mat          | Deprecated and replaced by "colSort" and "rowSort" |
>  | AddToNamespace        | remove unused option       |
>  | nth                   | fix bug				     |

</br>
</br>

###  **Version 1.9.3**
***
> <u> **Improved** </u>(_**by speed, correctness or options**_) 
>
>  |	    Function	    |           What's new!             |
>  | ------------------- | ----------------------------------------- |
>  | omp                   | Time improvement.                 |
>  | sort_mat              |  Replaced by "colSort" and "rowSort" and will be removed in the next udate. |
>  | checkUsage            |  Remove uneccessary option. |
>  | read.examples         |  Remove uneccessary option. |
>  | checkTF               |  Remove uneccessary option. |
>  | checkAliases          |  Remove uneccessary option. |
>  | comb_n                |  Add option for return list or matrix. |
>  | rownth                |  Fix of a bug. |
>
> 
> <u>**New**</u>
>
>  | 	     Function	      |          What's new!            |
>  | ------------------- | ----------------------------------------- |
>  | omp	                  |  Multinomial regression now added. |
>  | omp	                  |  Option to standardise the predictor variables. |
>  | cor.fbed	              |  Option to standardise the predictor variables. |
>  | cor.fsreg	              |  Option to standardise the predictor variables. |
>  | el.test2	              |  Empirical likelihood test for two sample means. |

</br>
</br>

###  **Version 1.9.2**
***
> <u> **Improved** </u>(_**by speed, correctness or options**_) 
>
>  |	    Function	    |           What's new!             |
>  | ------------------- | ----------------------------------------- |
>  | omp                   | Time improvement.                 |
>  | ompr                  | Time improvement.                 |
>  | cor.fbed              | Time improvement.                 |
>  | multinom.mle          | Time improvement.                 |
>  | weib.reg              | Fix a bug and time improvement.   |
>  | spatmed.reg           | Time improvement.                 |
>  | invgauss.reg          | Time improvement.                 |
>  | invgauss.regs         | Time improvement.                 |
>  | gammareg              | Time improvement.                 |
>  | gammaregs             | Time improvement.                 |
>  | gammacon              | Time improvement.                 |
>  | colvm.mle             | Time improvement.                 |
>  | gammaregs             | Time improvement.                 |
>  | el.test1              | Time improvement.                 |
>  | Norm                  | Fix a bug.                        |
>  | which_isFactor        |  Removed and replaced by "which.is". |
>  | checkUsage            |  Fix a bug. |
>  | colsums               |  Time improvement for integer matrices. |
>  | rowsums               |  Time improvement for integer matrices. |
>  | med                   |  Faster **version.|**
>  | sort.unique.length    |  slightly faster **version.|**
>  | sort.unique           |  slightly faster **version.|**
>  | Stack                 |  Fix a bug and add function clear for efficient reuse of the Stack. |
>  | read.example          |  Fix a bug. |
>
> 
> <u>**New**</u>
>
>  | 	     Function	      |          What's new!            |
>  | ------------------- | ----------------------------------------- |
>  | omp	                  |  Multinomial regression now added. |
>
> ### **Comments**
>***
> From now on the Rfast needs R **version >= 3.5**
>

</br>
</br>

###  **Version 1.9.1**
***
>
> <u> **Improved** </u>(_**by speed, correctness or options**_) 
>
>  |	      Function	    |     What's new!     |
>  | ------------------- | -------------------- |
>  |     ompr           	| Time improvement.    |
>  |     omp           	| Time improvement .   |
>  |   cholesky            | Time improvement. |
>  |   colnth, rownth      | Add extra options, "na.rm" and "descending" order and "index.return". |
>  |   colrow.zero         | deprecate it and replaced by "colrow.value". |
>
>
> 
> <u>**New**</u>
>
>  | 	   Function		    |                What's new!                |
>  | ----------------	| --------------------------------------- |
>  |     which.is          | The same with which_isFactor but general. Use this instead of which_isFactor. |
>  |     colrow.value      | Search if a column and row is filled with a specific value. |

</br>
</br>

###  **Version 1.9.0**
***
> <u> **Improved** </u>(_**by speed, correctness or options**_) 
>
>  |	      Function	    |     What's new!     |
>  | ----------------	| ----------------- |
>  | topological_sort      |  Time improvement.  |
>  | permcor               |  Fix a bug.         |
>  | omp                   |  Time improvement.  |
>  | normlog.regs          |  Fix a bug.                                    |
>  | cova                  |  Added an extra argument.                      |
>  | ompr                  |  Time improvement.                             |
>  |   as_integer          | Fix a bug.                                  |
>  |   Round               | Fix a bug. digit argument can be up to 15.  |
>  |   matrix.sum          | Deprecate it.                               |
>  |   countNA             | Deprecate it.                               |
>  |   sort_unique.length  | Deprecate it for numeric numbers.           |
>  |   Rank                | Deprecate method="first".                   |
>  |   Match               | Improved.                                   |
>  |   nth                 | Improved method for integers.               |
>  |   colshuffle          | Fix a bug.                                  |
>  |   transpose           | Now can handle generic matrix using parallel|
>  |   colshuffle          | Fix a bug.                                  |
>
> 
> <u>**New**</u>
>
>  | 	   Function		      |                           What's new!                           |
>  | ----------------	  | ------------------------------------------------------------- |
>  | betabinom.mle           |  MLE of the beta binomial distribution.                         |
>  | betageom.mle            |  MLE of the beta geometric distribution.                        |
>  | multivt.mle             |  MLE of the multivariate t distribution.                        |
>  | colpoisson.anovas       |  Column-wise ANOVA with Poisson distribution.                   |
>  | colquasipoisson.anovas  |  Column-wise ANOVA with quasi Poisson.                          |
>  | exact.ttest2            |  Exact permutations 2-sample t-test.                            |
>  | chi2Test                |  Chi-squared test of independence.                              |
>  | gchi2Test               |  G-square and Chi-square tests of indepdence.                   |
>  | chi2tests               |  Many chi-squared tests of independence.                        |
>  | chi2Test_univariate     |  Matrix with chi-square tests of indepedence.                   |
>  | mvlnorm.mle             |  MLE of the multivariate lognormal distribution.                |
>  | poly.cor                |  Polychoric correlation.                                        |
>  | pooled.cov              |  Pooled covariance matrix.                                      |
>  | spatmed.reg             |  Spatial median (multivariate) regression.                      |
>  | sscov                   |  Spatial sign covariance matrix.                                |
>  | trim.colmeans           |  Trimmed column-wise means.                                     |
>  | trim.rowmeans           |  Trimmed row-wise means.                                        |
>  | eigs.sym                |  Extract some principal components from a symmetric matrix.     |
>  | invgauss.regs           |  Many simple inverse Gaussian regressions with a log link.      |
>  | invgauss.reg            |  Inverse Gaussian regression with a log link.                   |
>  | gammaregs               |  Many simple Gamma regressions with a log link.                 |
>  | gammareg                |  Gamma regression with a log link.                              |
>  | gammacon                |  Gamma regression with a constant term only.                    |
>  |   Outer               | like R's outer.                           |
>  |   RemoveFromNamespace | remove exported functions.                |
>  |   Sort.int            | fast sorting integer.                     |
>  |   colCumMaxs          | apply cummax to column.                   |
>  |   colCumSums          | apply cumsum to column.                   |
>  |   colCumMins          | apply cummin to column.                   |
>  |   colCumProds         | apply cumprod to column.                  |
>  |   positive            | apply method to each positive value.      |
>  |   positive.negative   | apply method to each positive and negative value.   |
>  |   negative            | apply method to each negative value.                |
>  |   as.Rfast.function   | convert an R function to Rfast's equivalent         |
>  |   mat.mult            | Generic matrix multiplication using parallel.       |
>  |   checkUsage          | checking usage section in Rd files.                 |
>  |   Hash                | Create Hash object.                                 |
>  |   Hash.key            | Search key or multi key.                            |
>  |   apply.condition     | Apply method to each column using a condition. Only integers.    |
>  |   Stack               | Stack object. See man page.                               |
>  |   iterator            | iterator object. See man page.                            |
>  |   Elem                | access element of an iterator object.                     |
>  |   print.environment   | S3 method for printing environment.                       |
>  |   env.copy            | deep copy environment.                                    |
>  |   ufactor             | Untyped factor object. See man page.                      |

</br>
</br>

###  **Version 1.8.8**
***
> <u> **Improved** </u>(_**by speed, correctness or options**_) 
>
>  |	      Function	    |     What's new!     |
>  | ----------------	| ----------------- |
>  | rowMedians            | Fix a bug  |
>  | mad2                  | Add option handling NAs and fix a bug  |
>  | prop.regs             | Made it more stable.                   |
>  | logistic_only         | Made it more stable.                   |
>  | multinom.regs         | Fix a bug and removed parallel for safety reasons.   |
>  | rint.regs             | Fix a bug.                    |
>  | univglms              | Add the option for quasi Poisson regression     |
>  | hd.eigen              | Now returns eigen-vectors as well.   |
>  | ompr                  | Fix a bug.   |
>  | tobit.mle             | Made it faster.   |
>  | standardise           | Made it faster.   |
>  | auc                   | Made it faster.   |
>  | colaucs               | Made it faster.   |
>  | pc.skel               | Made it faster and added more utilities.  |
>  | cor.fsreg             | Made it faster and fix a bug.  |
>  | allttests             | Made it faster and fix a bug.  |
>  | data.frame.to_matrix  | Add option setting colnames and rownames. Fix a bug   |
>  | nth                   | Add option for handling NAs.  |
>  | Pmax                  | Add option for handling NAs.  |
>  | Pmin                  | Add option for handling NAs.  |
>  | Sort                  | Add option for handling NAs.  |
>  | Table                 | Add option for handling NAs, add option for second argument, deprecate argument as.vector and replaced from "names". |
>  | Round                 | Fix a bug.  |
>  | Norm                  | Fix a bug.  |
>  | colsums               | Add option for sum using specific indices.  |
>  | rowsums               | Add option for sum using specific indices and option for parallel.  |
>  | is_element            | Fix a bug.  |
>  | eachrow               | Deprecate argument suma and replaced from argument method.  |
>  | permutation           | Deprecate argument all and replaced from nperm.  |
>  | permutation.next      | Deprecate argument all.next and replaced from nperm.  |
>  | permutation.prev      | Deprecate argument all.prev and replaced from nperm.  |
>  | data.frame.to_matrix  | Fix a bug and time improvement.  |
>  | Rank                  | Fix a bug for method="first".    |
>  | Match                 | fix bug and time improvement.     |
>  | bc                  | Estimate the optimal lambda in the Box-Cox transformation.  |
>  | ar1                 | AR(1) model.  |
>  | colar1              | Many column-wise AR(1) models.  |
>  | bc                  | Estimate the optimal lambda in the Box-Cox transformation.  |
>  | rbing               | Random values generation from a special case of the Bingham distribution.   |
>  | rbingham            | Random values generation from the Bingham distribution.   |
>  | omp                 | Orthogonal Matching Pursuit allowing many regression models.   |
>  | yule                | Yule's coefficient of colligation.   |
>  | col.yule            | Many column-wise Yule's coefficients of colligation.   |
>  | cox.poisrat         | Test for the ratio of two Poisson means.   |
>  | col.coxpoisrat      | Many column-wise tests for the ratio of two Poisson means.   |
>  | colPmax             | column-wise parallel maxima  |
>  | colPmin             | column-wise parallel minima |
>  | freq.max            | maximum frequency of a number |
>  | freq.min            | minimum frequency of a number |
>  | Pmin_Pmax           | parallel minima-maxima values |
>  | Table.sign          | counting the positive, negative, zeros and NA values. |
>  | topological_sort    | Topological sort of a Directed Acyclic Graph (DAG)    |
>  | countNA             | count the NAs  |
>  | columns             | get specific columns from a matrix  |
>  | rows                | get specific rows from a matrix  |
>  | eachcol.apply       | apply a function to each col after the operation |
>  | checkTF             | checking man files for missing TRUE/FALSE values in examples |

