#include "Rfast.h"


Rfast::internal::NA_helper<string>::type Rfast::internal::NA_helper<string>::val = NA_STRING;
Rfast::internal::NA_helper<int>::type    Rfast::internal::NA_helper<int>::val    = NA_INTEGER;
Rfast::internal::NA_helper<bool>::type   Rfast::internal::NA_helper<bool>::val   = NA_LOGICAL;
Rfast::internal::NA_helper<double>::type Rfast::internal::NA_helper<double>::val = NA_REAL;