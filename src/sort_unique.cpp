//Author: Manos Papadakis

#include <RcppArmadillo.h>
#include <vector>
#include "Rfast.h"

using namespace Rcpp;

using std::vector;


vector<int> sort_unique_int(vector<int> x){
    int aa,mx,mn,count_not_zero=0;
    bool has_pos=false,has_neg=false;
    max_neg_pos<int>(&x[0],&x[x.size()-1]+1,mx,mn,has_pos,has_neg);
    vector<int> pos,f,neg;
    vector<int>::iterator a=x.begin(),F;
    const int init_pos=-1,init_neg=0;
    if(has_pos){
        pos.resize(mx+1,init_pos);
    }
    if(has_neg){
        neg.resize(1-mn,init_neg);
    }
    if(has_pos && has_neg){
        for(;a!=x.end();++a){
            aa=*a;
            if(aa<0 && neg[-aa]==init_neg){
                neg[-aa]=aa;
                ++count_not_zero;
            }else if(pos[aa]==init_pos){
                pos[aa]=aa;
                ++count_not_zero;
            }
        }
    }else if(has_pos){
        for(;a!=x.end();++a){
            aa=*a;
            if(pos[aa]==init_pos){
                pos[aa]=aa;
                ++count_not_zero;
            }
        }
    }else{ 
        for(;a!=x.end();++a){
            aa=*a;
            if(neg[-aa]==init_neg){
                neg[-aa]=aa;
                ++count_not_zero;
            }
        }
    }
    f.resize(count_not_zero);
    F=f.begin();
    if(has_neg){
        for(auto nr=neg.rbegin();nr!=neg.rend();++nr){
            if(*nr!=init_neg){
                *F++=*nr;
            }
        }
    }
    if(has_pos){
        for(a=pos.begin();a!=pos.end();++a){
            if(*a!=init_pos){
                *F++=*a;
            }
        }
    }
    return f;
}





/*//[[Rcpp::export]]
vector<int> sort_unique_int2(vector<int> x){
  int aa,mx,mn,count_not_zero=0;
  int has_pos=0;
  max_neg_pos<int>(&x[0],&x[x.size()-1]+1,mx,mn,has_pos);
  const int has_neg=x.size()-has_pos;
  vector<int> pos,f,neg;
  vector<int>::iterator a=x.begin(),F;
  const int init_pos=-1,init_neg=0;
  if(has_pos>0){
    pos.resize(mx+1,init_pos);
  }
  if(has_neg>0){
    neg.resize(1-mn,init_neg);
  }
  if(has_pos && has_neg){
    for(;a!=x.end();++a){
      aa=*a;
      aa<0 ? neg[-aa]=aa : pos[aa]=aa;
    }
  }else if(has_pos){
    for(;a!=x.end();++a){
      aa=*a;
      pos[aa]=aa;
    }
    
  }else{ 
    for(;a!=x.end();++a){
      aa=*a;
      neg[-aa]=aa;
    }
  }
  if(has_neg){
    count_not_zero=neg.size()-count(neg.begin(),neg.end(),init_neg);
  }
  if(has_pos){
    count_not_zero+=pos.size()-count(pos.begin(),pos.end(),init_pos);
  }
  f.resize(count_not_zero);
  std::copy_if(pos.begin(),pos.end(),f.begin(),[&](int& v){ return v!=init_pos;});
  return f;
}*/

RcppExport SEXP Rfast_sort_unique_int(SEXP xSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< vector<int> >::type x(xSEXP);
    __result = sort_unique_int(x);
    return __result;
END_RCPP
}


////////////////////////////////////////////////////////////////////////


//[[Rcpp::export]]
vector<double> sort_unique_double(vector<double> x){
  sort(x.begin(),x.end());
  x.erase( unique( x.begin(), x.end() ), x.end() );
  return x;
}

RcppExport SEXP Rfast_sort_unique_double(SEXP xSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< vector<double> >::type x(xSEXP);
    __result = sort_unique_double(x);
    return __result;
END_RCPP
}


////////////////////////////////////////////////////////////////////////




using std::vector;


//[[Rcpp::export]]
int len_sort_unique_int(IntegerVector x){
    int aa,mx,mn,count_not_zero=0;
    bool has_neg=false,has_pos=false;
    max_neg_pos<int>(&x[0],&x[x.size()-1]+1,mx,mn,has_pos,has_neg);
    vector<int> pos,f,neg;
    vector<int>::iterator pp,nn;
    IntegerVector::iterator a=x.begin();
    if(has_pos)
        pos.resize(mx+1,INT_MAX);
    if(has_neg)
        neg.resize(1-mn,INT_MAX);
    if(has_pos && has_neg){
        for(nn=neg.begin(),pp=pos.begin();a!=x.end();++a){
            aa=*a;
            if(aa<0 && nn[-aa]==INT_MAX){
                nn[-aa]=aa;
                ++count_not_zero;
            }else if(pp[aa]==INT_MAX){
                pp[aa]=aa;
                ++count_not_zero;
            }
        }
    }else if(has_pos){
        for(pp=pos.begin();a!=x.end();++a){
            aa=*a;
            if(pp[aa]==INT_MAX){
                pp[aa]=aa;
                ++count_not_zero;
            }
        }
    }else{ 
        for(nn=neg.begin();a!=x.end();++a){
            aa=*a;
            if(nn[-aa]==INT_MAX){
                nn[-aa]=aa;
                ++count_not_zero;
            }
        }
    }
    return count_not_zero;
}

RcppExport SEXP Rfast_len_sort_unique_int(SEXP xSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< IntegerVector >::type x(xSEXP);
    __result = len_sort_unique_int(x);
    return __result;
END_RCPP
}


//////////////////////////////////////////////////////////////////
