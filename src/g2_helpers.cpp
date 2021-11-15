// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <cstdlib>
#include "g2t.h"

using namespace std;
using namespace Rcpp;

/////////////// GEORGE ///////////////////////////

List g2Test(IntegerMatrix data, int x, int y, IntegerVector ics, IntegerVector dc,const bool parallel){
    IntegerVector cs = ics - 1;
    --x;
    --y;
    const int ncs = cs.size();
    double statistic = 0;
    int df = 0;
    if (ncs == 0) {
        auto result = g2Test(data, x, y, dc);
        statistic = result.stat;
        df = result.df;
    }else{
        int xdim = dc[x];
        int ydim = dc[y];
        int nsamples = data.nrow();
        vector<int> prod(ncs+1);
        prod[0] = 1;
        for (int i = 1; i <= ncs; ++i) {
            prod[i] = prod[i - 1] * dc[cs[i - 1]];
        }

        int size = prod[ncs];
        vector<vector<int>> counts(size,vector<int>(xdim * ydim));

        if(parallel){
          #pragma omp parallel for
          for (int i = 0; i < nsamples; ++i) {
              int key = 0;
              for (int j = 0; j < ncs; ++j) {
                  key += data(i, cs[j]) * prod[j];
              }
              int curx = data(i, x);
              int cury = data(i, y);
              /*(if (counts[key] == NULL) {
                  counts[key] = new int[xdim * ydim];
                  memset(counts[key], 0, sizeof(int) * xdim * ydim);
              }*/
              #pragma omp critical
              ++counts[key][cury * xdim + curx];
          }

          #pragma omp parallel for reduction(+: statistic)
          for (int i = 0; i < size; ++i) {
              statistic += g2Statistic(counts[i], xdim, ydim);
          }
        }else{
          for (int i = 0; i < nsamples; ++i) {
            int key = 0;
            for (int j = 0; j < ncs; ++j) {
                key += data(i, cs[j]) * prod[j];
            }
            int curx = data(i, x);
            int cury = data(i, y);
            /*(if (counts[key] == NULL) {
                counts[key] = new int[xdim * ydim];
                memset(counts[key], 0, sizeof(int) * xdim * ydim);
            }*/
            ++counts[key][cury * xdim + curx];
          }

          for (int i = 0; i < size; ++i) {
              statistic += g2Statistic(counts[i], xdim, ydim);
          }
        }
        df = (xdim - 1) * (ydim - 1) * prod[ncs];
    }

    return List::create(_["statistic"]=statistic,_["df"]=df);
}


List g2Test_univariate_perm(IntegerMatrix data, IntegerVector dc, int nperm) {
  int nvars = data.ncol();
  IntegerVector trash;
  
  int nout = nvars * (nvars - 1) / 2;
  IntegerVector xout(nout);
  IntegerVector yout(nout);
  NumericVector pvalues(nout);
  NumericVector statistics(nout);
  
  int idx = 0;
  for(int i = 0; i < nvars; ++i) {
    for(int j = i+1; j < nvars; ++j) {
      TestResult result = permG2Test(data, i, j, trash, dc, nperm); // pass trash instead of null
      xout[idx] = i + 1;
      yout[idx] = j + 1;
      pvalues[idx] = result.pvalue;
      statistics[idx] = result.stat;
      ++idx;
    }
  }
  
  List out;
  out["statistic"] = statistics;
  out["pvalue"] = pvalues;
  out["x"] = xout;
  out["y"] = yout;
  return out;
}


List g2tests_perm(IntegerMatrix data, IntegerVector x, int y, IntegerVector dc, int nperm) {
  IntegerVector trash;
	int nout = x.size();
	IntegerVector xout(nout);
	IntegerVector yout(nout);
	NumericVector pvalues(nout);
	NumericVector statistics(nout);

	y = y - 1;
	int xlen = x.size();
	for(int i = 0; i < xlen; ++i) {
		int curx = x[i] - 1;
		TestResult result = permG2Test(data, curx, y, trash, dc, nperm);// pass trash instead of null
		xout[i] = curx;
		yout[i] = y;
		pvalues[i] = result.pvalue;
		statistics[i] = result.stat;
	}

	List out;
	out["statistic"] = statistics;
	out["pvalue"] = pvalues;
	out["x"] = xout;
	out["y"] = yout;
	return out;
}


List g2Test_univariate(IntegerMatrix data, IntegerVector dc) {
    int nvars = data.ncol();

    int nout = nvars * (nvars - 1) / 2;
    IntegerVector xout(nout);
    IntegerVector yout(nout);
    NumericVector statistics(nout);
    IntegerVector df(nout);
  
    int idx = 0;
    for(int i = 0; i < nvars; ++i) {
        for(int j = i+1; j < nvars; ++j) {
            TestResult result = g2Test(data, i, j, dc);
            xout[idx] = i + 1;
            yout[idx] = j + 1;
            statistics[idx] = result.stat;
            df[idx] = (dc[i] - 1) * (dc[j] - 1);
            ++idx;
        }
    }
  
    List out;
    out["statistic"] = statistics;
    out["x"] = xout;
    out["y"] = yout;
    out["df"] = df;
    return out;
}



List g2tests(IntegerMatrix data, IntegerVector x, int y, IntegerVector dc) {
	int nout = x.size();
	int xlen = x.size();
	IntegerVector xout(nout);
	IntegerVector yout(nout);
	NumericVector statistics(nout);
	IntegerVector df(nout);

	y = y-1;
  const int dcy = dc[y];
	for(int i = 0; i < xlen; ++i) {
		int curx =  x[i] - 1;
		TestResult result = g2Test(data, curx, y, dc);
		xout[i] = curx;
		yout[i] = y;
		statistics[i] = result.stat;
		df[i] = (dc[curx] - 1) * (dcy - 1);
	}

	List out;
	out["statistic"] = statistics;
	out["x"] = xout;
	out["y"] = yout;
	out["df"] = df;
	return out;
}


List g2Test_perm(IntegerMatrix data, int x, int y, IntegerVector cs, IntegerVector dc, int nperm) {
  IntegerVector ics = cs - 1;
  TestResult result = permG2Test(data, x-1, y-1, ics, dc, nperm);
  
  List out;
  out["statistic"] = result.stat;
  out["pvalue"] = result.pvalue;
  out["x"] = x;
  out["y"] = y;
  out["df"] = result.df;
  return out;
}

// CHI2 TESTS

List chi2Test(IntegerMatrix data, int x, int y, IntegerVector cs, IntegerVector dc){
    IntegerVector ics = cs - 1;
    --x;
    --y;
    double statistic = 0;
    int df = 0;
    const int ncs = cs.size();
    if (ncs == 0) {
        auto result = chi2Test(data, x, y, dc);
        statistic = result.stat;
        df = result.df;
    }else{
        int xdim = dc[x];
        int ydim = dc[y];
        int nsamples = data.nrow();
        std::vector<int> prod(ncs + 1);
        prod[0] = 1;
        for (int i = 1; i <= ncs; ++i) {
            prod[i] = prod[i - 1] * dc[cs[i - 1]];
        }

        int size = prod[ncs];
        std::vector<std::vector<int>> counts(size,std::vector<int>(xdim*ydim));

        for (int i = 0; i < nsamples; ++i) {
            int key = 0;
            for (int j = 0; j < ncs; ++j) {
                key += data(i, cs[j]) * prod[j];
            }
            int curx = data(i, x);
            int cury = data(i, y);
            if (!counts[key].size()) {
                counts[key].resize(xdim * ydim);
            }
            ++counts[key][cury * xdim + curx];
        }

        for (int i = 0; i < size; ++i) {
            statistic += chi2Statistic(counts[i], xdim, ydim);
        }
        df = (xdim - 1) * (ydim - 1) * prod[ncs];
    }

    return List::create(_["statistic"]=statistic,_["df"]=df);
}


List chi2Test_univariate(IntegerMatrix data, IntegerVector dc) {
  int nvars = data.ncol();
  
  int nout = nvars * (nvars - 1) / 2;
  IntegerVector xout(nout);
  IntegerVector yout(nout);
  NumericVector statistics(nout);
  NumericVector df(nout);
  
  int idx = 0;
  for(int i = 0; i < nvars; ++i) {
    for(int j = i+1; j < nvars; ++j) {
      TestResult result = chi2Test(data, i, j, dc);
        xout[idx] = i + 1;
        yout[idx] = j + 1;
        statistics[idx] = result.stat;
        df[idx] = (dc[i] - 1) * (dc[j] - 1);
        ++idx;
    }
  }
  
  List out;
  out["statistic"] = statistics;
  out["x"] = xout;
  out["y"] = yout;
  out["df"] = df;
  return out;
}


List chi2tests(IntegerMatrix data, IntegerVector x, int y, IntegerVector dc) {

	int nout = x.size();
	IntegerVector xout(nout);
	IntegerVector yout(nout);
	NumericVector statistics(nout);
	NumericVector df(nout);

	y = y-1;
	int xlen = x.size();
	for(int i = 0; i < xlen; ++i) {
		int curx =  x[i] - 1;
		TestResult result = chi2Test(data, curx, y, dc);
		xout[i] = curx;
		yout[i] = y;
		statistics[i] = result.stat;
		df[i] = (dc[curx] - 1) * (dc[y] - 1);
	}

	List out;
	out["statistic"] = statistics;
	out["x"] = xout;
	out["y"] = yout;
	out["df"] = df;
	return out;
}

////////////////// GEORGE //////////////////////////////////////////////////