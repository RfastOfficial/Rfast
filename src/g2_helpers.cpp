// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <cstdlib>
#include "g2t.h"

using namespace std;
using namespace Rcpp;

/////////////// GEORGE ///////////////////////////

List g2Test(NumericMatrix data, int x, int y, IntegerVector cs, IntegerVector dc){
  IntegerVector ics = cs -1;
  TestResult result = g2Test(data, x-1, y-1, ics, dc);
  return List::create(_["statistic"]=result.stat,_["df"]=result.df);
}


List g2Test_univariate_perm(NumericMatrix data, IntegerVector dc, int nperm) {
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
      TestResult result = permG2Test(data, i, j, trash, 0, idc, nperm); // pass trash instead of null
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


List g2tests_perm(NumericMatrix data, IntegerVector x, int y, IntegerVector dc, int nperm) {
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
		TestResult result = permG2Test(data, curx, y, trash, 0, dc, nperm);// pass trash instead of null
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


List g2Test_univariate(NumericMatrix data, IntegerVector dc) {
    int nvars = data.ncol();

    int nout = nvars * (nvars - 1) / 2;
    IntegerVector xout(nout);
    IntegerVector yout(nout);
    NumericVector statistics(nout);
    NumericVector df(nout);
  
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



List g2tests(NumericMatrix data, IntegerVector x, int y, IntegerVector dc) {
	int nout = x.size();
	IntegerVector xout(nout);
	IntegerVector yout(nout);
	NumericVector statistics(nout);
	NumericVector df(nout);

	y = y-1;
	int xlen = x.size();
	for(int i = 0; i < xlen; ++i) {
		int curx =  x[i] - 1;
		TestResult result = g2Test(data, curx, y, dc);
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


List g2Test_perm(NumericMatrix data, int x, int y, IntegerVector cs, IntegerVector dc, int nperm) {
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

List chi2Test(NumericMatrix data, int x, int y, IntegerVector cs, IntegerVector dc){
  IntegerVector ics = cs - 1;
  TestResult result = chi2Test(data, x-1, y-1, ics, dc);
  return List::create(_["statistic"]=result.stat,_["df"]=result.df);;
}


List chi2Test_univariate(NumericMatrix data, IntegerVector dc) {
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


List chi2tests(NumericMatrix data, IntegerVector x, int y, IntegerVector dc) {

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