/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

 /*
  * File:   mxm.h
  * Author: borbudak
  *
  * Created on May 1, 2016, 6:52 PM
  */


#ifndef RFAST
#define RFAST

#include <Rcpp.h>
#include <random>
using namespace Rcpp;

class TestResult {
public:
	double pvalue;
	double logpvalue;
	double stat;
	int df;

	TestResult::TestResult(double _pvalue, double _stat, double _logpvalue, int _df) 
	: pvalue(_pvalue), stat(_stat), logpvalue(_logpvalue), df(_df) {}
};

void randomContigencyTable(int* matrix, const int* nrowt, const int* ncolt, const int nrow, const int ncol, const double *logfact, int *jwork, const int ntotal, std::mt19937 &rng);
TestResult g2Test(IntegerMatrix& data, int x, int y, IntegerVector& dc);
TestResult g2Test(IntegerMatrix& data, int x, int y, IntegerVector& cs, IntegerVector& dc);
TestResult chi2Test(NumericMatrix& data, int x, int y, IntegerVector& dc);
TestResult chi2Test(NumericMatrix& data, int x, int y, IntegerVector& cs, IntegerVector& dc);
TestResult permG2Test(NumericMatrix& data, int x, int y, IntegerVector& cs, IntegerVector& dc, int nperm);

#endif
