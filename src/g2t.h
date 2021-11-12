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
#include <vector>
using namespace Rcpp;
using std::vector;

struct TestResult {
	double pvalue;
	double stat;
	double logpvalue;
	int df;

	TestResult(double _pvalue, double _stat, double _logpvalue, int _df) 
	: pvalue(_pvalue), stat(_stat), logpvalue(_logpvalue), df(_df) {}
};

void randomContigencyTable(int*, const int*, const int*, const int, const int, const double *, int *, const int , std::mt19937&);
double g2Statistic(vector<int>&, int, int);
double chi2Statistic(vector<int>&, int, int);
TestResult g2Test(IntegerMatrix&, int, int, IntegerVector& );
TestResult chi2Test(IntegerMatrix&, int, int, IntegerVector& );
TestResult chi2Test2(IntegerMatrix&, int, int, IntegerVector& , IntegerVector&);
TestResult permG2Test(IntegerMatrix&, int, int, IntegerVector& , IntegerVector&, int);

#endif
