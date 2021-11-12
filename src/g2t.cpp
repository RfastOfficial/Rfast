#include <iostream>
#include <cstring>
#include <random>
#include <cmath>
#include "g2t.h"



double g2Statistic(int* counts, int xdim, int ydim) {
  if (counts == nullptr) {
    return 0;
  }
  double statistic = 0;
  int countsXY = 0;
  int* countsX = new int[xdim];
  int* countsY = new int[ydim];
  
  
  for (int x = 0; x < xdim; ++x) {
    for (int y = 0; y < ydim; ++y) {
      int curcounts = counts[y * xdim + x];
      countsXY += curcounts;
      countsX[x] += curcounts;
      countsY[y] += curcounts;
    }
  }
  
  for (int x = 0; x < xdim; ++x) {
    if (countsX[x] != 0) {
      for (int y = 0; y < ydim; ++y) {
        int curcounts = counts[y * xdim + x];
        if (countsY[y] != 0 && curcounts != 0) {
          statistic += curcounts * (log(((double)curcounts * countsXY) / ((double)countsX[x] * countsY[y])));
        }
      }
    }
  }
  
  delete[] countsX;
  delete[] countsY;
  return 2 * statistic;
}

double chi2Statistic(int* counts, int xdim, int ydim) {
  if (counts == nullptr) {
    return 0;
  }
  double statistic = 0;
  int countsXY = 0;
  int* countsX = new int[xdim];
  int* countsY = new int[ydim];
  
  
  for (int x = 0; x < xdim; ++x) {
    for (int y = 0; y < ydim; ++y) {
      int curcounts = counts[y * xdim + x];
      countsXY += curcounts;
      countsX[x] += curcounts;
      countsY[y] += curcounts;
    }
  }
  
  for (int x = 0; x < xdim; ++x) {
    if (countsX[x] != 0) {
      for (int y = 0; y < ydim; ++y) {
        int curcounts = counts[y * xdim + x];
		    double expected = ((double)(countsX[x] * countsY[y])) / countsXY;
		      statistic += ((curcounts - expected) * (curcounts - expected)) / expected;
      }
    }
  }
  
  delete[] countsX;
  delete[] countsY;
  return statistic;
}

static int totalCounts(int* counts, int xdim, int ydim) {
  if (counts == nullptr) {
    return 0;
  }
  int countsXY = 0;
  
  for (int x = 0; x < xdim; ++x) {
    for (int y = 0; y < ydim; ++y) {
      countsXY += counts[y * xdim + x];
    }
  }
  
  return countsXY;
}

static void rowCounts(int* counts, int xdim, int ydim, int* countsX) {
  if (counts == nullptr) {
    return;
  }
  for (int x = 0; x < xdim; ++x) {
    for (int y = 0; y < ydim; ++y) {
      countsX[x] += counts[y * xdim + x];
    }
  }
}

static void colCounts(int* counts, int xdim, int ydim, int* countsY) {
  if (counts == nullptr) {
    return;
  }
  for (int x = 0; x < xdim; ++x) {
    for (int y = 0; y < ydim; ++y) {
      countsY[y] += counts[y * xdim + x];
    }
  }
}

  
TestResult g2Test(IntegerMatrix& data, int x, int y, IntegerVector& dc) {
  int xdim = dc[x];
  int ydim = dc[y];
  int* counts = new int[xdim * ydim];
  
  for (int i = 0; i < data.nrow(); ++i) {
    int curx = data(i, x);
    int cury = data(i, y);
    ++counts[cury * xdim + curx];
  }
  int df = (xdim - 1) * (ydim - 1);
  double statistic = g2Statistic(counts, xdim, ydim);
  
  delete[] counts;
  return TestResult(0, statistic, 0, df);
}

TestResult chi2Test(IntegerMatrix& data, int x, int y, IntegerVector& dc) {
  int xdim = dc[x];
  int ydim = dc[y];
  int* counts = new int[xdim * ydim];
  
  for (int i = 0; i < data.nrow(); ++i) {
    int curx = data(i, x);
    int cury = data(i, y);
    ++counts[cury * xdim + curx];
  }
  int df = (xdim - 1) * (ydim - 1);
  double statistic = chi2Statistic(counts, xdim, ydim);
  
  delete[] counts;
  return TestResult(0, statistic, 0, df);
}

TestResult permG2Test(IntegerMatrix& data, int x, int y, IntegerVector& cs, IntegerVector& dc, int nperm) {
  int xdim = dc[x];
  int ydim = dc[y];
  const int ncs = cs.size();
  int nsamples = data.nrow();
  int* prod = new int[ncs + 1];
  prod[0] = 1;
  for (int i = 1; i <= ncs; ++i) {
    prod[i] = prod[i - 1] * dc[cs[i - 1]];
  }
  
  int size = prod[ncs];
  int **counts = new int*[size];
  
  for (int i = 0; i < nsamples; ++i) {
    int key = 0;
    for (int j = 0; j < ncs; ++j) {
      key += data(i, cs[j]) * prod[j];
    }
    int curx = data(i, x);
    int cury = data(i, y);
    if (counts[key] == nullptr) {
      counts[key] = new int[xdim * ydim];
    }
    ++counts[key][cury * xdim + curx];
  }
  
  double statistic = 0;
  for (int i = 0; i < size; ++i) {
    statistic += g2Statistic(counts[i], dc[x], dc[y]);
  }
  int df = (dc[x] - 1) * (dc[y] - 1) * prod[ncs];
  delete[] prod;
  
  if (nperm == 0) {
    for (int i = 0; i < size; ++i) {
      if (counts[i] != nullptr)
        delete[] counts[i];
    }
    delete[] counts;
    return TestResult(0, statistic, 0, df);
  }
  
  double* permstats = new double[nperm];
  
  std::random_device rd;
  std::mt19937 rng(rd());
  
  int* jwork = new int[ydim - 1];
  int* ct = new int[xdim * ydim];
  int* rowcounts = new int[xdim];
  int* colcounts = new int[ydim];
  int* totals = new int[size];
  double* nrc = new double[xdim * ydim];
  
  // Pre-computing total counts for each conditioning bin
  int maxtotal = 0;
  for (int i = 0; i < size; ++i) {
    totals[i] = totalCounts(counts[i], xdim, ydim);
    maxtotal = (totals[i] > maxtotal) ? totals[i] : maxtotal;
  }
  
  // Pre-computing logarithms and log factorials
  double* plog = new double[1 + maxtotal];
  double* logfact = new double[1 + maxtotal];
  plog[0] = 0;
  logfact[0] = 0;
  for (int j = 1; j <= maxtotal; ++j) {
    plog[j] = log((double)j);
    logfact[j] = logfact[j - 1] + plog[j];
  }
  
  for (int i = 0; i < size; ++i) {
    int ntotal = totals[i];
    if (ntotal > 0) {
      rowCounts(counts[i], xdim, ydim, rowcounts);
      colCounts(counts[i], xdim, ydim, colcounts);
      int ctr = 0;
      for (int x = 0; x < xdim; ++x) {
        for (int y = 0; y < ydim; ++y) {
          nrc[ctr++] = plog[ntotal] - plog[rowcounts[x]] - plog[colcounts[y]];
        }
      }
      
      for (int p = 0; p < nperm; ++p) {
        memcpy(jwork, colcounts, (ydim - 1) * sizeof(int));
        randomContigencyTable(ct, rowcounts, colcounts, xdim, ydim, logfact, jwork, ntotal, rng);
        
        double curstat = 0;
        int cti = 0;
        for (int x = 0; x < xdim; ++x) {
          if (rowcounts[x] != 0) {
            for (int y = 0; y < ydim; ++y) {
              // Indexing that way to avoid cache misses
              curstat += ct[cti] * (plog[ct[cti]] + nrc[cti]);
              ++cti;
            }
          }
          else {
            cti += ydim;
          }
        }
        permstats[p] += (2 * curstat);
      }
    }
  }
  delete[] jwork;
  delete[] ct;
  delete[] rowcounts;
  delete[] colcounts;
  delete[] totals;
  delete[] nrc;
  delete[] logfact;
  delete[] plog;
  for (int i = 0; i < size; ++i) {
    if (counts[i] != nullptr)
      delete[] counts[i];
  }
  delete[] counts;
  
  double pvalue = 1;
  for (int p = 0; p < nperm; ++p) {
    pvalue += (permstats[p] >= statistic);
  }
  pvalue /= (nperm + 1);
  delete[] permstats;
  
  return TestResult(pvalue, statistic, log(pvalue), df);
}

void randomContigencyTable(int* matrix, const int* nrowt, const int* ncolt, const int nrow, const int ncol, const double *logfact, int *jwork, const int ntotal, std::mt19937 &rng) {
  std::uniform_real_distribution<> dist(0, 1);
  int jc = ntotal;
  int ib = 0;
  
  for (int l = 0; l < nrow - 1; ++l) {
    int ia = nrowt[l];
    int ic = jc;
    jc -= ia;
    
    for (int m = 0; m < ncol - 1; ++m) {
      int id = jwork[m];
      int ie = ic;
      ib = ie - ia;
      int ii = ib - id;
      ic -= id;
      
      if (ie == 0) {
        ia = 0;
        break;
      }
      
      //  Compute the conditional expected value of MATRIX(L,M).
      bool done = false;
      int curnlm, nlm;
      curnlm = nlm = (((double)ia * id) / ie + 0.5);
      double x = exp(logfact[ia] + logfact[ib] + logfact[ic] + logfact[id] - logfact[ie] - logfact[nlm] - logfact[id - nlm] - logfact[ia - nlm] - logfact[ii + nlm]);
      for (double r = dist(rng), sumprb = x; !done && r > x; r = sumprb * dist(rng)) {
        bool lsp = false;
        int nll;
        double curx, y;
        curnlm = nll = nlm;
        curx = y = sumprb = x;
        
        //  Increment entry in row L, column M.
        do {
          int j = (id - curnlm) * (ia - curnlm);
          if (j == 0) {
            lsp = true;
            for (j = nll * (ii + nll); j != 0; j = nll * (ii + nll)) {
              --nll;
              y = (y * j) / ((id - nll) * (ia - nll));
              sumprb = sumprb + y;
              
              if (r <= sumprb) {
                curnlm = nll;
                done = true;
                break;
              }
            }
          }
          else {
            ++curnlm;
            curx = (curx * j) / (curnlm * (ii + curnlm));
            sumprb = sumprb + curx;
            
            if (r <= sumprb) {
              done = true;
              break;
            }
            for (j = nll * (ii + nll); j != 0; j = nll * (ii + nll)) {
              --nll;
              y = (y * j) / ((id - nll) * (ia - nll));
              sumprb = sumprb + y;
              
              if (r <= sumprb) {
                curnlm = nll;
                done = true;
                break;
              }
              if (!lsp) {
                break;
              }
            }
          }
        } while (!done && !lsp);
      }
      matrix[l * ncol + m] = curnlm;
      ia = ia - curnlm;
      jwork[m] = jwork[m] - curnlm;
    }
    matrix[l * ncol + ncol - 1] = ia;
  }
  memcpy(&matrix[(nrow - 1) * ncol], jwork, (ncol - 1) * sizeof(int));
  matrix[(nrow - 1) * ncol + ncol - 1] = ib - matrix[(nrow - 1) * ncol + ncol - 2];
}
