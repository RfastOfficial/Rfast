#include <iostream>
#include <cstring>
#include <random>
#include <cmath>
#include "g2t.h"

TestResult::TestResult(double _pvalue, double _stat, double _logpvalue, int _df){
	pvalue=_pvalue;
	stat=_stat;
	logpvalue=_logpvalue;
	df=_df;
}

static double g2Statistic(int* counts, int xdim, int ydim) {
  if (counts == NULL) {
    return 0;
  }
  double statistic = 0;
  int countsXY = 0;
  int* countsX = new int[xdim];
  int* countsY = new int[ydim];
  
  memset(countsX, 0, xdim * sizeof(int));
  memset(countsY, 0, ydim * sizeof(int));
  
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

static double chi2Statistic(int* counts, int xdim, int ydim) {
  if (counts == NULL) {
    return 0;
  }
  double statistic = 0;
  int countsXY = 0;
  int* countsX = new int[xdim];
  int* countsY = new int[ydim];
  
  memset(countsX, 0, xdim * sizeof(int));
  memset(countsY, 0, ydim * sizeof(int));
  
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
		  double expected = ((double)(countsX[x] * countsY[y])) / countsXY;
		  statistic += ((curcounts - expected) * (curcounts - expected)) / expected;
        }
      }
    }
  }
  
  delete[] countsX;
  delete[] countsY;
  return statistic;
}

static int totalCounts(int* counts, int xdim, int ydim) {
  if (counts == NULL) {
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
  if (counts == NULL) {
    return;
  }
  memset(countsX, 0, xdim * sizeof(int));
  for (int x = 0; x < xdim; ++x) {
    for (int y = 0; y < ydim; ++y) {
      countsX[x] += counts[y * xdim + x];
    }
  }
}

static void colCounts(int* counts, int xdim, int ydim, int* countsY) {
  if (counts == NULL) {
    return;
  }
  memset(countsY, 0, ydim * sizeof(int));
  for (int x = 0; x < xdim; ++x) {
    for (int y = 0; y < ydim; ++y) {
      countsY[y] += counts[y * xdim + x];
    }
  }
}

TestResult g2Test(NumericMatrix& data, int x, int y, int* dc) {
  int xdim = dc[x];
  int ydim = dc[y];
  int* counts = new int[xdim * ydim];
  memset(counts, 0, sizeof(int) * xdim * ydim);
  
  for (int i = 0; i < data.nrow(); ++i) {
    int curx = (int)data(i, x);
    int cury = (int)data(i, y);
    counts[cury * xdim + curx]++;
  }
  int df = (xdim - 1) * (ydim - 1);
  double statistic = g2Statistic(counts, xdim, ydim);
  
  delete[] counts;
  return TestResult(0, statistic, 0, df);
}

TestResult chi2Test(NumericMatrix& data, int x, int y, int* dc) {
  int xdim = dc[x];
  int ydim = dc[y];
  int* counts = new int[xdim * ydim];
  memset(counts, 0, sizeof(int) * xdim * ydim);
  
  for (int i = 0; i < data.nrow(); ++i) {
    int curx = (int)data(i, x);
    int cury = (int)data(i, y);
    counts[cury * xdim + curx]++;
  }
  int df = (xdim - 1) * (ydim - 1);
  double statistic = chi2Statistic(counts, xdim, ydim);
  
  delete[] counts;
  return TestResult(0, statistic, 0, df);
}

TestResult g2Test(NumericMatrix& data, int x, int y, int* cs, int ncs, int* dc) {
  if (ncs == 0) {
    return g2Test(data, x, y, dc);
  }
  int xdim = dc[x];
  int ydim = dc[y];
  int nsamples = data.nrow();
  int* prod = new int[ncs + 1];
  prod[0] = 1;
  for (int i = 1; i <= ncs; ++i) {
    prod[i] = prod[i - 1] * dc[cs[i - 1]];
  }
  
  int size = prod[ncs];
  int **counts = new int*[size];
  for (int i = 0; i < size; ++i) {
    counts[i] = new int[xdim * ydim];
    memset(counts[i], 0, sizeof(int) * xdim * ydim);
  }
  
  for (int i = 0; i < nsamples; ++i) {
    int key = 0;
    for (int j = 0; j < ncs; ++j) {
      key += (int)data(i, cs[j]) * prod[j];
    }
    int curx = (int)data(i, x);
    int cury = (int)data(i, y);
    if (counts[key] == NULL) {
      counts[key] = new int[xdim * ydim];
      memset(counts[key], 0, sizeof(int) * xdim * ydim);
    }
    counts[key][cury * xdim + curx]++;
  }
  
  double statistic = 0;
  for (int i = 0; i < size; ++i) {
    statistic += g2Statistic(counts[i], xdim, ydim);
  }
  int df = (xdim - 1) * (ydim - 1) * prod[ncs];
  
  delete[] prod;
  for (int i = 0; i < size; ++i) {
    if (counts[i] != NULL)
      delete[] counts[i];
  }
  delete[] counts;
  
  return TestResult(0, statistic, 0, df);
}

TestResult chi2Test(NumericMatrix& data, int x, int y, int* cs, int ncs, int* dc) {
  if (ncs == 0) {
    return chi2Test(data, x, y, dc);
  }
  int xdim = dc[x];
  int ydim = dc[y];
  int nsamples = data.nrow();
  int* prod = new int[ncs + 1];
  prod[0] = 1;
  for (int i = 1; i <= ncs; ++i) {
    prod[i] = prod[i - 1] * dc[cs[i - 1]];
  }
  
  int size = prod[ncs];
  int **counts = new int*[size];
  for (int i = 0; i < size; ++i) {
    counts[i] = new int[xdim * ydim];
    memset(counts[i], 0, sizeof(int) * xdim * ydim);
  }
  
  for (int i = 0; i < nsamples; ++i) {
    int key = 0;
    for (int j = 0; j < ncs; ++j) {
      key += (int)data(i, cs[j]) * prod[j];
    }
    int curx = (int)data(i, x);
    int cury = (int)data(i, y);
    if (counts[key] == NULL) {
      counts[key] = new int[xdim * ydim];
      memset(counts[key], 0, sizeof(int) * xdim * ydim);
    }
    counts[key][cury * xdim + curx]++;
  }
  
  double statistic = 0;
  for (int i = 0; i < size; ++i) {
    statistic += chi2Statistic(counts[i], xdim, ydim);
  }
  int df = (xdim - 1) * (ydim - 1) * prod[ncs];
  
  delete[] prod;
  for (int i = 0; i < size; ++i) {
    if (counts[i] != NULL)
      delete[] counts[i];
  }
  delete[] counts;
  
  return TestResult(0, statistic, 0, df);
}

TestResult permG2Test(NumericMatrix& data, int x, int y, int* cs, int ncs, int* dc, int nperm) {
  int xdim = dc[x];
  int ydim = dc[y];
  
  int nsamples = data.nrow();
  int* prod = new int[ncs + 1];
  prod[0] = 1;
  for (int i = 1; i <= ncs; ++i) {
    prod[i] = prod[i - 1] * dc[cs[i - 1]];
  }
  
  int size = prod[ncs];
  int **counts = new int*[size];
  memset(counts, 0, size * sizeof(int*));
  
  for (int i = 0; i < nsamples; ++i) {
    int key = 0;
    for (int j = 0; j < ncs; ++j) {
      key += (int)data(i, cs[j]) * prod[j];
    }
    int curx = (int)data(i, x);
    int cury = (int)data(i, y);
    if (counts[key] == NULL) {
      counts[key] = new int[xdim * ydim];
      memset(counts[key], 0, sizeof(int) * xdim * ydim);
    }
    counts[key][cury * xdim + curx]++;
  }
  
  double statistic = 0;
  for (int i = 0; i < size; ++i) {
    statistic += g2Statistic(counts[i], dc[x], dc[y]);
  }
  int df = (dc[x] - 1) * (dc[y] - 1) * prod[ncs];
  delete[] prod;
  
  if (nperm == 0) {
    for (int i = 0; i < size; ++i) {
      if (counts[i] != NULL)
        delete[] counts[i];
    }
    delete[] counts;
    return TestResult(0, statistic, 0, df);
  }
  
  double* permstats = new double[nperm];
  memset(permstats, 0, nperm * sizeof(double));
  
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
    if (counts[i] != NULL)
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
        memset(&matrix[l * ncol + m], 0, (ncol - m) * sizeof(int));
        break;
      }
      
      //  Compute the conditional expected value of MATRIX(L,M).
      bool done = false;
      int curnlm, nlm;
      curnlm = nlm = (int)(((double)ia * id) / ie + 0.5);
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
