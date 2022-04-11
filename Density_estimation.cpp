#include <Rcpp.h>
// #include <dqrng.h>
// #include <cmath>
using namespace Rcpp;
// [[Rcpp::plugins("cpp11")]]

double minimum(NumericVector x) {
  double m = x[0];
  for (auto v : x) {
    if (v < m)
      m = v;
  }
  return m;
}

double maximum(NumericVector x) {
  double m = x[0];
  for (auto v : x) {
    if (v > m)
      m = v;
  }
  return m;
}

double C_mean(NumericVector x) {
  double sum = 0;
  for(auto v : x)
    sum += v;
  return sum / x.size();
}

double C_sd(NumericVector x) {
  double mean = C_mean(x);
  double sum = 0;
  for(auto v : x)
    sum += pow(v - mean, 2);
  return sqrt(sum / (double)(x.size() - 1));
}

double C_which_max(NumericVector v) {
  int n = v.size();
  double tmp = v[0];
  int m = 0;
  for (int i = 1; i < n; ++i) {
    if (v[i] > tmp) {
      m = i;
      tmp= v[i];
    }
  }
  return m;
}

// [[Rcpp::export]]
double C_epanechnikov(double x) {
  if(abs(x) > 1) {
    return 0;
  } else {
    return 0.75 * (1 - pow(x, 2));
  }
}

double f_i(NumericVector x, double i, double bw) {
  int n = x.size();
  double f = 0;
  for(auto v : x) {
    f += C_epanechnikov((x[i] - v) / bw);
  }
  f -= C_epanechnikov(0);
  f /= bw * (n - 1);
  return f;
}

// [[Rcpp::export]]
NumericVector C_seq(double start, double end, double inc) {
  int n = floor((end - start) / inc) + 1;
  NumericVector ret(n);
  for(int i = 0; i < n; i++) {
    ret[i] = start + inc * i;
  }
  return ret;
}

std::vector<double> C_my_breaks(std::vector<double> x, int h = 5) {
  int n = x.size();
  std::sort(x.begin(), x.end());
  if(h >= n)
    return x;
  
  int n_breaks = n / h;
  
  std::vector<double> breaks(n_breaks + 1);
  breaks[0] = x[0];
  
  int k = 0;
  while(breaks[0] == x[k + 1]) {
    k += 1;
  }
  
  if(k == 0) {
    breaks[1] = x[k + h - 1];
    k += h - 1;
  } else {
    breaks[1] = x[k + h];
    k += h;
  }
  
  while(breaks[1] == x[k + 1]) {
    k += 1;
  }
  
  int j = 0;
  for(int i = 2; i < (n_breaks + 2); i++) {
    if(k >= n - h) {
      breaks[i-1] = x[n-1];
      j = i - 1;
      break;
    }
    
    breaks[i] = x[k + h];
    if(k + h == (n-1)) {
      j = i;
      break;
    }
    k += h;
    while((breaks[i] == x[k + 1]) && ((k+h) < n)) {
      k += 1;
    }
  }
  for(int i = (j+1); i < (n_breaks+1); i++){
    breaks.pop_back();
  }
  return breaks;
}

NumericVector C_my_breaks2(NumericVector x, int h = 5) {
  
  int n = x.size();
  std::sort(x.begin(), x.end());
  
  if(h >= n) {
    NumericVector ret = NumericVector::create(x[0], x[n-1]);
    return ret;
  }
  
  double len = x[n-1] - x[0];
  
  NumericVector breaks(n);
  breaks[0] = x[0];
  
  int k = 0;
  while(breaks[0] == x[k + 1]) {
    k += 1;
  }
  
  if((k + 1) >= h) {
    breaks[1] = x[k + 1];
    k += 1;
  } else {
    breaks[1] = x[h - 1];
    k = h - 1;
  }
  
  while(breaks[1] == x[k + 1]) {
    k += 1;
  }
  
  int j = 0;
  int i = 1;
  while(true) {
    i += 1;
    if(k + h >= n) {
      breaks[i] = x[n-1];
      j = i;
      break;
    }
    
    if(((x[k + h] - breaks[i-1]) / len) > 0.05) {
      double inc = (x[k + h] - breaks[i-1]) / 7;
      for(int s = 0; s < 7; s++)
        breaks[i + s] = breaks[i - 1] + inc * (s + 1);
      i = i + 6;
    } else {
      breaks[i] = x[k + h];
    }
    
    // breaks[i] = x[k + h];
    k += h;
    
    if((k + 1) >= n) {
      j = i;
      break;
    }
    
    while((breaks[i] == x[k + 1]) && ((k+1) < n)) {
      k += 1;
    }
  }
  
  breaks.erase(j+1, n);
  
  return breaks;
}

// [[Rcpp::export]]
NumericVector C_count_my_breaks(NumericVector xx, NumericVector breaks) {
  NumericVector x = clone(xx);
  std::sort(x.begin(), x.end());
  int n_breaks = breaks.size();
  int n_data = x.size();
  NumericVector counts(n_breaks - 1);
  int k = 0;

  for (int i = 1; i < n_breaks; i++) {
    while (x[k] <= breaks[i]) {
      counts[i - 1] += 1;
      if (k == (n_data - 1))
        return counts;
      k += 1;
    }
  }
  return counts;
}

// [[Rcpp::export]]
NumericVector C_my_midpoints(NumericVector b) {
  int n = b.size();
  NumericVector m(n-1);
  for (int i = 0; i < (n - 1); i++) {
    m[i] = (b[i + 1] - b[i]) / 2 + b[i];
  }
  return m;
}

// NumericMatrix C_count_my_breaks(NumericVector x, NumericVector breaks) {
//   std::sort(x.begin(), x.end());
//   int n_breaks = breaks.size();
//   int n_data = x.size();
//   NumericMatrix counts(n_breaks-1, 2);
//   int k = 0;
//   
//   for(int i = 1; i < n_breaks; i++) {
//     for(int j = k; j < n_data; j++) {
//       if(x[j] <= breaks[i]) {
//         counts(i-1, 0) += 1;
//         counts(i-1, 1) += x[j];
//       } else {
//         k = j;
//         break;
//       }
//     }
//     counts(i-1, 0) != 0 ? counts(i-1, 1) /= counts(i-1, 0) : 0;
//   }
//   return counts;
// }

NumericVector C_quartile(NumericVector x) {
  std::sort(x.begin(), x.end());
  int n = x.size();
  NumericVector ret(2);
  ret[0] = x[ceil(n * 0.25)];
  ret[1] = x[ceil(n * 0.75)];
  return ret;
}

// [[Rcpp::export]]
double C_amise(NumericVector x) {
  NumericVector q = C_quartile(x);
  double sig_hat = sd(x);
  double sig_tilde = (q[1] - q[0]) / 1.35;
  double sig = sig_hat < sig_tilde ? sig_hat : sig_tilde;
  
  return pow((8 * sqrt(3.14159265) * 0.6) / (3 * pow(0.2, 2)), 0.2) * sig * pow(x.size(), -0.2);
}

// // [[Rcpp::depends(dqrng)]]
// // [[Rcpp::export]]
// IntegerVector C_cv(NumericVector x, int bins, NumericVector search) {
//   int n = x.size();
//   IntegerVector s = dqrng::dqsample_int(n, n);
//   IntegerVector i(n);
//   int each = n / bins;
//   int v = 1;
//   for (int j = bins - 1; j > 0; j--) {
//     for (int k = 0; k < each; k++) {
//       i[n-v] = j;
//       v += 1;
//     }
//   }
//   NumericVector bw = C_seq(search[0], search[1], search[2]);
//   int n_bw = bw.size()
//   NumericVector l_CV(n_bw);
// 
//   for (int h = 0; h < n_bw; h++) {
//     for (int j : s) {
// 
//     }
//   }
// 
// }

// cv <- function(x, bins = 10) {
//   n <- length(x)
//   s <- sample.int(n, n)
//   i <- rep(1:bins, each = as.integer(n/bins))
//   bw <- seq(0.2, 2.1, 0.05)
//   l_CV <- numeric(length(bw))
//   
//   for (h in seq_along(bw)) {
//     for (j in s) {
//       k <- i[which(s==j)]
//       l_CV[h] <- l_CV[h] + 
//         log(1e-8 + sum(epanechnikov((x[j] - x[s[i!=k]]) / bw[h])) / 
//         (bw[h] * length(which(i!=k))))
//     }
//   }
//   return(bw[which.max(l_CV)])
// }

// [[Rcpp::export]]
double C_loocv(NumericVector x, NumericVector search) {
  int n = x.size();
  NumericVector bw = C_seq(search[0], search[1], search[2]);
  NumericVector l_cv(bw.size());
  for(int h = 0; h <  bw.size(); h++) {
    for(int i = 0; i < n; i++) {
      l_cv[h] += log(0.001 + f_i(x, i, bw[h]));
    }
  }
  return bw[C_which_max(l_cv)];
}

// [[Rcpp::export]]
NumericMatrix C_my_density2(NumericVector x, double bw, int m = 512, bool binning = true) {
  int n = x.size();
  double rg_min = minimum(x);
  double rg_max = maximum(x);
  NumericMatrix ret(m, 2);
  
  double xx_inc = (rg_max - rg_min + 6 * bw) / (double)(m - 1);
  ret(0, 0) = rg_min - 3 * bw;
  for(int i = 1; i < m; i++) {
    ret(i, 0) = ret(i-1, 0) + xx_inc;
  }
  
  bw *= 2.236068;
  
  NumericVector breaks = C_my_breaks2(x, 8);
  NumericVector midpoints = C_my_midpoints(breaks);
  NumericVector counts_means = C_count_my_breaks(x, breaks);
  int n_counts = breaks.size() - 1;
  
  if(binning) {
    for(int i = 0; i < m; i++) {
      for(int j = 0; j < n_counts; j++) {
        ret(i, 1) += counts_means[j] * C_epanechnikov((ret(i, 0) - midpoints[j]) / bw);
      }
      ret(i, 1) /= (bw * n);
    }
  } else {
    for(int i = 0; i < m; i++) {
      for(int j = 0; j < n; j++) {
        ret(i, 1) += C_epanechnikov((ret(i, 0) - x[j]) / bw);
      }
      ret(i, 1) /= (bw * n);
    }
  }
  
  return ret;
}

// [[Rcpp::export]]
List C_my_density(NumericVector x, double bw, int m = 512, bool binning = true) {
  int n = x.size();
  double rg_min = minimum(x);
  double rg_max = maximum(x);
  NumericVector xx(m);
  NumericVector y(m);

  xx = C_seq(rg_min - 3 * bw, rg_max + 3 * bw, (rg_max - rg_min + 6 * bw) / (double)(m - 1));
  
  bw *= 2.236068;
  
  int n_bins = (rg_max - rg_min) * pow(n, 0.3) / (3.5 * pow(0.2, 0.5));

  NumericVector breaks = C_seq(rg_min, rg_max, (rg_max - rg_min) / (n_bins - 1));
  NumericVector midpoints = C_my_midpoints(breaks);
  NumericVector counts_means = C_count_my_breaks(x, breaks);
  int n_counts = breaks.size() - 1;

  if(binning) {
    for(int i = 0; i < m; i++) {
      for(int j = 0; j < n_counts; j++) {
        y[i] += counts_means[j] * C_epanechnikov((xx[i] - midpoints[j]) / bw);
      }
      y[i] /= (bw * n);
    }
  } else {
    for(int i = 0; i < m; i++) {
      for(int j = 0; j < n; j++) {
        y[i] += C_epanechnikov((xx[i] - x[j]) / bw);
      }
      y[i] /= (bw * n);
    }
  }
  List L = List::create(Named("x") = xx, Named("y") = y);
  return L;
}



