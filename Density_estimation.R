epanechnikov <- function(x) {
  (abs(x) <= 1) * 0.75 * (1 - x^2)
}

amise <- function(x) {
  sig <- min(sd(x), sum(quantile(x, c(0.25, 0.75))*c(-1,1))/1.35)
  return(((8 * sqrt(pi) * 0.6) / (3 * 0.2^2))^0.2 * sig * length(x)^(-0.2))
}

loocv <- function(x) {
  n <- length(x)
  f_i <- function(i, bw) 
    1 / (bw * (n - 1)) * (sum(epanechnikov((x[i] - x) / bw)) - epanechnikov(0))
  bw <- seq(0.2, 2.1, 0.05)
  l_CV <- numeric(length(bw))
  for (h in seq_along(bw)) {
    for (i in seq_along(x)) {
      l_CV[h] <- l_CV[h] + log(0.00001 + f_i(i, bw[h]))
    }
  }
  return(bw[which.max(l_CV)])
}

cv <- function(x, bins = 10) {
  n <- length(x)
  s <- sample.int(n, n)
  i <- rep(bins, n)
  i[1:(as.integer(n/bins)*(bins-1))] <- rep(1:(bins-1), each = as.integer(n/bins))
  bw <- seq(0.2, 2.1, 0.05)
  l_CV <- numeric(length(bw))
  
  for (h in seq_along(bw)) {
    for (j in s) {
      k <- i[which(s==j)]
      l_CV[h] <- l_CV[h] + log(0.001 + sum(epanechnikov((x[j] - x[s[i!=k]]) / bw[h])) / (bw[h] * length(which(i!=k))))
    }
  }

  return(bw[which.max(l_CV)])
}

kernel_gen <- function(func, sigma_K, two_norm_sq) {
  structure(
    list(kernel = func, 
         sigma_K = sigma_K,
         two_norm_sq = two_norm_sq),
    class = "Kernel_class"
  )
}

my_breaks <- function(x, h = 5) {
  rg <- range(x)
  len <- rg[2] - rg[1]
  x <- sort(x)
  n <- length(x)
  if(h >= n)
    return(c(x[1], x[n]))
  
  breaks <- numeric(n)
  breaks[1] <- x[1]
  
  k <- 0
  while(breaks[1] == x[k + 1]) {
    k <- k + 1
  }
  
  if(k >= h) {
    breaks[2] <- x[k + 1]
  } else {
    breaks[2] <- x[h]
  }
  k <- h
  while(breaks[2] == x[k + 1]) {
    k  <- k + 1
  }
  
  i <- 2
  while(TRUE) {
    i <- i + 1
    if((k + h) > n) {
      breaks[i-1] <- x[n]
      i <- i - 1
      break
    }
    
    breaks[i] <- x[k + h]
    k <- k + h
    
    if((k + 1) > n) break
    if(breaks[i] == x[k + 1]) {
      k <- max(which(x == breaks[i]))
    }
    
  }
  breaks[1:i]
}
my_breaks2 <- function(x, h = 5) {
  rg <- range(x)
  len <- rg[2] - rg[1]
  x <- sort(x)
  n <- length(x)
  if(h >= n)
    return(c(x[1], x[n]))
  
  breaks <- numeric(n)
  breaks[1] <- x[1]
  
  k <- 0
  while(breaks[1] == x[k + 1]) {
    k <- k + 1
  }
  
  if(k >= h) {
    breaks[2] <- x[k + 1]
    k <- k + 1
  } else {
    breaks[2] <- x[h]
    k <- h
  }
  
  while(breaks[2] == x[k + 1]) {
    k  <- k + 1
  }
    
  i <- 2
  while(TRUE) {
    i <- i + 1
    if((k + h) > n) {
      breaks[i] <- x[n]
      break
    }
    if(((x[k + h] - breaks[i-1]) / len) > 0.05) {
      breaks[i:(i+6)] <- seq(breaks[i-1], x[k + h], length.out = 8)[2:8]
      i <- i + 6
    } else {
      breaks[i] <- x[k + h]
    }
    
    # breaks[i] <- x[k + h]
    k <- k + h
    
    if((k + 1) > n) break
    if(breaks[i] == x[k + 1]) {
      k <- max(which(x == breaks[i]))
    }
    
  }
  breaks[1:i]
}

count_my_breaks <- function(x, breaks) {
  x <- sort(x)
  n_breaks <- length(breaks)
  n_data <- length(x)
  counts <- numeric(n_breaks-1)
  k <- 1
  for (i in 1:(n_breaks-1)) {
    while (x[k] <= breaks[i + 1]) {
      counts[i] <- counts[i] + 1
      if (k == n_data)
        return(counts)
      k <- k + 1
    }
  }
  return(counts)
}

mean_in_breaks <- function(x, counts) {
  x <- sort(x)
  n <- length(counts)
  means <- numeric(n)
  means[1] <- mean(x[1:counts[1]])
  for (i in seq_along(counts)[-1]) {
    if(counts[i] == 0) next
    means[i] <- mean(x[(sum(counts[1:(i-1)])+1):sum(counts[1:i])])
  }
  means
}

my_density <- function(x, kern = NULL, bw = "amise", m = 512, binning = TRUE, mid = FALSE, equi_bin = TRUE) {
  if(class(kern) == "Kernel_class") {
    kern_func <- kern$kernel
    sigma_K <- kern$sigma_K
    two_norm_sq <- kern$two_norm_sq
  } else if(class(kern) == "function") {
    kern_func <- kern
    sigma_K <- integrate(function(x) x^2 * kern(x), lower = -Inf, upper = Inf)$value
    two_norm_sq <- integrate(function(x) kern(x)^2, lower = -Inf, upper = Inf)$value
  } else {
    kern_func <- function(x) (abs(x) <= 1) * 0.75 * (1 - x^2)
    sigma_K <- 0.2
    two_norm_sq <- 0.6
  }
  
  if (!is.numeric(bw)) {
    if (bw == "amise") {
      bw <- amise(x)
    } else if (bw == "cv") {
      bw <- cv(x)
    } else if (bw == "loocv") {
      bw <- loocv(x)
    }
  }
  
  rg <- range(x)
  xx <- seq(rg[1] - 3 * bw, rg[2] + 3 * bw, length.out = m)
  y <- numeric(m) 
  bw <- bw / sqrt(sigma_K)
  
  if(equi_bin) {
    breaks <- seq(rg[1], rg[2], length.out = 250)
  } else {
    breaks <- my_breaks2(x, 8)
  }
  counts <- count_my_breaks(x, breaks)
  
  if(mid) {
    midpoints <- breaks[1:(length(breaks)-1)] + abs(breaks[2:length(breaks)] - breaks[1:(length(breaks)-1)]) / 2
  } else {
    midpoints <- mean_in_breaks(lF12, counts)
  }
  
  if(binning) {
    for (i in seq_along(xx)) {
      y[i] <- sum(counts * kern_func((xx[i] - midpoints) / bw))
    }
  } else {
    for (i in seq_along(xx)) {
      y[i] <- sum(kern_func((xx[i] - x) / bw))
    }
  }
  y <- y / (bw * length(x))
  
  structure(list(x = xx, y = y), class = "My_density_class")
}

plot.My_density_class <- function(dens, add = FALSE, ...) {
  if (add) {
    lines(dens$x, dens$y, ...)
  } else {
    ggplot2::qplot(dens$x, dens$y, geom="line", ...) + 
      theme_bw() + xlab("x") + ylab("density") + theme(legend.position = "none")
  }
}
