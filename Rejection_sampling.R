envelope_vec <- function(p, func, limits, deriv = NULL, log_concave = FALSE) {
  if (is.null(deriv)) {
    deriv <- function(x)
      (func(x + 1e-10) - func(x - 1e-10)) / 2e-10
  }
  p <- unique(sort(p))
  p[p == limits[1]] <- p[p == limits[1]] + 1e-10
  p[p == limits[2]] <- p[p == limits[2]] - 1e-10
  n <- length(p)
  minimum <- min(limits)
  maximum <- max(limits)
  
  if (log_concave) {
    a <- deriv(p) / func(p)
    b <- log(func(p)) - a*p
    s <- if(n > 1) (b[1:(n-1)] - b[2:n]) / (a[2:n] - a[1:(n-1)]) else c()
  } else {
    a <- deriv(p)
    b <- func(p) - a*p
    s <- if(n > 1) (b[1:(n-1)] - b[2:n]) / (a[2:n] - a[1:(n-1)]) else c()
  }
  
  s <- c(minimum, s, maximum)
  
  env <- function(x, log = log_concave) {
    n_x <- length(x)
    interval <- numeric(n_x)
    interval[x == minimum] <- 1
    
    valid_range <- as.logical((x >= minimum) * (x <= maximum))
    
    for (j in 1:n) {
      i <- as.logical((x > s[j]) * (x <= s[j+1]))
      interval[i] <- j
    }
    
    ret <- numeric(n_x)
    if(log_concave) {
      ret[valid_range] <- exp(x[valid_range] * a[interval] +  b[interval])
    } else {
      ret[valid_range] <- x[valid_range] * a[interval] +  b[interval]
    }
    ret[!valid_range] <- 0
    return(ret)
    
  }
}

envelope_inv_vec <- function(p, func, limits, deriv = NULL, log_concave = FALSE) {
  if (is.null(deriv)) {
    deriv <- function(x) {
      (func(x + 0.001) - func(x - 0.001)) / 0.002
    }
  }
  p <- unique(sort(p))
  p[p == limits[1]] <- p[p == limits[1]] + 1e-10
  p[p == limits[2]] <- p[p == limits[2]] - 1e-10
  n <- length(p)
  minimum <- min(limits)
  maximum <- max(limits)
  
  if (log_concave) {
    a <- deriv(p) / func(p)
    b <- log(func(p)) - a*p
    s <- if(n > 1) (b[1:(n-1)] - b[2:n]) / (a[2:n] - a[1:(n-1)]) else c()
  } else {
    a <- deriv(p)
    b <- func(p) - a*p
    s <- if(n > 1) (b[1:(n-1)] - b[2:n]) / (a[2:n] - a[1:(n-1)]) else c()
  }
  
  if(sum(is.na(a)) > 0 || sum(sum(is.na(b))) > 0) {
    warning("Target density opfylder ikke krav for at være en tæthed, eller den valgte support er forkert")
    stop()
  }
  s <- c(minimum, s, maximum)
  
  integ <- function(a, b, from, to) {
    if(log_concave) {
      if(a != 0) {
        (exp(a * to + b) - exp(a * from + b)) / a
      } else {
        exp(b) * (to - from)
      }
    } else {
      (a / 2) * (to^2 - from^2) + b * (to - from)
    }
  }
  
  integrals <- numeric(n)
  for (i in 1:n) {
    integrals[i] <- integ(a[i], b[i], s[i], s[i+1])
  }
  
  env_inv <- function(x) {
    if(sum(x > 1) > 0 || sum(x < 0) > 0) warning(paste("One of more points outside the scope. Data must be in the range from 0 to 1."))
    n_x <- length(x)
    x <- x * sum(integrals)
    
    integrals_cum <- cumsum(c(0, integrals))
    
    interval <- numeric(n_x)
    interval[x == 0] <- 1
    for (j in 1:n) {
      i <- as.logical((x > integrals_cum[j]) * (x <= integrals_cum[j+1]))
      interval[i] <- j
    }
    
    x <- x - integrals_cum[interval]
    from <- s[-(n+1)]
    
    a_is_zero <- if(any(a == 0)) which(a == 0) else 0
    
    interval_a_n0 <- interval[interval != a_is_zero]
    interval_a_0 <- interval[interval == a_is_zero]
    x_a_n0 <- x[interval != a_is_zero]
    x_a_0 <- x[interval == a_is_zero]
    
    closed_form <- numeric(n_x)
    
    if(log_concave) {
      closed_form[interval != a_is_zero] <- (log(x_a_n0 * a[interval_a_n0] + exp(a[interval_a_n0] * from[interval_a_n0] + b[interval_a_n0])) - b[interval_a_n0]) / a[interval_a_n0]
      closed_form[interval == a_is_zero] <- x_a_0 * exp(-b[interval_a_0]) + from[interval_a_0]
    } else {
      c1 <- a[interval_a_n0] / 2
      c2 <- b[interval_a_n0]
      c3 <- - (a[interval_a_n0] / 2 * from[interval_a_n0]^2 + b[interval_a_n0] * from[interval_a_n0] + x_a_n0)
      
      closed_form[interval != a_is_zero] <- (-c2 + sqrt(c2^2 - 4 * c1 * c3)) / (2 * c1)
      closed_form[interval == a_is_zero] <- (x_a_0 + from[interval_a_0] * b[interval_a_0]) / b[interval_a_0]
    }
    
    return(closed_form)
  }
}

reverse_envelope <- function(p, func, limits, log_concave = FALSE) {
  minimum <- min(limits)
  maximum <- max(limits)
  p <- unique(sort(c(minimum, p, maximum)))
  n <- length(p)
  
  if (log_concave) {
    a <- (log(func(p[2:n])) - log(func(p[1:(n-1)]))) / (p[2:n] - p[1:(n-1)])
    b <- log(func(p[1:(n-1)])) - a * p[1:(n-1)]
  } else {
    a <- (func(p[2:n]) - func(p[1:(n-1)])) / (p[2:n] - p[1:(n-1)])
    b <- func(p[1:(n-1)]) - a * p[1:(n-1)]
  }
  
  
  env <- function(x, log = log_concave) {
    n_x <- length(x)
    interval <- numeric(n_x)
    interval[x == minimum] <- 1
    
    valid_range <- as.logical((x >= minimum) * (x <= maximum))
    
    for (j in 1:(n-1)) {
      i <- as.logical((x > p[j]) * (x <= p[j+1]))
      interval[i] <- j
    }
    ret <- numeric(n_x)
    if(log_concave) {
      ret[valid_range] <- exp(x[valid_range] * a[interval] +  b[interval])
    } else {
      ret[valid_range] <- x[valid_range] * a[interval] +  b[interval]
    }
    ret[!valid_range] <- 0
    return(ret)
    
  }
}

target_dens_gen <- function(target_dens, eval_points, limits, deriv, log_concave) {
  structure(list(
    target_dens = target_dens,
    eval_points = eval_points,
    limits = limits,
    deriv = deriv,
    log_concave = log_concave
  ), class = "rejection_target_dens_class")
}

rgen <- function(target_dens_class) {
  
  target_dens <- target_dens_class$target_dens
  limits <- target_dens_class$limits
  
  if(!is.numeric(target_dens_class$eval_points)) {
    eval_p <- numeric(8)
    
    if(limits[1] == -Inf && limits[2] == Inf) {
      eval_p[1] <- runif(1, min = -1000, max = 1000)
    } else if(limits[1] == -Inf) {
      eval_p[1] <- runif(1, min = -1000, max = limits[2])
    } else if(limits[2] == Inf) {
      eval_p[1] <- runif(1, min = limits[1], max = 1000)
    } else {
      eval_p[1] <- runif(1, min = limits[1], max = limits[2])
    }
    
    for (i in 2:8) {
      proposal_F <- envelope_inv_vec(
        p = eval_p[1:(i-1)],
        func = target_dens,
        limits = limits,
        deriv = target_dens_class$deriv,
        log_concave = target_dens_class$log_concave
      )
      eval_p[i] <- proposal_F(runif(1))
    }
    
  } else {
    eval_p <- target_dens_class$eval_points
  }
  
  # Invers af integrerede adaptive envelope, F^-1, til at simulere fra adaptive envelope
  proposal_F <- envelope_inv_vec(
    p = eval_p,
    func = target_dens,
    limits = limits,
    deriv = target_dens_class$deriv,
    log_concave = target_dens_class$log_concave
  )
  
  # Adaptive envelope
  proposal <- envelope_vec(
    p = eval_p,
    func = target_dens,
    limits = limits,
    deriv = target_dens_class$deriv,
    log_concave = target_dens_class$log_concave
  )
  
  target_int <- integrate(target_dens, lower = limits[1], upper = limits[2])$value
  alpha <- target_int / sum(environment(proposal_F)$integrals) 
  
  ran_gen <- function(n) {
    n_accept <- 0
    ret <- list()
    i <- 1
    while(n_accept < n) {
      u <- matrix(runif(10 + 2 * ceiling((n - n_accept) / alpha)), ncol = 2)
      y <- proposal_F(u[, 1])
      accepts <- u[, 2] * proposal(y) <= target_dens(y)
      n_accept <- n_accept + sum(accepts)
      ret[[i]] <- y[accepts]
      i <- i + 1
    }
    return(unlist(ret)[1:n])
  }
  
}
rgen_cache <- function(target_dens_class) {
  
  target_dens <- target_dens_class$target_dens
  limits <- target_dens_class$limits
  
  if(!is.numeric(target_dens_class$eval_points)) {
    eval_p <- numeric(8)
    
    if(limits[1] == -Inf && limits[2] == Inf) {
      eval_p[1] <- runif(1, min = -1000, max = 1000)
    } else if(limits[1] == -Inf) {
      eval_p[1] <- runif(1, min = -1000, max = limits[2])
    } else if(limits[2] == Inf) {
      eval_p[1] <- runif(1, min = limits[1], max = 1000)
    } else {
      eval_p[1] <- runif(1, min = limits[1], max = limits[2])
    }
    
    for (i in 2:8) {
      proposal_F <- envelope_inv_vec(
        p = eval_p[1:(i-1)],
        func = target_dens,
        limits = limits,
        deriv = target_dens_class$deriv,
        log_concave = target_dens_class$log_concave
      )
      eval_p[i] <- proposal_F(runif(1))
    }
    
  } else {
    eval_p <- target_dens_class$eval_points
  }
  
  # Invers af integrerede adaptive envelope, F^-1, til at simulere fra adaptive envelope
  proposal_F <- envelope_inv_vec(
    p = eval_p,
    func = target_dens,
    limits = limits,
    deriv = target_dens_class$deriv,
    log_concave = target_dens_class$log_concave
  )
  
  # Adaptive envelope
  proposal <- envelope_vec(
    p = eval_p,
    func = target_dens,
    limits = limits,
    deriv = target_dens_class$deriv,
    log_concave = target_dens_class$log_concave
  )
  
  # target_int <- integrate(target_dens, lower = limits[1], upper = limits[2])$value
  # alpha <- target_int / sum(environment(proposal_F)$integrals)
  
  u <- matrix(runif(10000), ncol = 2)
  y <- proposal_F(u[, 1])
  accepts <- u[, 2] * proposal(y) <= target_dens(y)
  alpha <- sum(accepts) / 5000
  
  cache <- y[accepts]
  cache_length <- length(cache)
  
  ran_gen <- function(n) {
    if (cache_length >= n) {
      ret <- cache[1:n]
      cache <<- cache[-(1:n)]
      cache_length <<- cache_length - n
      return(ret)
    }
    n_accept <- cache_length
    ret <- list(cache)
    i <- 1
    while(n_accept < n) {
      u <- matrix(runif(10 + 2 * ceiling((n - n_accept) / alpha)), ncol = 2)
      y <- proposal_F(u[, 1])
      accepts <- u[, 2] * proposal(y) <= target_dens(y)
      n_accept <- n_accept + sum(accepts)
      ret[[i]] <- y[accepts]
      i <- i + 1
    }
    unlisted_ret <- unlist(ret)
    cache <<- unlisted_ret[-(1:n)]
    cache_length <<- length(cache)
    return(unlisted_ret[1:n])
  }
  
}

gaussian_rejection_sampling_gen <- function(target_dens, sig_grid_search = seq(0.01, 0.21, 0.01)) {
  proposal_dens <- function(x, m, std, scale) {
    dnorm(x, mean = m, sd = std) / scale
  }
  m <- optimize(target_dens, c(0, 1), maximum = TRUE)$maximum
  std_v <- sig_grid_search
  
  obj_v <- numeric(length(std_v))
  for (i in seq_along(std_v)) {
    obj_v[i] <- optimize(function(x) proposal_dens(x, m, std_v[i], 1) / target_dens(x), c(0,m))$objective
    if(i != 1 && obj_v[i] < obj_v[i-1]) break
  }
  std <- std_v[which.max(obj_v)]
  obj <- obj_v[which.max(obj_v)]
  
  alpha <- integrate(target_dens, lower = -Inf, upper = Inf)$value * obj
  
  rejection_sampling <- function(n) {
    n_accept <- 0
    ret <- list()
    i <- 1
  
    while(n_accept < n) {
      k <- ceiling((n - n_accept) / alpha)
      u <- runif(k)
      y <- rnorm(k, mean = m, sd = std)
      accepts <- u * proposal_dens(y, m, std, obj) <= target_dens(y)
      n_accept <- n_accept + sum(accepts)
      ret[[i]] <- y[accepts]
      i <- i + 1
    }
    return(unlist(ret)[1:n])
  }
}
