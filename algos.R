## Gradient descent ------------------
GD <- function(
  par, 
  H,
  gr,
  d = 0.8,
  c = 0.2, 
  gamma0 = 0.1,
  epsilon = 1e-4, 
  maxiter = 3000,
  cb = NULL,
  ...
) {
  for(i in 1:maxiter) {
    value <- H(par, ...)
    grad <- gr(par, ...)
    h_prime <- sum(grad^2)
    if(!is.null(cb)) cb()
    # Convergence criterion based on gradient norm
    if(h_prime <= epsilon) break
    gamma <- gamma0
    # Proposed descent step
    par1 <- par - gamma * grad
    # Backtracking while descent is insufficient
    while(H(par1, ...) > value - c * gamma * h_prime) {
      gamma <- d * gamma
      par1 <- par - gamma * grad
    }
    par <- par1
  }
  par
}

## Variations of SG ---------------
SG <- function(
  par, 
  N,                 # Sample size
  gamma,             # Decay schedule or a fixed learning rate
  epoch = batch,     # Epoch update function
  ...,               # Other arguments passed to epoch updates
  maxiter = 1000,     # Max epoch iterations
  sampler = sample,  # How data is resampled. Default is a random permutation
  cb = NULL
) {
  gamma <- if (is.function(gamma)) gamma(1:maxiter) else rep(gamma, maxiter) 
  for(k in 1:maxiter) {
    samp <- sampler(N)
    if(!is.null(cb)) cb()
    par <- epoch(par, samp, gamma[k], ...)
  }
  par
}

## Mini-batch function
batch <- function(par,
                  samp,
                  gamma,
                  grad,              # Function of parameter and observation index
                  m = 100,            # Mini-batch size
                  ...) {
  M <- floor(length(samp) / m)
  for (j in 0:(M - 1)) {
    i <- samp[(j * m + 1):(j * m + m)]
    par <- par - gamma * grad(par, i, ...)
  }
  par
}

## Momentum function 
momentum <- function() {
  rho <- 0
  function(
    par, 
    samp,
    gamma,  
    grad,
    m = 100,             # Mini-batch size
    beta = 0.95,        # Momentum memory 
    ...
  ) {
    M <- floor(length(samp) / m)
    for(j in 0:(M - 1)) {
      i <- samp[(j * m + 1):(j * m + m)]
      # Using '<<-' assigns the value to rho in the enclosing environment
      rho <<- beta * rho + (1 - beta) * grad(par, i, ...)  
      par <- par - gamma * rho
    }
    par
  }
}

## Adam function
adam <- function() {
  rho <- v <- 0
  function(
    par, 
    samp,
    gamma,   
    grad,
    m = 100,         # Mini-batch size
    beta1 = 0.9,     # Momentum memory     
    beta2 = 0.9,     # Second moment memory 
    ...
  ) {
    M <- floor(length(samp) / m)
    
    for(j in 0:(M - 1)) {
      i <- samp[(j * m + 1):(j * m + m)]
      gr <- grad(par, i, ...) - rho*beta1
      rho <<- beta1 * rho + (1 - beta1) * gr
      v <<- beta2 * v + (1 - beta2) * gr^2
      par <- par - gamma * (rho / (sqrt(v) + 1e-8))  
    }
    par
  }
}