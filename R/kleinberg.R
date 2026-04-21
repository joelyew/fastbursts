#' Optimized Kleinberg Burst Detection
#'
#' @param offsets A numeric vector of event times.
#' @param s The base of the exponent used to determine event frequencies in a given state.
#' @param gamma A coefficient that modifies the cost of a transition to a higher state.
#' @return A data frame of class 'bursts'
#' @export
kleinberg <- function(offsets, s = 2, gamma = 1) {
  if (s <= 1) stop("s must be greater than 1!")
  if (gamma <= 0) stop("gamma must be positive!")
  if (length(offsets) < 1) stop("offsets must be non-empty!")

  if (length(offsets) == 1) {
    bursts <- data.frame(level = 0, start = offsets[1], end = offsets[1])
    class(bursts) <- c("bursts", "data.frame")
    return(bursts)
  }

  offsets <- sort(offsets)
  gaps <- as.numeric(diff(offsets))

  if (any(gaps == 0)) {
    stop("Input cannot contain events with zero time between!")
  }

  T_sum <- sum(gaps)
  n <- length(gaps)
  ghat <- T_sum / n
  
  k <- ceiling(1 + log(T_sum, s) + log(1/min(gaps), s))
  gammalogn <- gamma * log(n)
  
  alpha <- (s^(0:(k - 1))) / ghat
  
  tau_mat <- matrix(0, nrow = k, ncol = k)
  for (j in 1:k) {
    for (i in 1:k) {
      if (i < j) {
        tau_mat[i, j] <- (j - i) * gammalogn
      }
    }
  }

  C <- c(0, rep(Inf, k - 1))
  B <- matrix(NA_integer_, nrow = k, ncol = n)
  
  for (t in 1:n) {
    Cprime <- rep(NA_real_, k)
    log_f_val <- log(alpha) - alpha * gaps[t]
    
    for (j in 1:k) {
      costs_to_j <- C + tau_mat[, j]
      best_ell <- which.min(costs_to_j)
      B[j, t] <- best_ell
      Cprime[j] <- costs_to_j[best_ell] - log_f_val[j]
    }
    C <- Cprime
  }

  q <- integer(n)
  q[n] <- which.min(C)
  if (n > 1) {
    for (t in n:2) {
      q[t - 1] <- B[q[t], t]
    }
  }

  prev_q <- 0
  N_bursts <- 0
  for (t in 1:n) {
    if (q[t] > prev_q) N_bursts <- N_bursts + q[t] - prev_q
    prev_q <- q[t]
  }

  bursts <- data.frame(
    level = rep(NA_integer_, N_bursts), 
    start = rep(offsets[1], N_bursts), 
    end = rep(offsets[1], N_bursts)
  )

  burstcounter <- 0
  prev_q <- 0
  stack <- rep(NA_integer_, N_bursts)
  stackcounter <- 0
  
  for (t in 1:n) {
    if (q[t] > prev_q) {
      num_levels_opened <- q[t] - prev_q
      for (i in 1:num_levels_opened) {
        burstcounter <- burstcounter + 1
        bursts$level[burstcounter] <- prev_q + i
        bursts$start[burstcounter] <- offsets[t]
        stackcounter <- stackcounter + 1
        stack[stackcounter] <- burstcounter
      }
    } else if (q[t] < prev_q) {
      num_levels_closed <- prev_q - q[t]
      for (i in 1:num_levels_closed) {
        bursts$end[stack[stackcounter]] <- offsets[t]
        stackcounter <- stackcounter - 1
      }
    }
    prev_q <- q[t]
  }

  while (stackcounter > 0) {
    bursts$end[stack[stackcounter]] <- offsets[n + 1]
    stackcounter <- stackcounter - 1
  }

  class(bursts) <- c("bursts", "data.frame")
  return(bursts)
}

#' Plot Burst Model
#' 
#' @param x A 'bursts' object.
#' @param ... Additional arguments to plot.
#' @export
plot.bursts <- function(x, ...) {
  if (!inherits(x, "bursts")) stop("x must be a 'bursts' object")
  
  graphics::plot(c(min(x$start), max(x$end)), c(0, max(x$level) + 1), 
                 type = "n", xlab = "Time", ylab = "Level", ...)
  
  for (i in 1:nrow(x)) {
    graphics::rect(x$start[i], x$level[i] - 0.4, x$end[i], x$level[i] + 0.4, ...)
  }
}
