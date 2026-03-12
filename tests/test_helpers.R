# =============================================================================
# tests/test_helpers.R
# Minimal test harness: no external dependencies required.
# =============================================================================

.test_env <- new.env()
.test_env$results <- list()

# -- Assertions ---------------------------------------------------------------

expect_true <- function(condition, msg = "") {
  if (!isTRUE(condition)) stop(paste("Assertion failed:", msg))
}

expect_equal <- function(actual, expected, tol = 1e-6, msg = "") {
  if (abs(actual - expected) > tol)
    stop(sprintf("Expected %.6g, got %.6g. %s", expected, actual, msg))
}

expect_less_than <- function(actual, threshold, msg = "") {
  if (actual >= threshold)
    stop(sprintf("%.6g >= threshold %.6g. %s", actual, threshold, msg))
}

expect_greater_than <- function(actual, threshold, msg = "") {
  if (actual <= threshold)
    stop(sprintf("%.6g <= threshold %.6g. %s", actual, threshold, msg))
}

expect_no_na <- function(x, msg = "") {
  if (any(is.na(x))) stop(paste("Unexpected NAs found.", msg))
}

# -- Test runner --------------------------------------------------------------

test_that <- function(description, code) {
  result <- tryCatch({
    code
    list(name = description, passed = TRUE, message = "")
  }, error = function(e) {
    list(name = description, passed = FALSE, message = conditionMessage(e))
  })

  tag <- if (result$passed) "\033[32mPASS\033[0m" else "\033[31mFAIL\033[0m"
  cat(sprintf("  [%s] %s\n", tag, result$name))
  if (!result$passed)
    cat(sprintf("         -> %s\n", result$message))

  .test_env$results[[length(.test_env$results) + 1]] <- result
  invisible(result$passed)
}

describe <- function(suite_name, code) {
  cat(sprintf("\n=== %s ===\n", suite_name))
  code
}

# -- Summary ------------------------------------------------------------------

test_summary <- function() {
  results <- .test_env$results
  n_pass <- sum(sapply(results, `[[`, "passed"))
  n_fail <- length(results) - n_pass
  cat(sprintf(
    "\n--- Summary: %d passed, %d failed (total %d) ---\n",
    n_pass, n_fail, length(results)
  ))
  invisible(list(passed = n_pass, failed = n_fail))
}

# -- Gradient checker ---------------------------------------------------------
# Compares an analytical gradient against central finite differences.
# Returns a list with max_rel_error and cosine_similarity.

check_gradient <- function(fn, gr, par, eps = 1e-5, ...) {
  g_analytic  <- gr(par, ...)
  n           <- length(par)
  g_numerical <- numeric(n)

  for (i in seq_len(n)) {
    p_plus  <- par; p_plus[i]  <- par[i] + eps
    p_minus <- par; p_minus[i] <- par[i] - eps
    g_numerical[i] <- (fn(p_plus, ...) - fn(p_minus, ...)) / (2 * eps)
  }

  rel_err  <- abs(g_analytic - g_numerical) / (abs(g_numerical) + 1e-10)
  norm_a   <- sqrt(sum(g_analytic^2))
  norm_n   <- sqrt(sum(g_numerical^2))
  cos_sim  <- if (norm_a > 0 && norm_n > 0)
    sum(g_analytic * g_numerical) / (norm_a * norm_n) else NA

  list(
    max_rel_error      = max(rel_err),
    cosine_similarity  = cos_sim,
    g_analytic         = g_analytic,
    g_numerical        = g_numerical
  )
}
