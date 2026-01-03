# polyenclose.R  ---------------------------------------------------------------
# Tests / examples (symbolic only: feasibility of face-multisets; no embedding)
#
# A) Vertex-range sanity (external + internal ladders)
#   V_report(8)                    # 8-vertex genus-0 surfaces: E/F/S ranges + internal decomposition ladder: T/Ni/Si
#   decompose(8)                   # Internal decomposition ladder only: T/Ni/Si
#
# B) Canonical face-multisets (known polyhedra; should pass symbolic checks)
#   assess("6x4")                  # Cube: 6 squares
#   assess("12x5")                 # Dodecahedron: 12 pentagons
#
# C) “House / dome / capped” mixed faces (symbolically plausible combos; not unique geometrically)
#   assess("4x3+5x4")            # “Capped cube / house-like”: 5 squares + 4 triangles (brick-house vibe)
#                                # Interpretation: cube-like body with a pyramidal cap replacing one square
# D) 36 hexagons
#   assess("36x6")            
#                                # Interpretation: sets of hexahons are too flat to enclose
#
# E) Enumeration / exploration (face-type distributions for fixed V,S)
#   Pk_wizard(V=20, S=24)        # Count/enumerate Pk solutions at (V,S)
#   Pk_wizard(V=20, S=24, k_max=5) # Same, restricting to face degrees ≤ 5
#
#
# F) Summarizes the solution space of admissible face-type configurations
#   out <- Pk_wizard(V=20, S=18, k_max=6, show_solutions=FALSE)
#   Pk_summary(out)
#
# # G) Loose flatness-only upper bound on admissible face-type configurations as function of V (partition growth)
#   face_multiset_upperbound(20)
#
#
# =============================================================================
# Block A — Helpers: validate + parse Pk
# =============================================================================

#' Validate set of polygons can enclose genus-0 polyhedron (input hygiene)
pk_validate <- function(Pk) {
  if (length(Pk) == 0) stop("Pk is empty.")
  if (any(is.na(Pk))) stop("Pk contains NA.")
  if (any(Pk < 0)) stop("Pk contains negative counts.")
  if (any(Pk != as.integer(Pk))) stop("Pk counts must be integers.")
  k <- suppressWarnings(as.integer(names(Pk)))
  if (any(is.na(k))) stop("Pk names must be integer face degrees (e.g., '3','4','5').")
  if (any(k < 3)) stop("Face degrees must be >= 3.")
  Pk <- Pk[Pk != 0]
  if (length(Pk) == 0) stop("All Pk counts are zero.")
  Pk
}

#' Parse Pk from named vector, data.frame(k,Pk), or string like "4x3+6x4"
pk_parse <- function(x) {
  if (is.null(x)) stop("Pk input is NULL.")
  
  if (is.numeric(x) && !is.null(names(x))) {
    Pk <- as.integer(x); names(Pk) <- as.character(names(x))
    return(pk_validate(Pk))
  }
  
  if (is.data.frame(x)) {
    if (!all(c("k","Pk") %in% names(x))) stop("data.frame input must have columns 'k' and 'Pk'.")
    out <- stats::setNames(as.integer(x$Pk), as.character(as.integer(x$k)))
    return(pk_validate(out))
  }
  
  if (is.character(x) && length(x) == 1) {
    s <- gsub("\\s+", "", x)
    s <- gsub("^Pk=\\{", "", s, ignore.case = TRUE)
    s <- gsub("^\\{", "", s)
    s <- gsub("\\}$", "", s)
    s <- gsub(";", ",", s)
    
    pieces <- unlist(strsplit(s, "\\+|,"))
    pieces <- pieces[pieces != ""]
    
    k_list <- integer(0)
    c_list <- integer(0)
    
    for (p in pieces) {
      m1 <- regexec("^([0-9]+)(x|\\*|@)([0-9]+)$", p, ignore.case = TRUE)
      r1 <- regmatches(p, m1)[[1]]
      if (length(r1) == 4) {
        c_list <- c(c_list, as.integer(r1[2]))  # count
        k_list <- c(k_list, as.integer(r1[4]))  # degree
        next
      }
      m2 <- regexec("^([0-9]+)(:|=)([0-9]+)$", p)
      r2 <- regmatches(p, m2)[[1]]
      if (length(r2) == 4) {
        k_list <- c(k_list, as.integer(r2[2]))  # degree
        c_list <- c(c_list, as.integer(r2[4]))  # count
        next
      }
      stop(sprintf("Could not parse piece '%s'. Try '12x5+20x6' or '5:12,6:20'.", p))
    }
    
    out <- tapply(c_list, k_list, sum)
    
    # preserve names explicitly before coercion
    nm <- names(out)
    out <- as.integer(out)
    names(out) <- nm
    
    # normalize + sort names
    names(out) <- as.character(as.integer(names(out)))
    out <- out[order(as.integer(names(out)))]
    
    return(pk_validate(out))
  }
  
  stop("Unsupported Pk input. Use a named numeric vector, a data.frame(k,Pk), or a single string.")
}

# =============================================================================
# Block B — Core symbolic math: V_report (+ printer)
# =============================================================================

#' Report external (E,F,S) and internal (T,Ni,Si) ranges from vertex count V
#'
#' @param V integer >= 4
#' @return list of class 'V_report' with two data.frames: external and internal
#' @export
V_report <- function(V) {
  V <- as.integer(V)
  if (length(V) != 1 || is.na(V) || V < 4) stop("V must be a single integer >= 4.")
  
  # ---- External: triangulated maxima ----
  E_max <- 3*V - 6
  F_max <- 2*V - 4
  
  # ---- External: parity-aware minima via trivalence lower bound ----
  E_min <- ceiling(3*V/2)
  F_min <- 2 - V + E_min  # from Euler: V - E + F = 2
  
  # ---- External: S range induced by E = E_max - S ----
  S_max <- E_max - E_min
  if (S_max < 0) S_max <- 0L  # safety; should not happen for V>=4
  
  S_vals <- 0:S_max
  external <- data.frame(
    V = V,
    S = S_vals,
    E = E_max - S_vals,
    F = F_max - S_vals
  )
  
  # ---- Internal: SALT+MIE ranges + ladder (as in Section 5) ----
  # SALT+MIE empirical bound:
  #   0 <= S_i <= V-5
  # ladder relations:
  #   T = V-3 + S_i
  #   N_i = V-4 + 2 S_i
  
  Si_max <- V - 5L
  if (Si_max < 0L) {
    internal <- data.frame()
    internal_bounds <- list(
      T_min = NA_integer_, T_max = NA_integer_,
      Ni_min = NA_integer_, Ni_max = NA_integer_,
      Si_min = NA_integer_, Si_max = NA_integer_
    )
  } else {
    Si_vals <- 0L:Si_max
    
    # ladder family (exact admissible triples under SALT+MIE)
    T_vals  <- (V - 3L) + Si_vals
    Ni_vals <- (V - 4L) + 2L * Si_vals
    
    internal <- data.frame(
      V  = V,
      k  = Si_vals,     # ladder index (equals S_i)
      T  = T_vals,
      Ni = Ni_vals,
      Si = Si_vals
    )
    
    internal_bounds <- list(
      T_min  = V - 3L,
      T_max  = 2L * (V - 4L),
      Ni_min = V - 4L,
      Ni_max = 3L * (V - 4L) - 2L,
      Si_min = 0L,
      Si_max = Si_max
    )
  }
  
  out <- list(
    V = V,
    external = external,
    internal = internal,
    bounds = list(
      E_min = E_min, E_max = E_max,
      F_min = F_min, F_max = F_max,
      S_min = 0L,    S_max = S_max,
      internal = internal_bounds
    )
  )
  
  class(out) <- "V_report"
  out
}

#' @export
print.V_report <- function(x, ...) {
  cat(sprintf("V_report (symbolic ranges from V=%d)\n\n", x$V))
  
  b <- x$bounds
  cat("External bounds:\n")
  cat(sprintf("  E in [%d, %d],  F in [%d, %d],  S in [%d, %d]\n\n",
              b$E_min, b$E_max, b$F_min, b$F_max, b$S_min, b$S_max))
  
  ib <- b$internal
  cat("Internal bounds (SALT+MIE):\n")
  cat(sprintf("  T in [%s, %s],  Ni in [%s, %s],  Si in [%s, %s]\n\n",
              ib$T_min, ib$T_max, ib$Ni_min, ib$Ni_max, ib$Si_min, ib$Si_max))
  
  cat("External surface table (S, E(S), F(S)):\n")
  print(x$external, row.names = FALSE)
  
  cat("\nInternal (SALT+MIE ladder; ranges are coupled):\n")
  print(x$internal, row.names = FALSE)
  
  cat("\nRemark.\n")
  cat("The upper endpoint (T = T_max) corresponds only to bipyramidal configurations.\n")
  cat("This extremal case is valid regardless of coplanarity of equatorial vertices\n")
  cat("and regardless of the geometric alignment of the two apices.\n")
  
  invisible(x)
}

# =============================================================================
# Block B2 — Internal-only entry point: decompose(V)
# =============================================================================

#' Internal decomposition ladder only (SALT+MIE)
#'
#' Convenience wrapper around V_report(V) that prints ONLY the internal
#' decomposition ladder (T/Ni/Si) and its coupled bounds.
#'
#' @param V integer >= 4
#' @return object of class 'decompose' (invisibly printed by print.decompose)
#' @export
decompose <- function(V) {
  vr <- V_report(V)
  out <- list(
    V = vr$V,
    internal = vr$internal,
    internal_bounds = vr$bounds$internal
  )
  class(out) <- "decompose"
  out
}

#' @export
print.decompose <- function(x, ...) {
  cat(sprintf("decompose(V=%d)\n\n", x$V))
  
  ib <- x$internal_bounds
  cat("Internal bounds (SALT+MIE):\n")
  cat(sprintf("  T in [%s, %s],  Ni in [%s, %s],  Si in [%s, %s]\n\n",
              ib$T_min, ib$T_max, ib$Ni_min, ib$Ni_max, ib$Si_min, ib$Si_max))
  
  cat("Internal decomposition ladder (ranges are coupled):\n")
  print(x$internal, row.names = FALSE)
  
  cat("\nRemark.\n")
  cat("This ladder enumerates the admissible triples (T, Ni, Si) under SALT+MIE;\n")
  cat("it does not depend on any embedding, only the vertex count V.\n")
  cat("Maximal T corresponds ONLY to bipyramidal configurations, realizable at S=0.\n")
  
  invisible(x)
}


# =============================================================================
# Block C — Pk invariants + symbolic feasibility report (Section 7 filter)
# =============================================================================

pk_invariants <- function(Pk) {
  Pk <- pk_parse(Pk)
  k <- as.integer(names(Pk))
  F <- sum(Pk)
  N <- sum((k - 2L) * Pk)
  S_total <- sum((k - 3L) * Pk)
  
  V <- if (N %% 2L == 0L) (N/2L + 2L) else NA_integer_
  E <- if (!is.na(V)) (3L * (V - 2L) - S_total) else NA_real_
  
  list(k = k, Pk = Pk, N = N, V = V, F = F, S_total = S_total, E = E)
}

assess <- function(Pk) {
  
  # ---- dispatch: allow assess("V=8") ----
  if (is.character(Pk) && length(Pk) == 1) {
    s <- gsub("\\s+", "", Pk)
    if (grepl("^V=\\d+$", s, ignore.case = TRUE)) {
      V <- as.integer(sub("^V=", "", s, ignore.case = TRUE))
      return(V_report(V))
    }
  }
  
  inv <- pk_invariants(Pk)
  integer_ok <- !is.na(inv$V)
  
  euler_ok <- if (integer_ok) (inv$V - inv$E + inv$F == 2) else FALSE
  sum_kPk <- sum(inv$k * inv$Pk)
  incidence_ok <- if (integer_ok) (2 * inv$E == sum_kPk) else FALSE
  
  # ---- UNIFY S_max: use V_report parity-aware S_max ----
  vr <- if (integer_ok) V_report(inv$V) else NULL
  Smax <- if (integer_ok) vr$bounds$S_max else NA_real_
  Sdiff <- if (integer_ok) (Smax - inv$S_total) else NA_real_
  flatness_ok <- if (integer_ok) (inv$S_total >= 0 && inv$S_total <= Smax) else FALSE
  
  feasible_symbolic <- integer_ok && euler_ok && incidence_ok && flatness_ok
  
  out <- list(
    input = inv$Pk,
    invariants = inv,
    checks = list(
      integer_ok = integer_ok,
      euler_ok = euler_ok,
      incidence_ok = incidence_ok,
      flatness_ok = flatness_ok
    ),
    S_max = Smax,
    S_difference = Sdiff,
    feasible_symbolic = feasible_symbolic
  )
  class(out) <- "assess"
  out
}

print.assess <- function(x, ...) {
  cat("polyenclose report (Section 7 symbolic filter)\n")
  df <- data.frame(k = as.integer(names(x$input)), Pk = as.integer(x$input))
  df <- df[order(df$k), , drop = FALSE]
  cat("Input Pk:\n"); print(df, row.names = FALSE)
  
  inv <- x$invariants
  cat("\nInvariants:\n")
  cat(sprintf("  N = %s\n", inv$N))
  cat(sprintf("  V = %s\n", ifelse(is.na(inv$V), "NA (N odd)", inv$V)))
  cat(sprintf("  F = %s\n", inv$F))
  cat(sprintf("  S_total = %s\n", inv$S_total))
  cat(sprintf("  E = %s\n", ifelse(is.na(inv$E), "NA", inv$E)))
  cat(sprintf("  S_max(V) = %s\n", ifelse(is.na(x$S_max), "NA", x$S_max)))
  cat(sprintf("  S_max - S_total = %s\n", ifelse(is.na(x$S_difference), "NA", x$S_difference)))
  
  cat("\nChecks:\n")
  ch <- x$checks
  cat(sprintf("  integer_ok   : %s\n", ch$integer_ok))
  cat(sprintf("  euler_ok     : %s\n", ch$euler_ok))
  cat(sprintf("  incidence_ok : %s\n", ch$incidence_ok))
  cat(sprintf("  flatness_ok  : %s\n", ch$flatness_ok))
  
  cat("\nResult:\n")
  cat(sprintf("  feasible_symbolic = %s\n", x$feasible_symbolic))
  invisible(x)
}

# =============================================================================
# Block D — Enumeration module: counts + solutions + wizard (+ optional printer)
# =============================================================================

#' Count number of face-type solutions for fixed V,S
#' counts solutions to sum_{k>=4} (k-3)Q_k = S with P3 >= 0, k<=k_max
Pk_count_one <- function(V, S, k_max = V-1L) {
  V <- as.integer(V); S <- as.integer(S); k_max <- as.integer(k_max)
  if (V < 4L) stop("V must be >= 4.")
  if (S < 0L) return(0L)
  if (k_max < 3L) stop("k_max must be >= 3.")
  if (k_max > V-1L) k_max <- V-1L
  
  vr <- V_report(V)
  if (S > vr$bounds$S_max) return(0L)
  F <- (2L*V - 4L) - S
  
  w <- 1L:(k_max - 3L)
  if (length(w) == 0L) return(if (S == 0L) 1L else 0L)
  
  dp <- matrix(0, nrow = S + 1L, ncol = F + 1L)
  dp[1L, 1L] <- 1L  # sum=0, m=0
  
  for (wi in w) {
    dp_prev <- dp
    dp_new  <- dp_prev
    for (s0 in 0:S) {
      for (m0 in 0:F) {
        cur <- dp_prev[s0 + 1L, m0 + 1L]
        if (cur == 0L) next
        max_t <- min((S - s0) %/% wi, F - m0)
        if (max_t <= 0L) next
        for (t in 1:max_t) {
          dp_new[s0 + t*wi + 1L, m0 + t + 1L] <- dp_new[s0 + t*wi + 1L, m0 + t + 1L] + cur
        }
      }
    }
    dp <- dp_new
  }
  
  sum(dp[S + 1L, 1L:(F + 1L)])
}


#' Enumerate face-type solutions Pk for fixed V,S
#' Returns a data.frame with columns P3..Pkmax plus metadata.
Pk_solutions <- function(V, S, k_max = V-1L, max_solutions = 50000L,
                         require_kmax = FALSE) {
  V <- as.integer(V); S <- as.integer(S); k_max <- as.integer(k_max)
  if (V < 4L) stop("V must be >= 4.")
  if (k_max > V-1L) k_max <- V-1L
  if (k_max < 3L) stop("k_max must be >= 3.")
  
  vr <- V_report(V)
  if (S < 0L || S > vr$bounds$S_max) return(data.frame())
  
  E <- (3L*V - 6L) - S
  F <- (2L*V - 4L) - S
  
  # weights correspond to k=4..k_max: w = k-3
  w <- 1L:(k_max - 3L)
  K <- length(w)
  
  # only triangles allowed
  if (K == 0L) {
    if (S != 0L) return(data.frame())
    return(data.frame(V=V, S=S, E=E, F=F, P3=F))
  }
  
  sols <- vector("list", 0L)
  Q <- integer(K)
  
  rec <- function(i, S_rem, m_used) {
    if (length(sols) >= max_solutions) return(invisible(NULL))
    
    if (i > K) {
      if (S_rem != 0L) return(invisible(NULL))
      P3 <- F - m_used
      if (P3 < 0L) return(invisible(NULL))
      if (require_kmax && Q[K] == 0L) return(invisible(NULL))
      sols[[length(sols) + 1L]] <<- c(P3, Q)
      return(invisible(NULL))
    }
    
    wi <- w[i]
    max_t <- min(S_rem %/% wi, F - m_used)
    for (t in 0:max_t) {
      Q[i] <<- t
      rec(i + 1L, S_rem - t*wi, m_used + t)
      if (length(sols) >= max_solutions) break
    }
  }
  
  rec(1L, S, 0L)
  if (length(sols) == 0L) return(data.frame())
  
  mat <- do.call(rbind, sols)
  colnames(mat) <- c("P3", paste0("P", 4:k_max))
  
  out <- data.frame(V=V, S=S, E=E, F=F, mat, check.names = FALSE)
  rownames(out) <- NULL
  out
}

#' Wizard: pick V,S, count solutions, optionally enumerate, optionally print footprint
Pk_wizard <- function(V = NULL, S = NULL,
                      k_max = NULL,
                      auto_clamp = FALSE,
                      max_solutions = 50000L,
                      counts_only = FALSE,
                      show_solutions = TRUE,
                      show_footprint = FALSE) {
  
  if (is.null(V)) V <- 8L
  V <- as.integer(V)
  
  vr <- V_report(V)
  S_lo <- 0L
  S_hi <- vr$bounds$S_max
  
  if (is.null(k_max)) k_max <- V - 1L
  k_max <- as.integer(k_max)
  if (k_max > V - 1L) k_max <- V - 1L
  if (k_max < 3L) stop("k_max must be >= 3.")
  
  if (is.null(S)) {
    S <- as.integer(floor(S_hi / 3))
  } else {
    S <- as.integer(S)
  }
  
  if (S < S_lo || S > S_hi) {
    cat(sprintf("S not symbolically allowable for V=%d.\n", V))
    cat(sprintf("Allowed S range: [%d, %d]\n", S_lo, S_hi))
    S_suggest <- min(max(S, S_lo), S_hi)
    cat(sprintf("Your S=%d; suggested S=%d\n", S, S_suggest))
    if (!auto_clamp) {
      out <- list(V=V, S=S, S_suggest=S_suggest, V_report=vr, n=NA_integer_, solutions=NULL)
      class(out) <- "Pk_wizard"
      return(invisible(out))
    }
    S <- S_suggest
  }
  
  n <- Pk_count_one(V, S, k_max = k_max)
  cat(sprintf("Pk_wizard: V=%d, S=%d, k_max=%d => n_face_type_solutions=%s\n",
              V, S, k_max, n))
  
  if (!counts_only && !is.na(n) && n > max_solutions) {
    cat("Too many solutions to enumerate safely.\n")
    cat("Try smaller S, or set a smaller k_max, or use counts_only=TRUE.\n")
    out <- list(V=V, S=S, k_max=k_max, n=n, solutions=NULL, counts_only=TRUE)
    class(out) <- "Pk_wizard"
    return(invisible(out))
  }
  
  if (counts_only) {
    out <- list(V=V, S=S, k_max=k_max, n=n, solutions=NULL, counts_only=TRUE)
    class(out) <- "Pk_wizard"
    return(invisible(out))
  }
  
  sols <- Pk_solutions(V, S, k_max = k_max, max_solutions = max_solutions)
  
  poly_cols <- c("P3", paste0("P", 4:k_max))
  
  if (isTRUE(show_solutions)) {
    cat("\nFace-type solution matrix (rows = feasible Pk vectors):\n")
    print(sols[, poly_cols, drop = FALSE], row.names = FALSE)
  }
  
  if (isTRUE(show_footprint) && nrow(sols) > 0) {
    footprint <- as.data.frame((sols[, poly_cols, drop = FALSE] > 0L) * 1L)
    colnames(footprint) <- gsub("^P", "has_", colnames(footprint))
    cat("\nFootprint matrix (1 = present):\n")
    print(footprint, row.names = FALSE)
  }
  
  out <- list(V=V, S=S, k_max=k_max, n=n, solutions=sols, counts_only=FALSE)
  class(out) <- "Pk_wizard"
  invisible(out)
}

# Optional printer for returned object
print.Pk_wizard <- function(x, ...) {
  cat(sprintf("Pk_wizard(V=%d, S=%d, k_max=%d)\n", x$V, x$S, x$k_max))
  cat(sprintf("n_face_type_solutions=%s\n", x$n))
  if (!is.null(x$solutions) && nrow(x$solutions) > 0) {
    poly_cols <- c("P3", paste0("P", 4:x$k_max))
    print(x$solutions[, poly_cols, drop = FALSE], row.names = FALSE)
  }
  invisible(x)
}

# =============================================================================
# Block E — Solution-space summary (prints insight, not huge matrices)
# =============================================================================

#' Summarize face-type solution space from Pk_wizard() output or a solutions data.frame
#'
#' @param x A Pk_wizard object (from Pk_wizard()) or a data.frame of solutions
#' @param top_n Integer. Show top-N most common "footprints" (presence patterns).
#' @param digits Integer. Rounding digits for means.
#' @param quiet Logical. If TRUE, return list without printing.
#' @return A list with per-face-type stats, footprint frequencies, and entropy stats.
#' @export
Pk_summary <- function(x, top_n = 10L, digits = 3L, quiet = FALSE) {
  top_n <- as.integer(top_n); digits <- as.integer(digits)
  
  # ---- accept either Pk_wizard object or a raw solutions df ----
  if (inherits(x, "Pk_wizard")) {
    sols <- x$solutions
    meta <- list(V = x$V, S = x$S, k_max = x$k_max, n = x$n)
  } else if (is.data.frame(x)) {
    sols <- x
    meta <- list(V = if ("V" %in% names(sols)) unique(sols$V) else NA,
                 S = if ("S" %in% names(sols)) unique(sols$S) else NA,
                 k_max = NA, n = nrow(sols))
  } else {
    stop("Pk_summary() expects a Pk_wizard object or a solutions data.frame.")
  }
  
  if (is.null(sols) || nrow(sols) == 0L) {
    out <- list(meta = meta, n_solutions = 0L, face_stats = data.frame(),
                footprints = data.frame(), entropy = data.frame())
    if (!quiet) cat("Pk_summary: no solutions to summarize.\n")
    return(out)
  }
  
  # ---- detect Pk columns ----
  pk_cols <- grep("^P[0-9]+$", names(sols), value = TRUE)
  if (length(pk_cols) == 0L) stop("No Pk columns found (expected columns like P3, P4, ...).")
  
  # sort by degree
  deg <- as.integer(sub("^P", "", pk_cols))
  pk_cols <- pk_cols[order(deg)]
  deg <- sort(deg)
  
  X <- as.matrix(sols[, pk_cols, drop = FALSE])
  storage.mode(X) <- "integer"
  
  n_sol <- nrow(X)
  
  # ---- per-face stats ----
  present <- colSums(X > 0L) / n_sol
  face_stats <- data.frame(
    k = deg,
    col = pk_cols,
    present_pct = round(100 * present, digits),
    min = apply(X, 2, min),
    max = apply(X, 2, max),
    mean = round(colMeans(X), digits),
    stringsAsFactors = FALSE
  )
  
  # ---- footprint classes: presence/absence pattern ----
  fp <- (X > 0L) * 1L
  fp_key <- apply(fp, 1, paste0, collapse = "")
  tab <- sort(table(fp_key), decreasing = TRUE)
  top <- head(tab, top_n)
  
  # decode footprint key to "has_Pk" list
  decode_key <- function(key) {
    bits <- as.integer(strsplit(key, "", fixed = TRUE)[[1]])
    cols_present <- pk_cols[which(bits == 1L)]
    if (length(cols_present) == 0L) return("(none)")
    paste(cols_present, collapse = "+")
  }
  
  footprints <- data.frame(
    footprint_key = names(top),
    footprint = vapply(names(top), decode_key, character(1)),
    count = as.integer(top),
    pct = round(100 * as.integer(top) / n_sol, digits),
    stringsAsFactors = FALSE
  )
  
  # ---- entropy of face-type mix per solution (normalized by F) ----
  # H = -sum p_k log p_k with p_k = Pk/F, ignoring zeros.
  F_tot <- rowSums(X)
  safe_entropy <- function(row, F) {
    if (F <= 0) return(NA_real_)
    p <- row / F
    p <- p[p > 0]
    -sum(p * log(p))
  }
  H <- mapply(safe_entropy, split(X, row(X)), F_tot)  # vectorized-ish trick
  H <- as.numeric(H)
  
  entropy <- data.frame(
    H_min = min(H, na.rm = TRUE),
    H_max = max(H, na.rm = TRUE),
    H_mean = mean(H, na.rm = TRUE),
    H_median = stats::median(H, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  entropy[] <- lapply(entropy, function(z) round(z, digits))
  
  out <- list(
    meta = meta,
    n_solutions = n_sol,
    face_stats = face_stats,
    footprints = footprints,
    entropy = entropy
  )
  
  if (!quiet) {
    cat(sprintf("Pk_summary: V=%s, S=%s\n", meta$V, meta$S))
    if (!is.na(meta$k_max)) cat(sprintf("  k_max=%s\n", meta$k_max))
    cat(sprintf("  solutions=%d\n\n", n_sol))
    
    cat("Per-face-type stats:\n")
    print(face_stats[, c("k", "present_pct", "min", "max", "mean")], row.names = FALSE)
    
    cat("\nTop footprint classes (presence patterns):\n")
    print(footprints, row.names = FALSE)
    
    cat("\nEntropy of face-type mix (H):\n")
    print(entropy, row.names = FALSE)
    
    cat("\nInterpretation notes:\n")
    cat("• All statistics are computed across ALL admissible face-type solutions,\n")
    cat("  not for a single polyhedron.\n\n")
    cat("• Solutions include FALSE POSITIVES.\n\n")
    
    cat("• present_pct: percentage of solutions in which at least one k-gon appears.\n")
    cat("• min / max  : symbolic extremal counts of each face degree.\n")
    cat("• mean       : average number of k-gons per admissible solution.\n\n")
    
    cat("• Footprints record which face degrees are present (ignoring multiplicity).\n")
    cat("  They describe structural patterns and co-occurrence constraints.\n\n")
    
    cat("• Entropy H measures how mixed the face types are within a solution:\n")
    cat("    H = 0     → all faces have the same degree (single face type)\n")
    cat("    larger H → greater heterogeneity of face types\n")
    cat("  H_min / H_max bound symbolic flexibility,\n")
    cat("  while H_mean / H_median describe a typical admissible surface.\n")
    
  }
  
 invisible(out)
}

# =============================================================================
# Block F — Flatness-only loose upper bound (unrestricted partitions)
# =============================================================================

#' Loose upper bound on number of face-type multisets from flatness alone
#'
#' This computes a *very loose* upper bound on the number of distinct face-type
#' multisets compatible with a vertex count V, using unrestricted integer
#' partitions of the flatness budget S.
#'
#' Specifically, it returns sum_{S=0}^{S_max(V)} p(S), where p(S) is the
#' unrestricted partition function.
#'
#' Important: This ignores face-count constraints (F = 2V-4-S), ignores any
#' bound on face degrees (k_max), ignores vertex-type incidence constraints,
#' and ignores geometric realizability. It is intended only as a
#' "combinatorial explosion / search-space scale" diagnostic.
#'
#' @param V integer >= 4
#' @return A list with V, S_max, a by_S table (S, p_S, cumulative), and total.
#' @export
face_multiset_upperbound <- function(V) {
  V <- as.integer(V)
  if (length(V) != 1 || is.na(V) || V < 4L) stop("V must be a single integer >= 4.")
  
  vr <- V_report(V)
  S_max <- as.integer(vr$bounds$S_max)
  
  p <- partition_numbers_up_to(S_max)
  
  by_S <- data.frame(
    V = V,
    S = 0:S_max,
    p_S = p,
    cumulative = cumsum(p)
  )
  
  out <- list(
    V = V,
    S_max = S_max,
    by_S = by_S,
    total = sum(p)
  )
  class(out) <- "face_multiset_upperbound"
  out
}

#' Compute unrestricted partition numbers p(0..n) by dynamic programming
#'
#' p(0) = 1 and p(n) counts the number of integer partitions of n.
#' Uses an O(n^2) DP; fine for the S_max ranges used in polyenclose.
#'
#' @param n integer >= 0
#' @return integer vector length (n+1) giving p(0),...,p(n)
partition_numbers_up_to <- function(n) {
  n <- as.integer(n)
  if (length(n) != 1 || is.na(n) || n < 0L) stop("n must be a single integer >= 0.")
  
  # Use numeric to avoid integer overflow for moderate n; users can format later.
  p <- numeric(n + 1L)
  p[1L] <- 1  # p(0)=1
  
  # Standard partition DP: for each part size k, update counts
  for (k in 1L:n) {
    for (i in k:n) {
      p[i + 1L] <- p[i + 1L] + p[i - k + 1L]
    }
  }
  
  p
}

#' @export
print.face_multiset_upperbound <- function(x, ...) {
  cat(sprintf("face_multiset_upperbound(V=%d)\n", x$V))
  cat(sprintf("S_max(V) = %d\n", x$S_max))
  cat(sprintf("sum_{S=0}^{S_max} p(S) = %s\n\n", format(x$total, scientific = FALSE)))
  
  cat("By-S table (unrestricted partitions):\n")
  print(x$by_S, row.names = FALSE)
  
  cat("\nInterpretive note:\n")
  cat("• At S = 0 (fully triangulated surface), the external face-type diversity is minimal.\n")
  cat("• The maximal tetrahedral endpoint in the SALT+MIE ladder (T = T_max) occurs only in bipyramidal configurations,\n")
  cat("  and is attainable only at S = 0; for any S > 0 this extremal case is excluded.\n")
  cat("  is attainable only at S = 0 as well as the rest of the T/Ni/Si ladder; for any S > 0 this extremal case is excluded.\n")
  cat("• As S increases, face-type multiset diversity grows rapidly, while the\n")
  cat("  maximal achievable tetrahedral count drops from its S=0 extremum to minimal range.\n")
  
  invisible(x)
}
