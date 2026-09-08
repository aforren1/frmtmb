# Route 2: the two-dimensional Fokker-Planck solve for three alternatives.
#
# This is what gddm() would have to become to reach n alternatives in
# general: a grid solve on an (n-1)-dimensional domain with absorbing
# walls. Written for n = 3, where the domain is a triangle, so that the
# closed form of images.R is available as a reference.
#
# COORDINATES. Work in s_k = <v_k, z>, the evidence for alternative k
# relative to the mean, and keep s1 and s2 (s3 = -s1-s2). The continuation
# region is then {s1 < c, s2 < c, s1+s2 > -c}: a RIGHT triangle, two of
# whose sides are grid lines and whose third is a grid diagonal. No side is
# staircased, which is what keeps the wall flux second order. The price is a
# noise covariance with a cross term,
#
#     Sigma = [[2/3, -1/3], [-1/3, 2/3]],
#
# because the s_k are not orthogonal. The alternative, whitening to make the
# noise isotropic, turns the triangle equilateral and puts two of its sides
# off the grid, which is worse.
#
# SCHEME. Douglas ADI: one explicit stage carrying the whole operator
# including the cross term, then two implicit corrections, one per
# direction, each a set of independent tridiagonal systems. Second order in
# space; the explicit cross term makes it first order in the step, which is
# a real cost against the shipped 1-D Crank-Nicolson and is priced rather
# than hidden.
#
# PADDING. The lines are ragged: the domain is half a square. Rather than
# sweep ragged lines, the solve runs over the FULL (N-1) by (N-1) square and
# the rows outside the triangle carry the identity, (lo, di, up) = (0, 1, 0),
# with a zero right-hand side. The Thomas sweep then leaves them at zero and
# they decouple exactly, so every line has the same length and the whole
# recurrence vectorizes across lines. That is what keeps the tape build from
# being one R-level scalar iteration per node per step. It solves twice as
# many nodes as the triangle has, which is counted in the cost.

# Grid, masks and index tables. All data: nothing here sees a parameter.
gd_pde_grid <- function(N) {
  n <- N - 1L                       # interior candidates per axis
  ij <- expand.grid(i = seq_len(n), j = seq_len(n))
  inside <- (ij$i + ij$j) >= (N + 1L)
  msk <- as.numeric(inside)
  # flat index into the padded (N+1)^2 array, i and j from 0 to N
  pad <- function(i, j) i + j * (N + 1L) + 1L
  list(N = N, n = n, np = N + 1L, msk = msk,
       ij = ij, nin = sum(inside),
       ctr = pad(ij$i, ij$j),
       ip = pad(ij$i + 1L, ij$j), im = pad(ij$i - 1L, ij$j),
       jp = pad(ij$i, ij$j + 1L), jm = pad(ij$i, ij$j - 1L),
       pp = pad(ij$i + 1L, ij$j + 1L), pm = pad(ij$i + 1L, ij$j - 1L),
       mp = pad(ij$i - 1L, ij$j + 1L), mm = pad(ij$i - 1L, ij$j - 1L),
       # lines: element k of line l lives at these flat positions in the
       # length n^2 interior vector, i fastest
       li = lapply(seq_len(n), function(k) (seq_len(n) - 1L) * n + k),
       lj = lapply(seq_len(n), function(k) (k - 1L) * n + seq_len(n)),
       # flux layers
       f1a = which(ij$i == n - 0L & inside), f1b = which(ij$i == n - 1L & inside),
       f2a = which(ij$j == n - 0L & inside), f2b = which(ij$j == n - 1L & inside),
       f3a = which(ij$i + ij$j == N + 1L),
       # the inward neighbour of a diagonal-layer-1 node is not a grid node,
       # so pair it with the mean of the two that straddle it
       # one straddling neighbour at each end of the diagonal sits on an
       # axis wall, where p is zero; index 1 with weight 0 stands in for it,
       # because an NA index would poison the whole AD vector
       f3b1 = local({ m <- match(paste(ij$i[ij$i + ij$j == N + 1L] + 1L,
                                       ij$j[ij$i + ij$j == N + 1L]),
                                 paste(ij$i, ij$j))
                      list(i = ifelse(is.na(m), 1L, m),
                           w = as.numeric(!is.na(m))) }),
       f3b2 = local({ m <- match(paste(ij$i[ij$i + ij$j == N + 1L],
                                       ij$j[ij$i + ij$j == N + 1L] + 1L),
                                 paste(ij$i, ij$j))
                      list(i = ifelse(is.na(m), 1L, m),
                           w = as.numeric(!is.na(m))) }),
       # GHOST. The cross-derivative stencil at a node one layer inside the
       # DIAGONAL wall reads (i-1, j-1), which is a layer OUTSIDE the
       # domain. Masking it to zero is wrong: Dirichlet pins p on the wall,
       # not beyond it, and the smooth continuation is ODD, so the ghost is
       # minus the value at the mirror node. Reflecting (i-1, j-1) across
       # i + j = N gives (N-j+1, N-i+1), which for i + j = N + 1 is the node
       # itself. So the correction is exactly -p at the node. Left as a zero
       # the diagonal wall's flux comes out 24 percent low and does NOT
       # converge under refinement, while the two grid-aligned walls do:
       # measured in probe3. The two axis walls have no ghost, because every
       # interior node has i and j at least 2 and i, j at most N-1.
       ghost = which(ij$i + ij$j == N + 1L))
}

# Vectorized Thomas across all n lines at once. lo/di/up/rhs are length n^2
# with lo[p] = A[k,k-1] and up[p] = A[k,k+1] for the node at position p.
gd_pde_thomas <- function(lo, di, up, rhs, lines) {
  "[<-" <- RTMB::ADoverload("[<-")
  m <- length(lines)
  z <- rhs * 0
  cp <- z; dp <- z; out <- z
  i1 <- lines[[1L]]
  cp[i1] <- up[i1] / di[i1]
  dp[i1] <- rhs[i1] / di[i1]
  for (k in 2:m) {
    ik <- lines[[k]]; ip <- lines[[k - 1L]]
    q <- di[ik] - lo[ik] * cp[ip]
    cp[ik] <- up[ik] / q
    dp[ik] <- (rhs[ik] - lo[ik] * dp[ip]) / q
  }
  ie <- lines[[m]]
  out[ie] <- dp[ie]
  for (k in (m - 1L):1L) {
    ik <- lines[[k]]
    out[ik] <- dp[ik] - cp[ik] * out[lines[[k + 1L]]]
  }
  out
}

# One condition's three defective densities on the time grid.
#   b     drift in the s coordinates, b_k = mu_k - mean(mu), length 2
#   cc    threshold (evidence over the mean needed to stop)
#   s0    start in s coordinates, length 2 (the multi-alternative bias)
gd_pde_solve <- function(b, cc, s0, g, nt, dt, theta = 0.5,
                         ghost = TRUE) {
  "[<-" <- RTMB::ADoverload("[<-")
  N <- g$N; n <- g$n
  h <- 3 * cc / N
  s11 <- 2 / 3; s22 <- 2 / 3; s12 <- -1 / 3

  # start: bilinear spreading of a point mass over the four nearest nodes,
  # written with the same smooth-shift idea gddm() uses in time, so that no
  # grid index depends on a parameter. The node positions are
  # s = -2c + i h, so the fractional index is (s0 + 2c)/h, and c is a
  # parameter, so the WEIGHT moves and the index does not: spread over
  # every node with a B-spline in the fractional offset.
  ii <- g$ij$i; jj <- g$ij$j
  fi <- (s0[1L] + 2 * cc) / h
  fj <- (s0[2L] + 2 * cc) / h
  rho <- frmtmb.eam:::gd_b3(fi - ii) * frmtmb.eam:::gd_b3(fj - jj) * g$msk
  rho <- rho / (sum(rho) * h * h)

  # constant tridiagonal coefficients, masked to the identity outside
  c1lo <- -theta * dt * (b[1L] / (2 * h) + s11 / (2 * h * h))
  c1up <- -theta * dt * (-b[1L] / (2 * h) + s11 / (2 * h * h))
  c1di <- theta * dt * s11 / (h * h)
  c2lo <- -theta * dt * (b[2L] / (2 * h) + s22 / (2 * h * h))
  c2up <- -theta * dt * (-b[2L] / (2 * h) + s22 / (2 * h * h))
  c2di <- theta * dt * s22 / (h * h)
  msk <- g$msk
  lo1 <- msk * c1lo; up1 <- msk * c1up; di1 <- 1 + msk * c1di
  lo2 <- msk * c2lo; up2 <- msk * c2up; di2 <- 1 + msk * c2di
  # no neighbour beyond the ends of a line
  lo1[g$li[[1L]]] <- lo1[g$li[[1L]]] * 0
  up1[g$li[[n]]] <- up1[g$li[[n]]] * 0
  lo2[g$lj[[1L]]] <- lo2[g$lj[[1L]]] * 0
  up2[g$lj[[n]]] <- up2[g$lj[[n]]] * 0

  P <- numeric(g$np * g$np) + rho[1L] * 0
  pad <- function(v) { P[g$ctr] <- v; P }
  A1 <- function(Q) msk * (-b[1L] * (Q[g$ip] - Q[g$im]) / (2 * h) +
                           s11 / 2 * (Q[g$ip] - 2 * Q[g$ctr] + Q[g$im]) / (h * h))
  A2 <- function(Q) msk * (-b[2L] * (Q[g$jp] - Q[g$jm]) / (2 * h) +
                           s22 / 2 * (Q[g$jp] - 2 * Q[g$ctr] + Q[g$jm]) / (h * h))
  A0 <- function(Q, pv) {
    v <- msk * s12 * (Q[g$pp] - Q[g$pm] - Q[g$mp] + Q[g$mm]) / (4 * h * h)
    if (ghost) v[g$ghost] <- v[g$ghost] - s12 * pv[g$ghost] / (4 * h * h)
    v
  }

  fl <- vector("list", 3L)
  for (k in 1:3) fl[[k]] <- numeric(nt + 1L) + rho[1L] * 0
  p <- rho
  # The wall flux, read in the s coordinates, where the divergence theorem
  # is the ordinary Euclidean one. On the wall p = 0 and the tangential
  # derivative with it, so the flux is -(1/2) nu' Sigma nu times the normal
  # derivative. For the s1 wall that is -(1/2)(2/3) d_1 p integrated over
  # ds2, and for the diagonal wall nu' Sigma nu is also 2/3 while the normal
  # spacing h/sqrt(2) and the node spacing h sqrt(2) cancel, so both walls
  # come out as 1/6 times a one-sided second-order difference over 2h times
  # the node spacing h. Hence one constant, 1/6.
  fc <- 1 / 6
  for (st in seq_len(nt)) {
    Q <- pad(p)
    a1 <- A1(Q); a2 <- A2(Q); a0 <- A0(Q, p)
    y0 <- p + dt * (a0 + a1 + a2)
    y1 <- gd_pde_thomas(lo1, di1, up1, y0 - theta * dt * a1, g$li)
    y2 <- gd_pde_thomas(lo2, di2, up2, y1 - theta * dt * a2, g$lj)
    p <- y2 * msk
    fl[[1L]][st + 1L] <- fc * (4 * sum(p[g$f1a]) - sum(p[g$f1b]))
    fl[[2L]][st + 1L] <- fc * (4 * sum(p[g$f2a]) - sum(p[g$f2b]))
    fl[[3L]][st + 1L] <- fc * (4 * sum(p[g$f3a]) -
                               (sum(p[g$f3b1$i] * g$f3b1$w) +
                                sum(p[g$f3b2$i] * g$f3b2$w)) / 2)
  }
  fl
}
