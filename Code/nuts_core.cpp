// nuts_core.cpp
// ---------------------------------------------------------------------------
// Compiled backend for a No-U-Turn Sampler (NUTS).
//
// Implements the modern NUTS trajectory machinery in C++ (Rcpp + RcppEigen):
//   * gradient-threaded leapfrog (one gradient evaluation reused across steps),
//   * recursive tree doubling with MULTINOMIAL sampling via stable log-sum-exp
//     and BIASED PROGRESSIVE selection favouring the newer subtree,
//   * generalized (metric-weighted) U-turn checks on both trajectory ends,
//   * divergence detection.
//
// References:
//   Hoffman & Gelman (2014), "The No-U-Turn Sampler", JMLR 15:1593-1623.
//   Betancourt (2017), "A Conceptual Introduction to Hamiltonian Monte Carlo",
//     arXiv:1701.02434  (multinomial sampling, generalized U-turn, biased
//     progressive sampling; matches Stan's current implementation).
//
// Two entry points share all tree logic through a small gradient/log-density
// abstraction:
//   * nuts_step_callback : target supplied as R functions f, grad_f
//                          (works for ANY target; incurs R-callback overhead)
//   * nuts_step_xptr     : target supplied as compiled XPtr functions
//                          (zero R-callback overhead)
// ---------------------------------------------------------------------------

// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>
#include <functional>

using namespace Rcpp;
using Eigen::VectorXd;

// Signatures for compiled targets passed via XPtr.
typedef double     (*logdens_fn)(const VectorXd& theta);
typedef VectorXd   (*grad_fn)   (const VectorXd& theta);

// ---------------------------------------------------------------------------
// Target abstraction. A Target bundles a log-density and its gradient behind
// std::function so the tree code is written once and reused by both the
// R-callback tier and the compiled XPtr tier. We also count evaluations for
// diagnostics.
// ---------------------------------------------------------------------------
struct Target {
  std::function<double(const VectorXd&)>   logdens;
  std::function<VectorXd(const VectorXd&)> grad;
  // Optional FUSED evaluator: returns log-density AND fills the gradient in a
  // single call. When present (has_fused), the leapfrog/base-case hot path uses
  // it so each leaf costs ONE target evaluation instead of separate f and g
  // calls at the same theta. logdens/grad remain populated (derived from the
  // same fused function) for the rare code paths that need only one of them
  // (find_reasonable_epsilon, the once-per-iteration initial state).
  std::function<double(const VectorXd&, VectorXd&)> fused;
  bool has_fused = false;
  long n_logdens = 0;
  long n_grad    = 0;
  long n_fused   = 0;

  inline double f(const VectorXd& theta) { ++n_logdens; return logdens(theta); }
  inline VectorXd g(const VectorXd& theta) { ++n_grad; return grad(theta); }
  // Fused: returns logp, writes gradient into grad_out. One R callback.
  inline double fg(const VectorXd& theta, VectorXd& grad_out) {
    ++n_fused; return fused(theta, grad_out);
  }
};

// ---------------------------------------------------------------------------
// Leapfrog integrator, gradient-threaded.
//
// The incoming `grad` is grad_f(theta). We return the updated (theta, r) AND
// the gradient at the new theta, so the caller can feed it straight back in as
// the next step's incoming gradient. Over a trajectory this uses one gradient
// evaluation per step rather than two.
//
// inv_M is the diagonal of the inverse mass matrix (== 1 / M_diag), i.e. the
// momentum covariance we sample from is M = diag(1/inv_M).
// ---------------------------------------------------------------------------
// Returns the log-density at the new theta. When the target is NOT fused this
// is NaN (the caller then computes f() itself, preserving the original call
// pattern exactly); when it IS fused the gradient and log-density come from a
// single evaluation and logp is returned so the base case need not call f again.
static inline double leapfrog(Target& tgt,
                              VectorXd& theta, VectorXd& r, VectorXd& grad,
                              double eps, const VectorXd& inv_M) {
  r.noalias() += 0.5 * eps * grad;            // half momentum kick (reuses grad)
  theta.noalias() += eps * (inv_M.array() * r.array()).matrix();  // full drift
  double logp = std::numeric_limits<double>::quiet_NaN();
  if (tgt.has_fused) {
    logp = tgt.fg(theta, grad);               // one fused eval: grad + logp
  } else {
    grad = tgt.g(theta);                      // single new gradient
  }
  r.noalias() += 0.5 * eps * grad;            // half momentum kick
  return logp;
}

// Joint log-density (negative Hamiltonian): log pi(theta) - 0.5 r' M^{-1} r.
static inline double joint_logdens(Target& tgt, const VectorXd& theta,
                                   const VectorXd& r, const VectorXd& inv_M) {
  return tgt.f(theta) - 0.5 * (inv_M.array() * r.array().square()).sum();
}
// As above but with a precomputed log pi(theta) (avoids a redundant f call).
static inline double joint_from_logp(double logp, const VectorXd& r,
                                     const VectorXd& inv_M) {
  return logp - 0.5 * (inv_M.array() * r.array().square()).sum();
}

// ---------------------------------------------------------------------------
// Generalized U-turn criterion (Betancourt 2017). The trajectory keeps
// expanding while the momentum at both ends still points "outward" along the
// span (theta_plus - theta_minus), measured in the metric via M^{-1} r.
// ---------------------------------------------------------------------------
static inline bool no_u_turn(const VectorXd& theta_plus, const VectorXd& theta_minus,
                             const VectorXd& r_plus, const VectorXd& r_minus,
                             const VectorXd& inv_M) {
  VectorXd delta = theta_plus - theta_minus;
  double a = delta.dot((inv_M.array() * r_minus.array()).matrix());
  double b = delta.dot((inv_M.array() * r_plus.array()).matrix());
  return (a >= 0.0) && (b >= 0.0);
}

// ---------------------------------------------------------------------------
// Recursive subtree state.
// ---------------------------------------------------------------------------
struct Tree {
  VectorXd theta_minus, r_minus, grad_minus;   // leftmost leaf
  VectorXd theta_plus,  r_plus,  grad_plus;    // rightmost leaf
  VectorXd theta_prop;                         // selected proposal from subtree
  double   grad_prop_logp;                     // log pi at the proposal (cached)
  VectorXd grad_prop;                          // grad at the proposal (cached)
  double   log_w;        // log of summed multinomial weight over valid leaves
  bool     s;           // 1 = keep going (no U-turn AND no divergence)
  bool     divergent;   // TRUE if a leaf in this subtree diverged (energy blew up)
  double   sum_alpha;   // sum of Metropolis accept probs (for dual averaging)
  int      n_alpha;     // number of leaves contributing to sum_alpha
};

// RNG helper: uniform(0,1) drawn from R's stream so set.seed() is respected.
static inline double runif01() { return ::unif_rand(); }

// ---------------------------------------------------------------------------
// build_tree: recursive doubling. `v` in {-1,+1} is the direction of
// integration in fictitious time.  `log_u0` is the joint log-density of the
// initial (theta0, r0) state for this iteration -- passed in and thus computed
// only ONCE per iteration (not per leaf).
// ---------------------------------------------------------------------------
static Tree build_tree(Target& tgt,
                       const VectorXd& theta, const VectorXd& r, const VectorXd& grad,
                       double log_u0, int v, int j, double eps,
                       const VectorXd& inv_M, double Delta_max) {
  Tree t;
  if (j == 0) {
    // Base case: one leapfrog step in direction v.
    VectorXd th = theta, rr = r, gr = grad;
    double logp = leapfrog(tgt, th, rr, gr, v * eps, inv_M);
    // Fused targets return logp from the leapfrog (same theta, one eval);
    // otherwise compute it now -- identical to the original call pattern.
    if (!tgt.has_fused) logp = tgt.f(th);
    double joint = joint_from_logp(logp, rr, inv_M);

    t.theta_minus = th; t.r_minus = rr; t.grad_minus = gr;
    t.theta_plus  = th; t.r_plus  = rr; t.grad_plus  = gr;
    t.theta_prop  = th; t.grad_prop = gr; t.grad_prop_logp = logp;

    double delta_energy = joint - log_u0;           // = -(H_new - H_0)
    if (!std::isfinite(joint)) delta_energy = -std::numeric_limits<double>::infinity();

    t.log_w = std::isfinite(joint) ? joint : -std::numeric_limits<double>::infinity();
    t.divergent = !(delta_energy > -Delta_max);     // energy blew up
    t.s = !t.divergent;                             // base case stops only on divergence
    // Metropolis accept prob of this leaf vs the initial state.
    double a = std::exp(std::min(0.0, delta_energy));
    if (!std::isfinite(a)) a = 0.0;
    t.sum_alpha = a;
    t.n_alpha = 1;
    return t;
  }

  // Recurse: build the first half-subtree.
  Tree t0 = build_tree(tgt, theta, r, grad, log_u0, v, j - 1, eps, inv_M, Delta_max);
  t = t0;
  if (!t0.s) { return t; }   // early stop: divergence/U-turn already inside (t.divergent carried)

  // Build the second half-subtree from the appropriate boundary.
  Tree t1;
  if (v == -1) {
    t1 = build_tree(tgt, t0.theta_minus, t0.r_minus, t0.grad_minus,
                    log_u0, v, j - 1, eps, inv_M, Delta_max);
    t.theta_minus = t1.theta_minus; t.r_minus = t1.r_minus; t.grad_minus = t1.grad_minus;
  } else {
    t1 = build_tree(tgt, t0.theta_plus, t0.r_plus, t0.grad_plus,
                    log_u0, v, j - 1, eps, inv_M, Delta_max);
    t.theta_plus = t1.theta_plus; t.r_plus = t1.r_plus; t.grad_plus = t1.grad_plus;
  }

  // Biased progressive multinomial selection (Betancourt 2017, eq. for
  // choosing between subtrees): accept t1's proposal with probability
  // w1 / (w0 + w1) computed stably in log-space.
  double log_sum;
  if (t0.log_w > t1.log_w)
    log_sum = t0.log_w + std::log1p(std::exp(t1.log_w - t0.log_w));
  else if (std::isfinite(t1.log_w))
    log_sum = t1.log_w + std::log1p(std::exp(t0.log_w - t1.log_w));
  else
    log_sum = t0.log_w;

  double p_choose_t1 = std::exp(t1.log_w - log_sum);   // in [0,1]
  if (std::isfinite(p_choose_t1) && runif01() < p_choose_t1) {
    t.theta_prop     = t1.theta_prop;
    t.grad_prop      = t1.grad_prop;
    t.grad_prop_logp = t1.grad_prop_logp;
  }
  t.log_w = log_sum;

  // Combine statistics.
  t.sum_alpha = t0.sum_alpha + t1.sum_alpha;
  t.n_alpha   = t0.n_alpha + t1.n_alpha;

  // Divergence is monotone: either half carrying one taints the combined tree.
  t.divergent = t0.divergent || t1.divergent;
  // Stop if the sub-subtree diverged/U-turned, or a U-turn now spans the ends.
  t.s = t1.s && no_u_turn(t.theta_plus, t.theta_minus, t.r_plus, t.r_minus, inv_M);
  return t;
}

// ---------------------------------------------------------------------------
// One full NUTS transition. Given the current theta and step size eps, samples
// a new theta and returns it along with diagnostics. `inv_M` = 1 / M_diag.
// The momentum is drawn r ~ N(0, M) = N(0, diag(1/inv_M)).
// ---------------------------------------------------------------------------
static List nuts_transition(Target& tgt, const VectorXd& theta0,
                            double eps, const VectorXd& inv_M,
                            int max_treedepth, double Delta_max) {
  const int d = theta0.size();

  // Sample momentum r ~ N(0, M). sd of component i is sqrt(M_ii) = 1/sqrt(inv_M_i).
  VectorXd r0(d);
  for (int i = 0; i < d; ++i) r0[i] = ::norm_rand() / std::sqrt(inv_M[i]);

  double logp0 = tgt.f(theta0);
  double log_u0 = joint_from_logp(logp0, r0, inv_M);   // computed ONCE per iter
  VectorXd grad0 = tgt.g(theta0);

  VectorXd theta_minus = theta0, theta_plus = theta0;
  VectorXd r_minus = r0, r_plus = r0;
  VectorXd grad_minus = grad0, grad_plus = grad0;

  VectorXd theta = theta0;
  double theta_logp = logp0;
  double log_w = log_u0;             // running log summed weight of whole tree
  bool s = true;
  int depth = 0;
  int n_leapfrog = 0;
  bool divergent = false;
  double sum_alpha = 0.0; int n_alpha = 0;

  while (s && depth < max_treedepth) {
    int v = (runif01() < 0.5) ? -1 : 1;
    Tree t;
    if (v == -1) {
      t = build_tree(tgt, theta_minus, r_minus, grad_minus,
                     log_u0, v, depth, eps, inv_M, Delta_max);
      theta_minus = t.theta_minus; r_minus = t.r_minus; grad_minus = t.grad_minus;
    } else {
      t = build_tree(tgt, theta_plus, r_plus, grad_plus,
                     log_u0, v, depth, eps, inv_M, Delta_max);
      theta_plus = t.theta_plus; r_plus = t.r_plus; grad_plus = t.grad_plus;
    }
    n_leapfrog += (1 << depth);

    // Progressive multinomial acceptance of the new subtree's proposal against
    // the existing tree: accept with prob min(1, w_new / w_old).
    if (t.s) {
      double log_accept = t.log_w - log_w;   // may be > 0
      if (std::log(runif01()) < log_accept) {
        theta = t.theta_prop;
        theta_logp = t.grad_prop_logp;
      }
    }

    sum_alpha += t.sum_alpha;
    n_alpha   += t.n_alpha;

    // Update total weight (log-sum-exp of old tree and new subtree).
    if (log_w > t.log_w)
      log_w = log_w + std::log1p(std::exp(t.log_w - log_w));
    else if (std::isfinite(t.log_w))
      log_w = t.log_w + std::log1p(std::exp(log_w - t.log_w));

    if (t.divergent) divergent = true;   // a leaf's energy blew up

    // Continue only if the new subtree is valid AND no U-turn spans the ends.
    s = t.s && no_u_turn(theta_plus, theta_minus, r_plus, r_minus, inv_M);
    depth += 1;
  }

  double accept_stat = (n_alpha > 0) ? (sum_alpha / n_alpha) : 0.0;
  // Report the potential energy -log pi(theta) at the accepted point. (A full
  // Hamiltonian-energy series for E-BFMI would need the momentum at the
  // accepted state, which multinomial NUTS does not retain; the potential
  // energy series is sufficient for the diagnostics we expose.)
  double energy = -theta_logp;

  return List::create(
    _["theta"] = theta,
    _["accept_stat"] = accept_stat,
    _["tree_depth"] = depth,
    _["n_leapfrog"] = n_leapfrog,
    _["divergent"] = divergent,
    _["energy"] = energy,
    _["logp"] = theta_logp
  );
}

// ===========================================================================
// R-callback tier
// ===========================================================================

static Target make_callback_target(Function f, Function grad_f) {
  Target tgt;
  tgt.logdens = [f](const VectorXd& theta) -> double {
    NumericVector th(theta.data(), theta.data() + theta.size());
    NumericVector out = f(th);
    return out[0];
  };
  tgt.grad = [grad_f](const VectorXd& theta) -> VectorXd {
    NumericVector th(theta.data(), theta.data() + theta.size());
    NumericVector g = grad_f(th);
    return Eigen::Map<VectorXd>(g.begin(), g.size());
  };
  return tgt;
}

// [[Rcpp::export]]
List nuts_step_callback(NumericVector theta0, Function f, Function grad_f,
                        double eps, NumericVector inv_M,
                        int max_treedepth = 10, double Delta_max = 1000.0) {
  Target tgt = make_callback_target(f, grad_f);
  VectorXd th = Eigen::Map<VectorXd>(theta0.begin(), theta0.size());
  VectorXd im = Eigen::Map<VectorXd>(inv_M.begin(), inv_M.size());
  return nuts_transition(tgt, th, eps, im, max_treedepth, Delta_max);
}

// Build a FUSED callback target from a single R function fg(theta) returning a
// list(logp=<scalar>, grad=<numeric vector>). The hot path (leapfrog + base
// case) uses one fused call per leaf; logdens/grad are also derived from fg so
// the rare paths that need just one work unchanged. This is numerically a
// drop-in for make_callback_target(f, grad_f) provided fg returns exactly the
// same logp/grad as f/grad_f -- so the NUTS trajectory and RNG-driven draws are
// bitwise identical, at ~half the per-leaf cost.
static Target make_fused_callback_target(Function fg) {
  Target tgt;
  auto call_fg = [fg](const VectorXd& theta, VectorXd* grad_out) -> double {
    NumericVector th(theta.data(), theta.data() + theta.size());
    List out = fg(th);
    if (grad_out) {
      NumericVector g = out["grad"];
      *grad_out = Eigen::Map<VectorXd>(g.begin(), g.size());
    }
    NumericVector lp = out["logp"];
    return lp[0];
  };
  tgt.fused = [call_fg](const VectorXd& theta, VectorXd& grad_out) -> double {
    return call_fg(theta, &grad_out);
  };
  tgt.has_fused = true;
  // Derived single-purpose accessors (used off the hot path). These still cost
  // one fused evaluation each, which is acceptable: the initial-state logp/grad
  // (once per iteration) and find_reasonable_epsilon (once per adaptation step).
  tgt.logdens = [call_fg](const VectorXd& theta) -> double {
    return call_fg(theta, nullptr);
  };
  tgt.grad = [call_fg](const VectorXd& theta) -> VectorXd {
    VectorXd g; call_fg(theta, &g); return g;
  };
  return tgt;
}

// Run a whole chain in C++ with an R-function target. Keeping the iteration
// loop in C++ (rather than calling nuts_step_callback from R n_iter times)
// avoids per-iteration R interpreter overhead; adaptation state is passed in
// as a fixed eps / inv_M (adaptation itself is orchestrated from R between
// calls to this sampler for the sampling phase).
// [[Rcpp::export]]
List nuts_sample_callback(NumericVector theta0, Function f, Function grad_f,
                          int n_iter, NumericVector eps_seq, NumericMatrix inv_M,
                          int max_treedepth = 10, double Delta_max = 1000.0) {
  Target tgt = make_callback_target(f, grad_f);
  const int d = theta0.size();
  VectorXd th = Eigen::Map<VectorXd>(theta0.begin(), theta0.size());

  NumericMatrix trace(n_iter, d);
  NumericVector accept(n_iter), depth(n_iter), nleap(n_iter), div(n_iter), energy(n_iter);

  bool inv_M_fixed = (inv_M.nrow() == 1);
  for (int it = 0; it < n_iter; ++it) {
    VectorXd im(d);
    int row = inv_M_fixed ? 0 : it;
    for (int i = 0; i < d; ++i) im[i] = inv_M(row, i);
    double eps = (eps_seq.size() == 1) ? eps_seq[0] : eps_seq[it];

    List step = nuts_transition(tgt, th, eps, im, max_treedepth, Delta_max);
    NumericVector nt = step["theta"];
    th = Eigen::Map<VectorXd>(nt.begin(), nt.size());
    for (int i = 0; i < d; ++i) trace(it, i) = th[i];
    accept[it] = step["accept_stat"];
    depth[it]  = step["tree_depth"];
    nleap[it]  = step["n_leapfrog"];
    div[it]    = (bool)step["divergent"] ? 1.0 : 0.0;
    energy[it] = step["energy"];
  }
  return List::create(_["theta"] = trace, _["accept_stat"] = accept,
                      _["tree_depth"] = depth, _["n_leapfrog"] = nleap,
                      _["divergent"] = div, _["energy"] = energy,
                      _["n_logdens"] = (double)tgt.n_logdens,
                      _["n_grad"] = (double)tgt.n_grad);
}

// Single NUTS transition with a FUSED R target (fg -> list(logp, grad)).
// [[Rcpp::export]]
List nuts_step_fused_callback(NumericVector theta0, Function fg,
                              double eps, NumericVector inv_M,
                              int max_treedepth = 10, double Delta_max = 1000.0) {
  Target tgt = make_fused_callback_target(fg);
  VectorXd th = Eigen::Map<VectorXd>(theta0.begin(), theta0.size());
  VectorXd im = Eigen::Map<VectorXd>(inv_M.begin(), inv_M.size());
  return nuts_transition(tgt, th, eps, im, max_treedepth, Delta_max);
}

// Whole-chain sampler with a FUSED R target. Mirrors nuts_sample_callback but
// each leaf makes ONE fused evaluation (logp+grad) instead of separate f and g.
// [[Rcpp::export]]
List nuts_sample_fused_callback(NumericVector theta0, Function fg,
                                int n_iter, NumericVector eps_seq,
                                NumericMatrix inv_M,
                                int max_treedepth = 10, double Delta_max = 1000.0) {
  Target tgt = make_fused_callback_target(fg);
  const int d = theta0.size();
  VectorXd th = Eigen::Map<VectorXd>(theta0.begin(), theta0.size());

  NumericMatrix trace(n_iter, d);
  NumericVector accept(n_iter), depth(n_iter), nleap(n_iter), div(n_iter), energy(n_iter);

  bool inv_M_fixed = (inv_M.nrow() == 1);
  for (int it = 0; it < n_iter; ++it) {
    VectorXd im(d);
    int row = inv_M_fixed ? 0 : it;
    for (int i = 0; i < d; ++i) im[i] = inv_M(row, i);
    double eps = (eps_seq.size() == 1) ? eps_seq[0] : eps_seq[it];

    List step = nuts_transition(tgt, th, eps, im, max_treedepth, Delta_max);
    NumericVector nt = step["theta"];
    th = Eigen::Map<VectorXd>(nt.begin(), nt.size());
    for (int i = 0; i < d; ++i) trace(it, i) = th[i];
    accept[it] = step["accept_stat"];
    depth[it]  = step["tree_depth"];
    nleap[it]  = step["n_leapfrog"];
    div[it]    = (bool)step["divergent"] ? 1.0 : 0.0;
    energy[it] = step["energy"];
  }
  return List::create(_["theta"] = trace, _["accept_stat"] = accept,
                      _["tree_depth"] = depth, _["n_leapfrog"] = nleap,
                      _["divergent"] = div, _["energy"] = energy,
                      _["n_fused"] = (double)tgt.n_fused);
}

// ===========================================================================
// XPtr (fully compiled) tier
// ===========================================================================

static Target make_xptr_target(SEXP f_ptr, SEXP g_ptr) {
  XPtr<logdens_fn> xf(f_ptr);
  XPtr<grad_fn>    xg(g_ptr);
  logdens_fn fp = *xf;
  grad_fn    gp = *xg;
  Target tgt;
  tgt.logdens = [fp](const VectorXd& theta) -> double { return fp(theta); };
  tgt.grad    = [gp](const VectorXd& theta) -> VectorXd { return gp(theta); };
  return tgt;
}

// [[Rcpp::export]]
List nuts_sample_xptr(NumericVector theta0, SEXP f_ptr, SEXP g_ptr,
                      int n_iter, NumericVector eps_seq, NumericMatrix inv_M,
                      int max_treedepth = 10, double Delta_max = 1000.0) {
  Target tgt = make_xptr_target(f_ptr, g_ptr);
  const int d = theta0.size();
  VectorXd th = Eigen::Map<VectorXd>(theta0.begin(), theta0.size());

  NumericMatrix trace(n_iter, d);
  NumericVector accept(n_iter), depth(n_iter), nleap(n_iter), div(n_iter), energy(n_iter);

  bool inv_M_fixed = (inv_M.nrow() == 1);
  for (int it = 0; it < n_iter; ++it) {
    VectorXd im(d);
    int row = inv_M_fixed ? 0 : it;
    for (int i = 0; i < d; ++i) im[i] = inv_M(row, i);
    double eps = (eps_seq.size() == 1) ? eps_seq[0] : eps_seq[it];

    List step = nuts_transition(tgt, th, eps, im, max_treedepth, Delta_max);
    NumericVector nt = step["theta"];
    th = Eigen::Map<VectorXd>(nt.begin(), nt.size());
    for (int i = 0; i < d; ++i) trace(it, i) = th[i];
    accept[it] = step["accept_stat"];
    depth[it]  = step["tree_depth"];
    nleap[it]  = step["n_leapfrog"];
    div[it]    = (bool)step["divergent"] ? 1.0 : 0.0;
    energy[it] = step["energy"];
  }
  return List::create(_["theta"] = trace, _["accept_stat"] = accept,
                      _["tree_depth"] = depth, _["n_leapfrog"] = nleap,
                      _["divergent"] = div, _["energy"] = energy,
                      _["n_logdens"] = (double)tgt.n_logdens,
                      _["n_grad"] = (double)tgt.n_grad);
}

// find_reasonable_epsilon (Hoffman & Gelman Alg. 4) in C++, callback tier.
// [[Rcpp::export]]
double find_reasonable_epsilon_cpp(NumericVector theta0, Function f, Function grad_f,
                                   NumericVector inv_M, double eps = 1.0) {
  Target tgt = make_callback_target(f, grad_f);
  const int d = theta0.size();
  VectorXd theta = Eigen::Map<VectorXd>(theta0.begin(), theta0.size());
  VectorXd im = Eigen::Map<VectorXd>(inv_M.begin(), inv_M.size());

  VectorXd r(d);
  for (int i = 0; i < d; ++i) r[i] = ::norm_rand() / std::sqrt(im[i]);
  VectorXd grad = tgt.g(theta);

  double logp0 = tgt.f(theta);
  double joint0 = joint_from_logp(logp0, r, im);

  VectorXd th = theta, rr = r, gr = grad;
  leapfrog(tgt, th, rr, gr, eps, im);
  double joint1 = joint_logdens(tgt, th, rr, im);

  double log_ratio = joint1 - joint0;
  double a = (log_ratio > std::log(0.5)) ? 1.0 : -1.0;
  int count = 1;
  while (!std::isfinite(log_ratio) || a * log_ratio > -a * std::log(2.0)) {
    eps = std::pow(2.0, a) * eps;
    th = theta; rr = r; gr = grad;
    leapfrog(tgt, th, rr, gr, eps, im);
    joint1 = joint_logdens(tgt, th, rr, im);
    log_ratio = joint1 - joint0;
    if (++count > 100) stop("find_reasonable_epsilon: no reasonable eps in 100 steps");
  }
  return eps;
}
