#ifndef FG_H
#define FG_H

#include "Eigen-3.3/Eigen/Core"
#include <cmath>
#include <cppad/cppad.hpp>

using CppAD::AD;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::max;

using ADVec = CPPAD_TESTVECTOR(AD<double>);

class FG {

private:
  int N_, nx_, nu_, delay_;
  vector<int> starts_;

  virtual AD<double> Cost(int t, const ADVec &xt, const ADVec &ut,
                          const ADVec &utp1) = 0;
  virtual AD<double> TerminalCost(const ADVec &xN) = 0;
  virtual ADVec DynamicsF(int t, const ADVec &xt, const ADVec &ut) = 0;

  ADVec x_t(int t, const ADVec &vars) {
    ADVec xt(nx_);
    for (int i = 0; i < nx_; ++i)
      xt[i] = vars[starts_[i] + t];

    return xt;
  }

  ADVec u_t(int t, const ADVec &vars) {
    ADVec ut(nu_);
    for (int i = 0; i < nu_; ++i)
      ut[i] = vars[starts_[nx_ + i] + t];

    return ut;
  }

public:
  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;

  FG(int N, int nx, int nu, int delay)
      : N_(N), nx_(nx), nu_(nu), delay_(delay), starts_{0} {
    for (int i = 1; i <= nx_; ++i)
      starts_.push_back(starts_[i - 1] + N);

    for (int i = 0; i < nu_ - 1; ++i)
      starts_.push_back(starts_[nx_ + i] + N - 1);
  }

  int nx() { return nx_; }
  int nu() { return nu_; }
  int N() { return N_; }

  void operator()(ADVec &fg, const ADVec &vars) {

    // start by setting the cost
    fg[0] = 0;
    for (int t = 0; t < N_ - 2; ++t) {
      fg[0] += Cost(t, x_t(t, vars), u_t(t, vars), u_t(t + 1, vars));
    }
    fg[0] +=
        Cost(N_ - 2, x_t(N_ - 2, vars), u_t(N_ - 2, vars), u_t(N_ - 2, vars));
    fg[0] += TerminalCost(x_t(N_ - 1, vars));

    // set the constraint used to enforce x(0) = x0
    ADVec x0 = x_t(0, vars);
    for (int i = 0; i < nx_; ++i)
      fg[1 + starts_[i]] = vars[starts_[i]];

    // set the constraints used to enforce dynamics;
    for (int t = 0; t < N_ - 1; ++t) {
      ADVec Fxtut = DynamicsF(t, x_t(t, vars), u_t(max(0, t - delay_), vars));
      ADVec xtp1 = x_t(t + 1, vars);

      // enforce things elementwise
      for (int i = 0; i < nx_; ++i)
        fg[2 + starts_[i] + t] = xtp1[i] - Fxtut[i];
    }
  }
  virtual ~FG() {}
};

class FGBikeModel : public FG {
private:
  double dt_, vref_;
  VectorXd coeffs_;
  enum State { X, Y, PSI, V, CTE, EPSI };
  enum Input { DELTA, A };

  const double Lf_ = 2.67;

  virtual AD<double> Cost(int t, const ADVec &xt, const ADVec &ut,
                          const ADVec &utp1) override {
    AD<double> cost(0);

    cost += 3000 * CppAD::pow(xt[CTE], 2);
    cost += 3000 * CppAD::pow(xt[EPSI], 2);
    cost += CppAD::pow(xt[V] - vref_, 2);

    cost += 5 * CppAD::pow(ut[DELTA], 2);
    cost += 5 * CppAD::pow(ut[DELTA], 2);

    cost += 700 * CppAD::pow(ut[DELTA] * xt[V], 2);

    cost += 200 * CppAD::pow(utp1[DELTA] - ut[DELTA], 2);
    cost += 10 * CppAD::pow(utp1[A] - ut[A], 2);

    return cost;
  }

  virtual AD<double> TerminalCost(const ADVec &xN) override {
    AD<double> cost(0);

    cost += 3000 * CppAD::pow(xN[CTE], 2);
    cost += 3000 * CppAD::pow(xN[EPSI], 2);
    cost += CppAD::pow(xN[V] - vref_, 2);

    return cost;
  }
  virtual ADVec DynamicsF(int t, const ADVec &xt, const ADVec &ut) override {
    AD<double> x = xt[X];
    AD<double> y = xt[Y];
    AD<double> psi = xt[PSI];
    AD<double> v = xt[V];
    AD<double> cte = xt[CTE];
    AD<double> epsi = xt[EPSI];
    AD<double> a = ut[A];
    AD<double> delta = ut[DELTA];

    AD<double> f = coeffs_[0] + coeffs_[1] * x + coeffs_[2] * CppAD::pow(x, 2) +
                   coeffs_[3] * CppAD::pow(x, 3);

    AD<double> psides = CppAD::atan(coeffs_[1] + 2 * coeffs_[2] * x +
                                    3 * coeffs_[3] * CppAD::pow(x, 2));

    ADVec fxtut(nx());
    fxtut[X] = x + v * CppAD::cos(psi) * dt_;
    fxtut[Y] = y + v * CppAD::sin(psi) * dt_;
    fxtut[PSI] = psi - v / Lf_ * delta * dt_;
    fxtut[V] = v + a * dt_;
    fxtut[CTE] = (f - y) + v * CppAD::sin(epsi) * dt_;
    fxtut[EPSI] = (psi - psides) - v / Lf_ * delta * dt_;

    return fxtut;
  }

public:
  FGBikeModel(int N, int nx, int nu, int delay, double dt, double vref,
              VectorXd coeffs)
      : FG(N, nx, nu, delay), dt_(dt), vref_(vref), coeffs_{coeffs} {}

  virtual ~FGBikeModel() {}
};

#endif
