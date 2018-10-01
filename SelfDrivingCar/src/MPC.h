#ifndef MPC_H
#define MPC_H

#include "Eigen-3.3/Eigen/Core"
#include "FG.h"
#include <string>
#include <vector>

using Eigen::VectorXd;
using std::string;
using std::vector;
using Model = FG;
using Dvector = CPPAD_TESTVECTOR(double);

struct Bounds {
  vector<double> x_up_, x_low_, u_up_, u_low_;
};

class MPC {
protected:
  size_t N_; // horizon
  Model &model_;
  size_t nvars_, nconstraints_;
  Bounds bounds_;

  CppAD::vector<double> warmstart_;

  virtual void SetupOptions(string &options);

public:
  MPC(size_t N, Model &model, size_t nvars, size_t nconstraints, Bounds bounds)
      : N_{N}, model_(model), nvars_(nvars), nconstraints_(nconstraints),
        bounds_(bounds), warmstart_(0) {}

  // default val for n_vars and n_constraints if we no additional constraints
  // beyond dynamics are imposed

  size_t get_horizon() { return N_; }
  void set_horizon(int n) { N_ = n; }

  void set_model(const FG &new_model) { model_ = new_model; }

  virtual ~MPC() {}

  // Solve the model given an initial state and a reference trajectory.
  // Return the first actuations.
  virtual vector<double> Solve(const VectorXd &state, const VectorXd &ref) = 0;
};

class NlpMPC : public MPC {

private:
  double dt_;

public:
  NlpMPC(size_t N, FG &model, double dt, size_t nvars, size_t nconstraints,
         Bounds bounds)
      : MPC(N, model, nvars, nconstraints, bounds), dt_(dt) {}

  // the reference trajectory is specified by polynomial coefficients contained
  // in coeffs (named ref in base class)
  virtual vector<double> Solve(const VectorXd &state,
                               const VectorXd &coeffs) override;
};

#endif /* MPC_H */
