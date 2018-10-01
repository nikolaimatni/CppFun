#include "Eigen-3.3/Eigen/Core"
#include "FG.h"
#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include <vector>

using CppAD::AD;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// TODO: Set the

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.

//
// MPC class definition implementation.
//
// MPC::MPC() {}
// MPC::~MPC() {}

void MPC::SetupOptions(std::string &options) {
  // options for IPOPT solver

  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";
}

vector<double> NlpMPC::Solve(const VectorXd &state, const VectorXd &coeffs) {

  // vector containing starting points of the different stacked variables x(0:N)
  // and u(0:N-1) in vars (this is obtained from the model to ensure
  // consistency)
  const vector<int> starts = model_.starts();

  // vector of decision variables
  Dvector vars(nvars_);

  // warm start using previous solution, if it exists
  if (warmstart_.size()) {
    // warm start x[t] = x_old[t+1] u[t] = u_old[t+1] for t<N-2
    for (int t = 0; t < N_ - 2; t++) {
      for (int i = 0; i < starts.size(); i++) {
        vars[starts[i] + t] = warmstart_[starts[i] + t + 1];
      }
    }
    // warm start x[N-2] = x_old[N-1] and x[N-1] = 0, and u[N-2] = 0
    for (int i = 0; i < starts.size(); i++) {
      if (i < model_.nx()) {
        vars[starts[i] + N_ - 2] = warmstart_[starts[i] + N_ - 1];
        vars[starts[i] + N_ - 1] = 0;
      } else {
        vars[starts[i] + N_ - 2] = 0;
      }
    }
  } else { // just start with everything equal to zero except x[0] = x0;
    for (int i = 0; i < model_.nx(); i++) {
      vars[starts[i]] = state[i];
    }
    for (int i = model_.nx(); i < nvars_; i++) {
      vars[i] = 0;
    }
  }

  // setup upper and lower box constraints on state and input
  Dvector vars_lowerbound(nvars_);
  Dvector vars_upperbound(nvars_);

  for (int t = 0; t < N_ - 1; t++) {
    for (int i = 0; i < model_.nx(); i++) {
      vars_lowerbound[starts[i] + t] = bounds_.x_low_[i];
      vars_upperbound[starts[i] + t] = bounds_.x_up_[i];
    }
    for (int i = 0; i < model_.nu(); i++) {
      vars_lowerbound[starts[model_.nx() + i] + t] = bounds_.u_low_[i];
      vars_upperbound[starts[model_.nx() + i] + t] = bounds_.u_up_[i];
    }
  }

  for (int i = 0; i < model_.nx(); i++) {
    vars_lowerbound[starts[i] + N_ - 1] = bounds_.x_low_[i];
    vars_upperbound[starts[i] + N_ - 1] = bounds_.x_up_[i];
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(nconstraints_);
  Dvector constraints_upperbound(nconstraints_);
  for (int i = 0; i < nconstraints_; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  for (int i = 0; i < model_.nx(); i++) {
    constraints_lowerbound[starts[i]] = state[i];
    constraints_upperbound[starts[i]] = state[i];
  }

  std::string options;
  SetupOptions(options);

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG>(options, vars, vars_lowerbound,
                                   vars_upperbound, constraints_lowerbound,
                                   constraints_upperbound, model_, solution);

  // Check some of the solution values
  bool ok = true;
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;

  // TODO: Return the first actuator values. The variables can be accessed
  // with `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  vector<double> result;

  // push back first acutation values
  for (int i = model_.nx(); i < starts.size(); i++) {
    result.push_back(solution.x[starts[i]]);
  }

  // result.push_back(solution.x[starts[4]]);
  // result.push_back(solution.x[starts[5]]);

  // push back predicted trajectory to plot in simulator
  for (int t = 1; t < N_ - 1; t++) {
    result.push_back(solution.x[starts[0] + t]);
    result.push_back(solution.x[starts[1] + t]);
  }

  warmstart_ = solution.x;

  return result;
}
