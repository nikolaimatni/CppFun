#ifndef MPC_H
#define MPC_H

#include "Eigen-3.3/Eigen/Core"
#include <vector>

using namespace std;

class MPC {
public:
  MPC();

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuations.
  virtual vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
};

class NlpMPC : public MPC {};

#endif /* MPC_H */
