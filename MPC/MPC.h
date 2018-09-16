#ifndef MPC_H
#define MPC_H

#include <vector>
#include "../Eigen-3.3/Eigen/Core"

using namespace std;
using namespace Eigen;

class MPC {
 private:
  size_t m_N;  // horizon

 public:
  MPC();

  virtual ~MPC();

  // Solve the model given an initial state and (optional) reference trajectory
  VectorXd Solve(const VectorXd& x0, int N, const MatrixXd& A,
                 const MatrixXd& B, const MatrixXd& Q, const MatrixXd& R,
                 const MatrixXd& P, const MatrixXd& Fx, const MatrixXd& Fu,
                 const VectorXd& bx, const VectorXd& bu);
};

#endif /* MPC_H */
