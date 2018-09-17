#ifndef MPC_H
#define MPC_H

#include "../Eigen-3.3/Eigen/Core"
#include <vector>

using namespace std;
using namespace Eigen;

/*! \brief MPC solver for linear system, polytopic constraints, and quadratic
cost

This MPC solver solves the problem

min_{x,u} sum_{t=0}^{N-1} x_t' Q x_t + u_t' R u_t

s.t. x_{t+1} = Ax_t + Bu_t, x_0 = x0

     F_x x_t \leq b_x, F_u u_t \leq b_u

using IPOPT and CppAD.  It uses Eigen to represent vectors and matrices.
 */
class MPC {
private:
  size_t m_N; ///< MPC horizon

public:
  MPC();

  virtual ~MPC();

  /*! \brief Solve the model given an initial state and (optional) reference
    trajectory*/
  VectorXd Solve(const VectorXd &x0, int N, const MatrixXd &A,
                 const MatrixXd &B, const MatrixXd &Q, const MatrixXd &R,
                 const MatrixXd &P, const MatrixXd &Fx, const MatrixXd &Fu,
                 const VectorXd &bx, const VectorXd &bu);
};

#endif /* MPC_H */
