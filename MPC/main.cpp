#include <random>
#include "../Eigen-3.3/Eigen/Dense"
#include "MPC.h"
#include "Simulator.h"

using namespace Eigen;
using LtiParams = LtiParameters<MatrixXd>;
using MpcParams = MpcParameters<MatrixXd, MPC>;
using namespace std;

/*
Using this as a simple way of testing my homespun dynamics simulator and MPC
solver
 */

// Simple approach to computing solution to Ric equation via the Ric recursion;
void computeRic(MatrixXd& P, const MatrixXd& Q, const MatrixXd& R,
                const LtiParams& lti, double precision = 1e-6,
                int maxIter = 1000) {
  int count = 0;
  MatrixXd Pm1 = P;

  do {
    P = Pm1;
    ++count;
    MatrixXd Psi = (lti.m_B.transpose() * P * lti.m_B + R).inverse();
    Pm1 = Q + lti.m_A.transpose() * P * lti.m_A -
          lti.m_A.transpose() * P * lti.m_B * Psi * lti.m_B.transpose() * P *
              lti.m_A;
  } while (((Pm1 - P).norm() > precision) && (count < maxIter));

  assert((Pm1 - P).norm() <= precision &&
         "Ric recursion failed to converge in maxIter");
}

// ---------------- LTI Simulator Functions ---------------------- //
void ltiDynamics(VectorXd& nextState, const VectorXd& state,
                 const VectorXd& input, const VectorXd& disturbance,
                 LtiParams& lti) {
  nextState = VectorXd(3);
  nextState = lti.m_A * state + lti.m_B * input + lti.m_H * disturbance;
}

void ltiController(VectorXd& controlAction, const VectorXd& state,
                   LtiParams& lti) {
  controlAction = VectorXd(1);
  controlAction = lti.m_K * state;
}

void ltiDisturbance(VectorXd& w, LtiParams& lti) {
  w = VectorXd(3);
  w << Simulator<VectorXd, MpcParams>::randn(),
      Simulator<VectorXd, MpcParams>::randn(),
      Simulator<VectorXd, MpcParams>::randn();
}

// ------------------ MPC Simulator Functions ---------------//

void mpcController(VectorXd& controlAction, const VectorXd& state,
                   MpcParams& mpc) {
  controlAction = VectorXd(1);
  controlAction =
      mpc.m_mpc.Solve(state, 10, mpc.m_A, mpc.m_B, mpc.m_Q, mpc.m_R, mpc.m_P,
                      mpc.m_Fx, mpc.m_Fu, mpc.m_bx, mpc.m_bu);
}

void mpcDynamics(VectorXd& nextState, const VectorXd& state,
                 const VectorXd& input, const VectorXd& disturbance,
                 MpcParams& mpc) {
  nextState = VectorXd(3);
  nextState = mpc.m_A * state + mpc.m_B * input + mpc.m_H * disturbance;
}

void mpcDisturbance(VectorXd& w, MpcParams& mpc) {
  w = VectorXd(3);
  w << Simulator<VectorXd, MpcParams>::randn(),
      Simulator<VectorXd, MpcParams>::randn(),
      Simulator<VectorXd, MpcParams>::randn();
}

//-------------------- Some simple tests ---------------//

int main() {
  // Set up the system
  MatrixXd A = 1.2 * MatrixXd::Random(3, 3);
  MatrixXd B = 2 * MatrixXd::Random(3, 1);
  MatrixXd H = Matrix3d::Identity() * .2;
  MatrixXd P = Matrix3d::Identity();  // IC for Ric recursion
  MatrixXd Q = Matrix3d::Identity();
  MatrixXd R = Matrix<double, 1, 1>::Identity();
  MatrixXd K(1, 3);  // will assign this after we comput P = DARE(A,B,Q,R)
  int N = 10;        // MPC horizon

  MatrixXd Fx = A * 0;             // no state constraints;
  VectorXd bx = Vector3d::Zero();  // no state constraints
  MatrixXd Fu = VectorXd(2);       // want to enforce -.5 <= u <= .5;
  Fu << 1, -1;
  VectorXd bu = Vector2d::Constant(.5);

  LtiParams ltiP{A, B, H, K};

  MPC mpc;
  MpcParams mpcP{A, B, H, Q, R, P, Fx, Fu, bx, bu, N, mpc};

  // Compute solution P = DARE(A,B,Q,R)
  computeRic(P, Q, R, ltiP);

  // Set controller used by lti simulator using P = DARE(A,B,Q,R)
  ltiP.m_K = -(ltiP.m_B.transpose() * P * ltiP.m_B + R).inverse() *
             ltiP.m_B.transpose() * P * ltiP.m_A;

  VectorXd x0 = VectorXd(3);
  x0 << 1, 2, 3;

  // Set up an LTI simulator and and MPC simulator using the same dynamics
  // The MPC controller has box constraints -.5 <= u <= .5
  // Not that the simulators are driven by different noise processes with
  // covariance given by HH'

  Simulator<VectorXd, MpcParams> simMPC{
      100, x0, mpcDynamics, mpcController, mpcDisturbance, mpcP};

  Simulator<VectorXd, LtiParams> simLTI{
      100, x0, ltiDynamics, ltiController, ltiDisturbance, ltiP};

  // make simulators quiet
  simMPC.setVerbose(false);
  simLTI.setVerbose(false);

  // Run the simulation
  simMPC.simulate();
  simLTI.simulate();

  // Write results to files -- plot these by running python plot_traces.py
  simLTI.writeToFile("lti_xlog.tr", "lti_ulog.tr");
  simMPC.writeToFile("mpc_xlog.tr", "mpc_ulog.tr");

  // Sanity check: MPC with Q_N = P = DARE(A,B,Q,R) and no constraints should
  // give same control action as Kx1
  cout << "Sanity check, these two should be the same: "
       << mpc.Solve(x0, 10, A, B, Q, R, P, A * 0, Fu * 0, x0 * 0, bu) << ", "
       << ltiP.m_K * x0 << "\n";

  return 0;
}
