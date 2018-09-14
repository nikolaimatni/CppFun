#include "Simulator.h"
#include "../Eigen-3.3/Eigen/Dense"
#include "MPC.h"
#include <random>


using namespace Eigen;
using Lti = LtiParameters<MatrixXd>;
using Sim = Simulator<VectorXd,Lti>;
using namespace std;

void computeRic(MatrixXd& P, const MatrixXd& Q, const MatrixXd& R, const Lti& lti, double precision = 1e-4, int maxIter = 1000)
{
  int count = 0;
  MatrixXd Pm1 = P;
  
  do{
    ++count;
    Pm1 = P;
    MatrixXd Psi =(lti.m_B.transpose()*P*lti.m_B + R).inverse();
    P = Q + lti.m_A.transpose()*P*lti.m_A -
      lti.m_A.transpose()*P*lti.m_B*Psi*lti.m_B.transpose()*P*lti.m_A;
  }while(((Pm1-P).norm()>precision) &&( count < maxIter));

  cout << "Ric precision is: " << (Pm1-P).norm() << "\n";
}

void dynamics (VectorXd& nextState, const VectorXd& state, const VectorXd& input, const VectorXd& disturbance, const Lti& lti)
{
  nextState = VectorXd(3);
  nextState = lti.m_A*state + lti.m_B*input + lti.m_H*disturbance;
}

void controller (VectorXd& controlAction, const VectorXd& state, const Lti& lti)
{
  controlAction = VectorXd(1);
  controlAction = lti.m_K*state;
}

void disturbance (VectorXd& w, const Lti& lti)
{
  w = VectorXd(1);
  w << Sim::randn();
}



int main()
{
  
  MatrixXd A = 1.2*MatrixXd::Random(3,3);
  MatrixXd B = 2*MatrixXd::Random(3,1);
  MatrixXd H = MatrixXd::Random(3,1);
  MatrixXd P = Matrix3d::Identity();
  MatrixXd Q = Matrix3d::Identity();
  MatrixXd R = Matrix<double,1,1>::Identity()*.01;
  MatrixXd K(1,3);
  VectorXd mean(1);
  mean << 0;
   
  Lti lti {A,B,H,K};
  
  cout << P << "\n";
  computeRic(P,Q,R,lti);
  cout << P << "\n";

  lti.m_K = -(lti.m_B.transpose()*P*lti.m_B + R).inverse()
                 *lti.m_B.transpose()*P*lti.m_A;


  
  VectorXd x0;
  x0 = VectorXd(3);
  x0 << 1, 4, 1;
 
  Sim sim {100, x0, dynamics, controller, disturbance, lti};
    
  sim.simulate();
  sim.writeToFile("xlog.tr", "ulog.tr");

  MPC mpc;

  vector<double> sol = mpc.Solve(x0,10,A,B,Q,R,P,A,A,x0,x0);

  for (auto x : sol)
    cout << x << ", ";
  cout << "\n";

  
  return 0;
}
