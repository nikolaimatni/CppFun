#include "Simulator.h"
#include "../Eigen-3.3/Eigen/Dense"

using namespace Eigen;
using Lti = LtiParameters<MatrixXd>;

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

void disturbance (VectorXd& w)
{
  w = VectorXd(1);
  w = VectorXd::Random(1)*0;
}

int main()
{
  
  /* RNG stuff
  std::default_random_engine generator;
  std::poisson_distribution<int> distribution(4.1);
  auto poisson = [&] (int) {return distribution(generator);};

  RowVectorXi v = RowVectorXi::NullaryExpr(10, poisson );
  std::cout << v << "\n";
  */
  
  MatrixXd A = 1*MatrixXd::Random(3,3);
  MatrixXd B = 10*MatrixXd::Random(3,1);
  MatrixXd H = MatrixXd::Random(3,1);
  MatrixXd P = Matrix3d::Identity();
  MatrixXd Q = Matrix3d::Identity();
  MatrixXd R = Matrix<double,1,1>::Identity()*.01;
  MatrixXd K(1,3);
 
  Lti lti {A,B,H,K};
  
  computeRic(P,Q,R,lti);

  lti.m_K = -(lti.m_B.transpose()*P*lti.m_B + R).inverse()
                 *lti.m_B.transpose()*P*lti.m_A;


  
  VectorXd x0;
  x0 = VectorXd(3);
  x0 << 4, 5, 6;
  
  Simulator<VectorXd,Lti> sim {100, x0, dynamics, controller, disturbance, lti};
  
  sim.simulate();
  
  return 0;
}
