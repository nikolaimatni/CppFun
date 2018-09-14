#include "Simulator.h"
#include "../Eigen-3.3/Eigen/Dense"
#include "MPC.h"
#include <random>


using namespace Eigen;
using LtiParams = LtiParameters<MatrixXd>;
using MpcParams = MpcParameters<MatrixXd,MPC>;
using namespace std;

void computeRic(MatrixXd& P, const MatrixXd& Q, const MatrixXd& R, const LtiParams& lti, double precision = 1e-6, int maxIter = 1000)
{
  int count = 0;
  MatrixXd Pm1 = P;
  
  do{
    P = Pm1;
    ++count;    
    MatrixXd Psi =(lti.m_B.transpose()*P*lti.m_B + R).inverse();
    Pm1 = Q + lti.m_A.transpose()*P*lti.m_A -
      lti.m_A.transpose()*P*lti.m_B*Psi*lti.m_B.transpose()*P*lti.m_A;
  }while(((Pm1-P).norm()>precision) &&( count < maxIter));

  cout << "Ric precision is: " << (Pm1-P).norm() << "\n";
}

void ltiDynamics (VectorXd& nextState, const VectorXd& state, const VectorXd& input, const VectorXd& disturbance, LtiParams& lti)
{
  nextState = VectorXd(3);
  nextState = lti.m_A*state + lti.m_B*input + lti.m_H*disturbance;
}

void ltiController (VectorXd& controlAction, const VectorXd& state, LtiParams& lti)
{
  controlAction = VectorXd(1);
  controlAction = lti.m_K*state;
}

void ltiDisturbance (VectorXd& w, LtiParams& lti)
{
  w = VectorXd(1);
  w << Simulator<VectorXd,LtiParams>::randn();
}




void mpcController (VectorXd& controlAction, const VectorXd& state, MpcParams& mpc)
{
  controlAction = VectorXd(1);
  controlAction = mpc.m_mpc.Solve(state,10,mpc.m_A,mpc.m_B,mpc.m_Q,mpc.m_R,mpc.m_P,mpc.m_A,mpc.m_A,state,state);
  }


void mpcDynamics (VectorXd& nextState, const VectorXd& state, const VectorXd& input, const VectorXd& disturbance, MpcParams& mpc)
{
  nextState = VectorXd(3);
  nextState = mpc.m_A*state + mpc.m_B*input + mpc.m_H*disturbance;
}


void mpcDisturbance (VectorXd& w, MpcParams& mpc)
{
  w = VectorXd(1);
  w << Simulator<VectorXd,MpcParams>::randn();
}


int main()
{
  
  MatrixXd A = 1.2*MatrixXd::Random(3,3);
  MatrixXd B = 2*MatrixXd::Random(3,1);
  MatrixXd H = MatrixXd::Random(3,1)*0;
  MatrixXd P = Matrix3d::Identity();
  MatrixXd Q = Matrix3d::Identity();
  MatrixXd R = Matrix<double,1,1>::Identity();
  int N = 10;
  MatrixXd K(1,3);
  VectorXd mean(1);
  mean << 0;
   
  LtiParams ltiP {A,B,H,K};

  MPC mpc;
  MpcParams mpcP {A,B,H,Q,R,P,N,mpc};
  
  
  cout << P << "\n";
  computeRic(P,Q,R,ltiP);
  cout << P << "\n";

  ltiP.m_K = -(ltiP.m_B.transpose()*P*ltiP.m_B + R).inverse()
    *ltiP.m_B.transpose()*P*ltiP.m_A;


  
  VectorXd x0;
  x0 = VectorXd(3);
  x0 << 1, 4, 1;
 
  Simulator<VectorXd,MpcParams> simMPC {30, x0, mpcDynamics, mpcController, mpcDisturbance, mpcP};

   Simulator<VectorXd,LtiParams> simLTI  {30, x0, ltiDynamics, ltiController, ltiDisturbance, ltiP};

   
  simLTI.simulate();
  simMPC.simulate();

  simMPC.writeToFile("mpc_xlog.tr", "mpc_ulog.tr");
  simLTI.writeToFile("lti_xlog.tr","lti_ulog.tr");

  //VectorXd sol = mpc.Solve(x0,10,A,B,Q,R,P,A,A,x0,x0);
  
  
  
  
  

  
  return 0;
}
