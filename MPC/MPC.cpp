#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "../Eigen-3.3/Eigen/Core"
#include <cassert>
#include <vector>

using CppAD::AD;
using namespace std;
using namespace Eigen;

class FG_eval {

private:
  //Pretty self-explanatory.  Define state space system (A,B)
  //LQR cost defined by (Q,R) with terminal cost P
  //Set initial condition to x0, and horizon to N
  //Enforce polytopic constraints on state Fx * x(t) <= bx
  //Enforce polytopic constraints on input Fu * u(t) <= bu
  
  const MatrixXd m_A, m_B, m_Q, m_R, m_P, m_Fx, m_Fu;
  const VectorXd m_x0, m_bx, m_bu;
  const int m_N;
  
  
public:
  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  
  FG_eval(const VectorXd& x0, int N,
	  const MatrixXd& A, const MatrixXd& B,
	  const MatrixXd& Q, const MatrixXd& R, const MatrixXd& P,
	  const MatrixXd& Fx, const MatrixXd& Fu, const VectorXd& bx, const VectorXd bu) : m_x0{x0}, m_N{N}, m_A{A}, m_B{B}, m_Q{Q}, m_R{R}, m_P{P}, m_Fx{Fx}, m_Fu{Fu}, m_bx{bx}, m_bu{bu} {}

  
  void operator()(ADvector& fg, const ADvector& z) {
    //fg[0] contains the objective function
    //fg[1+i], i>=1 contains constraint function that is upper/lower bounded by constraints[i]
    //z contains our optimization variables.  It is partitioned as
    // z = [x(0); x(1); ... ; x(N); u(0); u(1); ... ; u(N-1)];
    
    int nx = m_A.rows();
    int nu = m_B.cols();
    int mx = m_Fx.rows();
    int mu = m_Fu.rows();
    
    int xStart = 0;
    int uStart = xStart + (m_N+1)*nx;

    assert(z.size() == (m_N+1)*nx + m_N*nu && "not enough variables!");
    assert(fg.size()-1 == (m_N)*nx + (m_N+1)*mx + m_N*mu && "we don't have enough constraints!");

    for (int t = 0; t < m_N; t++){
      //unfortunatley ADvectors don't seem to have any of the standard iterator objects, so need to manually extract subvectors
      ADvector xt { extract(z,xStart + t*nx,nx) };
      ADvector xtp1 {extract(z,xStart + (t+1)*nx, nx) };
      ADvector ut { extract(z,uStart + t*nu,nu) };

      // add x(t)' Q x(t) + u(t)' R u(t) to the cost fcn
      fg[0] +=  quadForm(xt,m_Q) + quadForm(ut,m_R);
     
      // ADvector don't have +/- operators, so compute terms needed
      ADvector Ax = matVecMultiply(m_A,xt);
      ADvector Bu = matVecMultiply(m_B,ut);

      
      // Starting indices for dynamics x(t+1) = Ax(t) + Bu(t)
      // Fx x(t) <= bx, and Fu u(t) <= bu, respectively
      int cstIdx = t*nx;
      int stIdx = nx*m_N + t*mx;
      int inIdx = nx*m_N + (m_N+1)*mx + t*mu;

      // x(t+1) = Ax(t) + Bu(t)
      for (int i = 0; i<nx; ++i) {
	fg[1 + t*nx + i] = xtp1[i] - Ax[i] - Bu[i];
      }

      
      // Fx x(t) <= bx
      for (int i=0; i<mx; ++i) {
	fg[1 + stIdx + i] = innerProduct(m_Fx.row(i),xt);
      }

      // Fu u(t) <= bu
      for (int i=0; i<mu; ++i) {
	  fg[1 + inIdx + i] = innerProduct(m_Fu.row(i),ut);
      }
   

    }

    // Need to deal with x(N) explicitly
    ADvector xN = extract(z,m_N*nx,nx);

    // Fx x(N) <= bx
    for (int i = 0; i < mx; i++) {
      fg[1 + nx*m_N + mx*m_N + i] = innerProduct(m_Fx.row(i),xN);
   }
    
    //add x(N)' P x(N) to objective -- terminal cost
    fg[0] += quadForm(xN,m_P);

  }

  AD<double> quadForm(const ADvector& x, const MatrixXd& M)
  {
    //returns x'Mx
    return innerProduct(x,matVecMultiply(M,x));
    }

  ADvector extract(const ADvector& z, int start, int size)
  {
    ADvector subvec(size);
    for (int i = 0; i < size; ++i)
      subvec[i] = z[start+i];

    return subvec;
  }

  AD<double> innerProduct(const ADvector& a, const ADvector& x)
  {
    assert(a.size() == x.size() && "inner product failed: dim(a) != dim(x)");
    AD<double> ip = 0;
    for (int i = 0; i < a.size(); ++i) {
      ip += a[i]*x[i];
    }
    return ip;
  }


  AD<double> innerProduct(const VectorXd& a, const ADvector& x)
  {
    assert(a.size() == x.size() && "inner product failed: dim(a) != dim(x)");
    AD<double> ip = 0;
    for (int i = 0; i < a.size(); ++i) {
      ip += a[i]*x[i];
    }
    return ip;
  }
  
 ADvector matVecMultiply(const MatrixXd& A, const ADvector& x)
  {
    assert(A.cols() == x.size() && "A*x failed because A and x are incompatible");
    int m = A.rows();
    ADvector y(m);

    for (int i = 0; i < m; ++i)
      y[i] = innerProduct(A.row(i),x);

    return y;
  }
 

void setAxEqB( const MatrixXd& A, const ADvector& x, const VectorXd& b, ADvector& fg, int start)
  {
    assert(b.size() == A.rows() && "b = Ax failed becasue b and A are incompatible");

    ADvector Ax = matVecMultiply(A,x);
    
    for (int i = 0; i < A.rows(); ++i)
      fg[start+i] = b[i] - Ax[i];

  }
  
};


//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

VectorXd MPC::Solve(const VectorXd& x0, int N,
			  const MatrixXd& A, const MatrixXd& B,
			  const MatrixXd& Q, const MatrixXd& R, const MatrixXd& P,
			  const MatrixXd& Fx, const MatrixXd& Fu, const VectorXd& bx, const VectorXd& bu) {
  
  bool ok = true;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  int nx = A.rows();
  int nu = B.cols();
  int mx = Fx.rows();
  int mu = Fu.rows();

  // Set the number of model variables (includes both states and inputs). nVars = nx * (N+1) + nu * (N)
  size_t nVars = nx * (N+1) + nu * (N);
  
  // Set the number of constraints
  size_t nDynamicsConstraints = nx * N;
  size_t nStateConstraints = mx * (N+1);
  size_t nInputConstraints = mu * N;
  size_t nConstraints = nDynamicsConstraints + nStateConstraints + nInputConstraints;


   // TODO: warm start with last iteration's solution time shifted by 1
  // for now just set the to zero because hopefully we're doing good control
  Dvector vars(nVars);
  for (int i = 0; i < nVars; i++) {
    vars[i] = 0;
  }

  // warm start initial conditions at their value
  for (int i = 0; i < nx; i++){
    vars[i] = x0(i);
  } 
  
  //No upper lower bounds on variables, we'll enforce those using polytopic constraints Fx*x(t) <= bx, Fu*u(t) <= bu
  Dvector varsLowerBound(nVars);
  Dvector varsUpperBound(nVars);
  for (int i = 0; i < nVars; ++i) {
    varsLowerBound[i] = -1e19;
    varsUpperBound[i] = 1e19;
    
  }

  //constraint x[0] = x0;
   for (int i = 0; i < nx; i++) {
    varsLowerBound[i] = x0(i);
    varsUpperBound[i] = x0(i);
  }
 
  // Lower and upper limits for the constraints -- if upper == lower, then IpOpt converts it to an equality constraint.

   Dvector constraintsLowerBound(nConstraints);
   Dvector constraintsUpperBound(nConstraints);

   for (int t = 0; t < N; t++) {
     int cstIdx = t*nx;
     int stIdx = nx*N + t*mx;
     int inIdx = nx*N + (N+1)*mx + t*mu;

     // 0<= x(t+1) - Ax(t) - Bu(t) <=0
     for (int i = 0; i < nx; i++) {
       constraintsLowerBound[cstIdx+i] = 0;
       constraintsUpperBound[cstIdx+i] = 0;
     }

     // Fx * x(t) <= bx
     for (int i = 0; i < mx; i++) {
       constraintsLowerBound[stIdx + i] = -1e19;
       constraintsUpperBound[stIdx + i] = bx(i);
     }

     // Fu * u(t) <= bu
     for (int i = 0; i < mu; i++) {
       constraintsLowerBound[inIdx + i] = -1e19;
       constraintsUpperBound[inIdx + i] = bu(i);     
     } 
    
   }

   // Need to take care of Fx* x(N) <= bx explicitly
   for (int i = 0; i < mx; i++) {
       constraintsLowerBound[nx*N + mx*N + i] = -1e19;
       constraintsUpperBound[nx*N + mx*N + i] = bx(i);
   }

  // object that computes objective and constraints
  FG_eval fg_eval(x0, N, A, B, Q, R, P, Fx, Fu, bx, bu);

  //
  std::string options;
  // Change number to 12 for most verbose, 0 to turn off, 5 for normal output
  options += "Integer print_level 0\n";
  // This speeds things up
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // uncomment next line if you want 
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
					options, vars, varsLowerBound, varsUpperBound, constraintsLowerBound,
					constraintsUpperBound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  //TODO: need to handle infeasibiliyt more gracefully, right now we just break -- could also shift state/input constraints to soft constraints
  //for now, just kill things
  assert(ok && "MPC solve failed");
  
  // Build a VectorXd out of the first control input and return that to be applied -- TODO: add option to return Nu control inputs to compensate for comp delays  
  VectorXd control(nu);
  for (int i = 0; i < nu; i++)
    control(i) = solution.x[(N+1)*nx + i];
  
  return  control;
}
