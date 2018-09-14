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
    int nx = m_A.rows();
    int nu = m_B.cols();
    int xStart = 0;
    int uStart = xStart + (m_N+1)*nx;

    assert(z.size() == (m_N+1)*nx + m_N*nu && "not enough variables!");
    assert(fg.size()-1 == (m_N)*nx && "we don't have enough constraints to enforce dynamics!");

    /*
    fg[0] = 0;
    ADvector x0 = extract(z,0,3);
    ADvector x1 = extract(z,3,3);
    ADvector u0 = extract(z,6,3);

    ADvector Ax = matVecMultiply(m_A,x0);
    ADvector Bu = matVecMultiply(m_B,u0);

    fg[0] += quadForm(x0,m_Q) + quadForm(u0,m_R) + quadForm(x1,m_P);

    
    fg[1] = x1[0] - Ax[0] - Bu[0];
    fg[2] = x1[1] - Ax[1] - Bu[1];
    fg[3] = x1[2] - Ax[2] - Bu[2];*/

    for (int t = 0; t < m_N; t++){
      //unfortunatley ADvectors don't seem to have any of the standard iterator objects, so need to manually extract subvectors
       ADvector xt { extract(z,xStart + t*nx,nx) };
      ADvector xtp1 {extract(z,xStart + (t+1)*nx, nx) };
      ADvector ut { extract(z,uStart + t*nu,nu) };
      
      fg[0] +=  quadForm(xt,m_Q) + quadForm(ut,m_R);
     
      
      ADvector Ax = matVecMultiply(m_A,xt);
      ADvector Bu = matVecMultiply(m_B,ut);

      for (int m = 0; m<nx; ++m) {
	fg[1 + t*nx + m] = xtp1[m] - Ax[m] - Bu[m];
      }
    }

      //TODO: add polytope constraints Fx * x <= bx, Fu * u <= bu; 

    //terminal cost on state
    fg[0] += quadForm(extract(z,xStart + (m_N)*nx,nx),m_P);

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
			  const MatrixXd& Fx, const MatrixXd& Fu, const VectorXd& bx, const VectorXd bu) {
  
  bool ok = true;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  int nx = A.rows();
  int nu = B.cols();

  // Set the number of model variables (includes both states and inputs). n_vars = nx * (N+1) + nu * (N)
  size_t n_vars = nx * (N+1) + nu * (N);
  // Set the number of constraints (excluding the dynamics). This counts the equality constraints enforcing dynamics = nx*N + other stuff that cant be enforced as box constraints (not implemented yet)
  size_t n_constraints = nx * N;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  // TODO: WARM START
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  // warm start initial conditions at their value
  for (int i = 0; i < nx; i++){
    vars[i] = x0(i);
    cout << "warm start x0("<<i<<")="<<x0(i)<<"\n";
  } 

  //No upper lower bounds
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  for (int i = 0; i < n_vars; ++i) {
    vars_lowerbound[i] = -1e19;
    vars_upperbound[i] = 1e19;
    
  }

  //constraint x[0] = x0;
   for (int i = 0; i < nx; i++) {
    vars_lowerbound[i] = x0(i);
    vars_upperbound[i] = x0(i);
  }
 
  // Lower and upper limits for the constraints
  // Should all be 0 for dynamics, and bx,bu for other stuff (not implemented yet).
   
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

 

  // object that computes objective and constraints
  FG_eval fg_eval(x0, N, A, B, Q, R, P, Fx, Fu, bx, bu);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
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

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
					options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
					constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  assert(ok && "not ok");
  
  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  
  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.

  

  VectorXd control(nu);
  for (int i = 0; i < nu; i++)
    control(i) = solution.x[(N+1)*nx + i];
  
  return  control;
}
