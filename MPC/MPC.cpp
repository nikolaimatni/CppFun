#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "../Eigen-3.3/Eigen/Core"
#include <cassert>

using CppAD::AD;
using namespace std;
using namespace Eigen;

// TODO: Set the timestep length and duration
size_t N = 0;
double dt = 0;

class FG_eval {

private:
  const MatrixXd m_A, m_B, m_Q, m_R, m_P, m_Fx, m_Fu;
  const VectorXd m_x0, m_bx, m_bu;
  int m_N;
  
  
public:
  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  
  FG_eval(const VectorXd& x0, int N,
	  const MatrixXd& A, const MatrixXd& B,
	  const MatrixXd& Q, const MatrixXd& R, const MatrixXd& P,
	  const MatrixXd& Fx, const MatrixXd& Fu, const VectorXd& bx, const VectorXd bu) : m_x0{x0}, m_N{N}, m_A{A}, m_B{B}, m_Q{Q}, m_R{R}, m_P{P}, m_Fx{Fx}, m_Fu{Fu}, m_bx{bx}, m_bu{bu} {}

  
  void operator()(ADvector& fg, const ADvector& x) {
    /*
    fg[0] = x[0]*x[0] + x[1]*x[1];
    int m = m_A.rows();
    int n = m_A.cols();
    assert(m == m_b.size() && "y and A are not compatible");

    setAxEqB(m_A,x,m_b,fg,1);*/
    
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.
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

vector<double> MPC::Solve(const VectorXd& x0, int N,
			  const MatrixXd& A, const MatrixXd& B,
			  const MatrixXd& Q, const MatrixXd& R, const MatrixXd& P,
			  const MatrixXd& Fx, const MatrixXd& Fu, const VectorXd& bx, const VectorXd bu) {
  
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // Set the number of model variables (includes both states and inputs). n_vars = nx * N + nu * (N-1)
  size_t n_vars = 3;
  // Set the number of constraints (excluding the dynamics). This is weird, seems to be (nx+nu) * N, but we'll see.
  size_t n_constraints = 3;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  for (int i = 0; i < n_vars; ++i) {
    vars_lowerbound[i] = -1e19;
    vars_upperbound[i] = 1e19;
    
  }
  // TODO: Set lower and upper limits for variables.

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
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
  options += "Integer print_level  5\n";
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
  
  return  {solution.x[0], solution.x[1], solution.x[2]};
}
