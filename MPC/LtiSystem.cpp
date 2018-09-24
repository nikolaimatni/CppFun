#include "System.h"
#include <iostream>

void LtiSystem::disturber(VectorXd &disturbance) {
  disturbance = VectorXd(nw_);
  for (int i = 0; i < nw_; ++i)
    disturbance(i) = randn();
}

void LtiSystem::controller(VectorXd &input, const VectorXd &state,
                           const VectorXd &disturbance) {
  input = VectorXd(nu_);
  input = k_ * state;
}

void LtiSystem::dynamics(VectorXd &next_state, const VectorXd &state,
                         const VectorXd &input, const VectorXd &disturbance) {
  next_state = VectorXd(nx_);
  next_state = a_ * state + b_ * input + h_ * disturbance;
}

// Simple approach to computing solution to Ric equation via the Ric recursion;
void LtiSystem::ComputeController(const int nx, double precision,
                                  int max_iter) {
  int count = 0;
  MatrixXd p = MatrixXd(nx, nx);
  for (int i = 0; i < nx; ++i)
    p(i, i) = 1;

  MatrixXd pm1 = p;

  do {
    p = pm1;
    ++count;
    MatrixXd psi = (b_.transpose() * p * b_ + r_).inverse();
    pm1 = q_ + a_.transpose() * p * a_ -
          a_.transpose() * p * b_ * psi * b_.transpose() * p * a_;
  } while (((pm1 - p).norm() > precision) && (count < max_iter));

  assert((pm1 - p).norm() <= precision &&
         "Ric recursion failed to converge in maxIter");
  k_ = -(b_.transpose() * p * b_ + r_).inverse() * b_.transpose() * p * a_;
  std::cout << "controller K is: " << k_ << "\n";
}
