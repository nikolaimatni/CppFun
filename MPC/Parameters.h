#ifndef PARAMETERS_H
#define PARAMTERS_H
// Structs used to store parameters needed by dynamics, controller, disturbance
// functions used by Simulator

/*! Parameters specifying an LTI system x(k+1) = m_Ax(k) + m_Bu(k) + m_Hw(k),
  u(k) = m_Kx(k), w(k)~N(0,I);*/
template <class Matrix> struct LtiParameters {
public:
  Matrix m_A, m_B, m_H, m_K;
};

/*! Paramters specifying an LTI system as above, except now u(k) is computed
 using the MPC controller mpc, which is specified by parameters weights (Q,R)
 and terminal cost P, polytopic constraint Fx *x <= bx, Fu * u <= bu */
template <class Matrix, class Solver> struct MpcParameters {
public:
  Matrix m_A, m_B, m_H, m_Q, m_R, m_P, m_Fx, m_Fu, m_bx, m_bu;
  int m_N;
  Solver m_mpc;
};
#endif
