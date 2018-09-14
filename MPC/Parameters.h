#ifndef PARAMETERS_H
#define PARAMTERS_H

//Structs used to store parameters needed by dynamics, controller, disturbance functions used by Simulator

//Parameters specifying an LTI system x(k+1) = m_Ax(k) + m_Bu(k) + m_Hw(k), u(k) = m_Kx(k), w(k)~N(0,I);
template <class Matrix>
struct LtiParameters
{
public:
  Matrix m_A, m_B, m_H, m_K;
};

#endif
