#ifndef PARAMETERS_H
#define PARAMTERS_H

//Base Class, purposefully empty
class Parameters
{
 public:
  Parameters() {}
  ~Paramters() {}
};

template <class T>
class LtiParameters : public Parameters
{
 public:
  T& m_A, m_B, m_H;
  LtiParameters(T& A, T& B, T& H) : m_A {A}, m_B {B}, m_H {H} {}
}
#endif
