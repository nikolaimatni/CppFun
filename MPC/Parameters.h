#ifndef PARAMETERS_H
#define PARAMTERS_H

//Base Class, purposefully empty
/*class Parameters
{
 public:
  Parameters() {}
  ~Parameters() {};
  };*/

template <class T>
struct LtiParameters// : public Parameters
{
 public:
  T m_A, m_B, m_H, m_K;
  // LtiParameters() {};
  //LtiParameters(const T& A, const T& B, const T& H, const T& K);
};

/*
template<typename T>
LtiParameters<T>::LtiParameters(const T& A, const T& B, const T& H, const T& K)
{
  m_A = A; m_B = B; m_H = H; m_K = K;
  }*/
#endif
