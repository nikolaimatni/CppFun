#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <vector>
#include <iostream>
#include "Parameters.h"
#include "../Eigen-3.3/Eigen/Core"

using namespace std;
using Eigen::MatrixXd;

//Here T is what you are using to represent vectors, S is your Parameter class containing information about dynamics, controller, disturbances
template<typename T, typename S>
using dynamicsFcn = void(*)(T&, const T&, const T&, const T&, const S&);

template<typename T, typename S>
using controlFcn = void(*)(T&, const T&, const S&);

template<typename T>
using disturbFcn = void(*)(T&);

template<class T, class S>
class Simulator
{
  
  using dynamicsFcn = void(*)(T&, const T&, const T&, const T&, const S&);
  using controlFcn = void(*)(T&, const T&, const S&);
  using disturbFcn = void(*)(T&);
  
 private:
  int m_horizon;
  int m_time;
  vector<T> m_state;
  vector<T> m_input;
  vector<T> m_disturbance;
  dynamicsFcn m_dynamics;
  controlFcn m_controller;
  disturbFcn m_disturber;
  S m_params;

 public:
  Simulator(int horizon, const T& x0, dynamicsFcn dynamics=nullptr, controlFcn controller=nullptr, disturbFcn disturber=nullptr, const S& params=0 );
  void step();
  void setController(controlFcn controller) { m_controller = controller; }
  void setDisturber(disturbFcn disturber) { m_disturber = disturber; }
  void simulate();
 
};

template<typename T, typename S>
  Simulator<T,S>::Simulator(int horizon, const T& x0, dynamicsFcn dynamics, controlFcn controller, disturbFcn disturber, const S& params)
{
  cout << "x0 is " << x0 << '\n';
  
  m_time = 0;
  m_params = params;
  m_dynamics = dynamics;
  m_controller = controller;
  m_disturber = disturber;
  m_horizon = horizon;
  
  m_state.reserve(m_horizon+1);
  m_input.reserve(m_horizon);
  m_disturbance.reserve(m_horizon);

  m_state.push_back(x0);
}

template<typename T, typename S>
  void Simulator<T,S>::step()
{
  int t = m_time;
  
  m_disturbance[t] = T();
  m_state[t + 1] = T();
  m_input[t] = T();

  m_disturber(m_disturbance[t]);
  m_controller(m_input[t],m_state[t],m_params);
  m_dynamics(m_state[t + 1], m_state[t],m_input[t], m_disturbance[t], m_params);

  ++m_time;
  cout << "x(" << t+1 << ") = (" << m_state[t+1] << ")\n"; 
  cout << "u(" << t << ") = " << m_input[t] << "\n";
  cout << "w(" << t << ") = " << m_disturbance[t] << "\n";
}

template<typename T, typename S>
  void Simulator<T,S>::simulate()
{
  for (int i = 0; i < m_horizon; ++i)
    {
      step();
    }

  int N = m_horizon;
}
  

#endif
