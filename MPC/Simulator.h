#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <vector>
#include <iostream>
#include <random>
#include <fstream>
#include <ctime>

#include "Parameters.h"
#include "../Eigen-3.3/Eigen/Core"

using namespace std;

/*

This is a class used to simulate dynamics driven by a controller and a disturbance process.  The dynamics are specified by dynamicsFcn, the controller by controllerFcn, and the disturbance process by disturbFcn.  It stores the state, input, and disturbance history in vectors, which can then be printed out to separate files (currently only state and input files are created).

T is class used to represent vectors, S is a class or struct containing information about dynamics, controller, disturbances.  Must be compatible with T.  For example, if T is Eigen::VectorXd, then S could be used to specify state-space parameters that are represented as Eigen::MatrixXd;

We assume that T has a default construct T() -- this is needed to instantiate state, input, disturbance before they are populated in their respective functions, as described below

dynamicsFcn(T& nextState, T& state, T& input, T& disturbance, S& params): usese state, input, disturbance, pararams to compute nextState

controlFcn(T& input, T& state, S& params): uses state, params to compute input

disturbFcn(T& disturbance, S& params): uses params to compute disturbance

*/

template<typename T, typename S>
using dynamicsFcn = void(*)(T&, const T&, const T&, const T&, const S&);

template<typename T, typename S>
using controlFcn = void(*)(T&, const T&, const S&);

template<typename T, typename S>
using disturbFcn = void(*)(T&, const S&);

template<class T, class S>
class Simulator
{
  
  using dynamicsFcn = void(*)(T&, const T&, const T&, const T&, const S&);
  using controlFcn = void(*)(T&, const T&, const S&);
  using disturbFcn = void(*)(T&, const S&);
  
 private:
  static default_random_engine s_generator;
  static normal_distribution<double> s_distribution;
  
  int m_horizon; // how long to run the simulation for
  int m_time; // simulation time
  bool m_verbose = true;
  vector<T> m_state; // vector of states
  vector<T> m_input; // vector of inputs
  vector<T> m_disturbance; // vector of disturbances
  dynamicsFcn m_dynamics; // pointer to dynamics function
  controlFcn m_controller; // pointer to controller function
  disturbFcn m_disturber; // pointer to disturber function
  S m_params; // struct with parameters needed by above functions

 public:
  
  static double randn();
  
  Simulator(int horizon, const T& x0, dynamicsFcn dynamics=nullptr, controlFcn controller=nullptr, disturbFcn disturber=nullptr, const S& params=0 );

  //Advance the dynamics one time-step
  void step();

  //Run the simulation from m_time to m_horizon
  void simulate();

  void setDynamics(dynamicsFcn dynamics) { m_dynamics = dynamics; }
  void setController(controlFcn controller) { m_controller = controller; }
  void setDisturber(disturbFcn disturber) { m_disturber = disturber; }
  void setVerbose(bool verbose) {m_verbose = verbose; }
  bool getVerbose() const {return m_verbose; }
  void writeToFile(const string& stateFile, const string& inputFile);

  
 
};

template<class T, class S>
  default_random_engine Simulator<T,S>::s_generator = default_random_engine {time(NULL)};

template<class T, class S>
normal_distribution<double> Simulator<T,S>::s_distribution = normal_distribution<double>(0,1);

template<typename T, typename S>
double Simulator<T,S>::randn()  { return s_distribution(s_generator); }

template<typename T, typename S>
  Simulator<T,S>::Simulator(int horizon, const T& x0, dynamicsFcn dynamics, controlFcn controller, disturbFcn disturber, const S& params)
{  
  m_time = 0;
  m_params = params;
  m_dynamics = dynamics;
  m_controller = controller;
  m_disturber = disturber;
  m_horizon = horizon;

  //Pre-allocate vectors for state, input, disturbance
  m_state.reserve(m_horizon+1);
  m_input.reserve(m_horizon);
  m_disturbance.reserve(m_horizon);

  //Push our initial condition onto the m_state
  m_state.push_back(x0);
}

template<typename T, typename S>
  void Simulator<T,S>::step()
{
  int t = m_time;

  //instantiate the vectors T() to be populated below
  m_disturbance[t] = T();
  m_state[t + 1] = T();
  m_input[t] = T();

  //Populate m_disturbance[t], m_input[t] and m_state[t+1]
  m_disturber(m_disturbance[t], m_params);
  m_controller(m_input[t],m_state[t],m_params);
  m_dynamics(m_state[t + 1], m_state[t],m_input[t], m_disturbance[t], m_params);

  //advance time
  ++m_time;

  //spit stuff out
  if (m_verbose){    
    cout << "x(" << t+1 << ") = (" << m_state[t+1] << ")\n"; 
    cout << "u(" << t << ") = " << m_input[t] << "\n";
    cout << "w(" << t << ") = " << m_disturbance[t] << "\n";
  }
}

template<typename T, typename S>
  void Simulator<T,S>::simulate()
{
  for (int i = m_time; i < m_horizon; ++i)
    {
      step();
    }

  int N = m_horizon;
}

template<typename T, typename S>
  void Simulator<T,S>::writeToFile(const string& stateFile, const string& inputFile)
{
  cout << "writing state and input sequence to files: " << stateFile << " and " << inputFile << "\n";
  ofstream state;
  ofstream input;

  state.open(stateFile);
  input.open(inputFile);

  for (int i = 0; i < m_horizon; ++i){
    state << m_state[i] << ";\n";
    input << m_input[i] << ";\n";
  }

  state << m_state[m_horizon] << ";\n";
  state.close();
  input.close();
  
  
}

#endif
