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

  Vec is class used to represent vectors, Param is a class or struct containing information about dynamics, controller, disturbances.  Must be compatible with Vec.  For example, if Vec is Eigen::VectorXd, then Param could be used to specify state-space parameters that are represented as Eigen::MatrixXd;

  We assume that Vec has a default construct Vec() -- this is needed to instantiate state, input, disturbance before they are populated in their respective functions, as described below

  dynamicsFcn(Vec& nextState, Vec& state, Vec& input, Vec& disturbance, Param& params): usese state, input, disturbance, pararams to compute nextState

  controlFcn(Vec& input, Vec& state, Param& params): uses state, params to compute input

  disturbFcn(Vec& disturbance, Param& params): uses params to compute disturbance

*/

#define VecParamT template<typename Vec, typename Param>
#define VecParamC template<class Vec, class Param>

VecParamT
using dynamicsFcn = void(*)(Vec&, const Vec&, const Vec&, const Vec&, const Param&);

VecParamT
using controlFcn = void(*)(Vec&, const Vec&, const Param&);

VecParamT
using disturbFcn = void(*)(Vec&, const Param&);

VecParamC
class Simulator
{
  
  using dynamicsFcn = void(*)(Vec&, const Vec&, const Vec&, const Vec&, const Param&);
  using controlFcn = void(*)(Vec&, const Vec&, const Param&);
  using disturbFcn = void(*)(Vec&, const Param&);
  
private:
  static default_random_engine s_generator;
  static normal_distribution<double> s_distribution;
  
  int m_horizon; // how long to run the simulation for
  int m_time; // simulation time
  bool m_verbose = true;
  vector<Vec> m_state; // vector of states
  vector<Vec> m_input; // vector of inputs
  vector<Vec> m_disturbance; // vector of disturbances
  dynamicsFcn m_dynamics; // pointer to dynamics function
  controlFcn m_controller; // pointer to controller function
  disturbFcn m_disturber; // pointer to disturber function
  Param m_params; // struct with parameters needed by above functions

public:
  
  static double randn();
  
  Simulator(int horizon, const Vec& x0, dynamicsFcn dynamics=nullptr, controlFcn controller=nullptr, disturbFcn disturber=nullptr, const Param& params=0 );

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

VecParamC 
default_random_engine Simulator<Vec,Param>::s_generator = default_random_engine {time(NULL)};

VecParamC
normal_distribution<double> Simulator<Vec,Param>::s_distribution = normal_distribution<double>(0,1);

VecParamT
double Simulator<Vec,Param>::randn()  { return s_distribution(s_generator); }

VecParamT
Simulator<Vec,Param>::Simulator(int horizon, const Vec& x0, dynamicsFcn dynamics, controlFcn controller, disturbFcn disturber, const Param& params)
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

VecParamT
void Simulator<Vec,Param>::step()
{
  int t = m_time;

  //instantiate the vectors T() to be populated below
  m_disturbance[t] = Vec();
  m_state[t + 1] = Vec();
  m_input[t] = Vec();

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

VecParamT
void Simulator<Vec,Param>::simulate()
{
  for (int i = m_time; i < m_horizon; ++i)
    {
      step();
    }

  int N = m_horizon;
}

VecParamT
void Simulator<Vec,Param>::writeToFile(const string& stateFile, const string& inputFile)
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
