#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <ctime>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

#include "../Eigen-3.3/Eigen/Core"
#include "Parameters.h"

using namespace std;

#define VecParamT template <typename Vec, typename Param>
#define VecParamC template <class Vec, class Param>

VecParamT using dynamicsFcn =
    void(*)(Vec &, const Vec &, const Vec &, const Vec &, Param &);

VecParamT using controlFcn = void(*)(Vec &, const Vec &, Param &);

VecParamT using disturbFcn = void(*)(Vec &, Param &);

/*! \brief A class used to simulate a dynamical system driven by a controller
  and a disturbance process

  This is a class used to simulate dynamics driven by a controller and a
  disturbance process.  The dynamics are specified by dynamicsFcn, the
  controller by controllerFcn, and the disturbance process by disturbFcn.  It
  stores the state, input, and disturbance history in vectors, which can then be
  printed out to separate files (currently only state and input files are
  created).

  Vec is class used to represent vectors, Param is a class or struct containing
  information about dynamics, controller, disturbances.  Must be compatible with
  Vec.  For example, if Vec is Eigen::VectorXd, then Param could be used to
  specify state-space parameters that are represented as Eigen::MatrixXd;

  We assume that Vec has a default construct Vec() -- this is needed to
  instantiate state, input, disturbance before they are populated in their
  respective functions, as described below

  dynamicsFcn(Vec& nextState, Vec& state, Vec& input, Vec& disturbance, Param&
  params): uses state, input, disturbance, pararams to compute nextState

  controlFcn(Vec& input, Vec& state, Param& params): uses state, params to
  compute input

  disturbFcn(Vec& disturbance, Param& params): uses params to compute
  disturbance
*/
VecParamC class Simulator {
  using dynamicsFcn = void (*)(Vec &, const Vec &, const Vec &, const Vec &,
                               Param &);
  using controlFcn = void (*)(Vec &, const Vec &, Param &);
  using disturbFcn = void (*)(Vec &, Param &);

private:
  static default_random_engine s_generator;          //!< random engine
  static normal_distribution<double> s_distribution; //! << normal distribution

  int m_horizon;         //!< how long to run the simulation for
  int m_time;            //!< simulation time
  bool m_verbose = true; //!< Print out state/input/disturbances at each step
  vector<Vec> m_state;   //!< vector of states
  vector<Vec> m_input;   //!< vector of inputs
  vector<Vec> m_disturbance; //!< vector of disturbances
  dynamicsFcn m_dynamics;    //!< pointer to dynamics function
  controlFcn m_controller;   //!< pointer to controller function
  disturbFcn m_disturber;    //!< pointer to disturber function
  Param m_params; //!< struct with parameters needed by above functions

public:
  static double randn(); //!< static helper function to generate N(0,1) r.v.

  Simulator(int horizon, const Vec &x0, dynamicsFcn dynamics = nullptr,
            controlFcn controller = nullptr, disturbFcn disturber = nullptr,
            Param &params =
                0); /*!< \brief Constructer assigning horizon, IC, dynamics,
                      controller, disturbance functions and Parameter object*/

  void step(); //!< Advance the dynamics one time-step

  void simulate(); //!< Run the simulation from m_time to m_horizon

  void reset(int horizon, const Vec &x0, dynamicsFcn dynamics,
             controlFcn controller, disturbFcn disturber,
             Param &params); /*!< \brief reset the simulator (does the same
                                thing as constructor) */

  void setDynamics(dynamicsFcn dynamics) {
    m_dynamics = dynamics;
  } //!< set dynamics function
  void setController(controlFcn controller) {
    m_controller = controller;
  } //!< set controller function
  void setDisturber(disturbFcn disturber) {
    m_disturber = disturber;
  } //!< set disturbance function
  void setParams(Param &params) {
    m_params = params;
  }                                        //!< set Parameters object
  void setX0(Vec &x0) { m_state[0] = x0; } //!< set initial condition
  void setVerbose(bool verbose) {
    m_verbose = verbose;
  } //!< set whether Simulator verbosity
  bool getVerbose() const { return m_verbose; } //!< get Simulator verbosity
  void
  writeToFile(const string &stateFile,
              const string &inputFile); /*!< \brief write simulation state and
                                          input trajectories to stateFile
                                          and inputFile, respectively */
};

VecParamC default_random_engine Simulator<Vec, Param>::s_generator =
    default_random_engine{(unsigned long)time(0)};

VecParamC normal_distribution<double> Simulator<Vec, Param>::s_distribution =
    normal_distribution<double>(0, 1);

VecParamT double Simulator<Vec, Param>::randn() {
  return s_distribution(s_generator);
}

/*!
Creates a new simulator object; assigns appropriate function pointrers, and
reserves vector<Vec>s for state, input and disturbance of size m_horizon +1,
m_horizo, and m_horizon, respectively.
 */
VecParamT Simulator<Vec, Param>::Simulator(int horizon, const Vec &x0,
                                           dynamicsFcn dynamics,
                                           controlFcn controller,
                                           disturbFcn disturber,
                                           Param &params) {
  m_time = 0;
  m_params = params;
  m_dynamics = dynamics;
  m_controller = controller;
  m_disturber = disturber;
  m_horizon = horizon;

  // Pre-allocate vectors for state, input, disturbance
  m_state.reserve(m_horizon + 1);
  m_input.reserve(m_horizon);
  m_disturbance.reserve(m_horizon);

  // Push our initial condition onto the m_state
  m_state.push_back(x0);
}

VecParamT void Simulator<Vec, Param>::reset(int horizon, const Vec &x0,
                                            dynamicsFcn dynamics,
                                            controlFcn controller,
                                            disturbFcn disturber,
                                            Param &params) {
  m_time = 0;
  m_params = params;

  m_dynamics = dynamics;
  m_controller = controller;
  m_disturber = disturber;
  m_horizon = horizon;

  m_state.clear();
  m_input.clear();
  m_disturbance.clear();
  // Pre-allocate vectors for state, input, disturbance
  m_state.reserve(m_horizon + 1);
  m_input.reserve(m_horizon);
  m_disturbance.reserve(m_horizon);

  // Push our initial condition onto the m_state
  m_state.push_back(x0);
}

/*!
This function advances the simulator one time-step.  If m_time = t, it
initializes m_disturbance[t], m_state[t+1], and m_input[t] with empty Vectors
(by calling their default constructor), then runs, in order, the disturbance
function, the controller function, and the dynamics funciton to populate
m_disturbance[t], m_input[t], and m_state[t+1], respectively.  It then advances
time by one, and if m_verbose is set to true, spits out x(t+1), u(t), and w(t).
 */
VecParamT void Simulator<Vec, Param>::step() {
  int t = m_time;

  // instantiate the vectors T() to be populated below

  m_disturbance[t] = Vec();
  m_state[t + 1] = Vec();
  m_input[t] = Vec();

  // Populate m_disturbance[t], m_input[t] and m_state[t+1]
  m_disturber(m_disturbance[t], m_params);
  m_controller(m_input[t], m_state[t], m_params);
  m_dynamics(m_state[t + 1], m_state[t], m_input[t], m_disturbance[t],
             m_params);

  // advance time
  ++m_time;

  // spit stuff out
  if (m_verbose) {
    cout << "x(" << t + 1 << ") = (" << m_state[t + 1] << ")\n";
    cout << "u(" << t << ") = " << m_input[t] << "\n";
    cout << "w(" << t << ") = " << m_disturbance[t] << "\n";
  }
}

VecParamT void Simulator<Vec, Param>::simulate() {
  for (int i = m_time; i < m_horizon; ++i) {
    step();
  }
}

VecParamT void Simulator<Vec, Param>::writeToFile(const string &stateFile,
                                                  const string &inputFile) {
  cout << "writing state and input sequence to files: " << stateFile << " and "
       << inputFile << "\n";
  ofstream state;
  ofstream input;

  state.open(stateFile);
  input.open(inputFile);

  for (int i = 0; i < m_horizon; ++i) {
    state << m_state[i] << ";\n";
    input << m_input[i] << ";\n";
  }

  state << m_state[m_horizon] << ";\n";
  state.close();
  input.close();
}

#endif
