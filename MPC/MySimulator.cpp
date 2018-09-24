#include "MySimulator.h"
#include <fstream>
#include <iostream>

using Eigen::VectorXd;
using namespace std;

void MySimulator::Step() {
  int t = time_;

  // instantiate the vectors T() to be populated below

  disturbance_[t] = VectorXd();
  state_[t + 1] = VectorXd();
  input_[t] = VectorXd();

  // Populate m_disturbance[t], m_input[t] and m_state[t+1]
  system_->disturber(disturbance_[t]);
  system_->controller(input_[t], state_[t], disturbance_[t]);
  system_->dynamics(state_[t + 1], state_[t], input_[t], disturbance_[t]);

  // advance time
  ++time_;

  // spit stuff out
  if (verbose_) {
    cout << "x(" << t + 1 << ") = (" << state_[t + 1] << ")\n";
    cout << "u(" << t << ") = " << input_[t] << "\n";
    cout << "w(" << t << ") = " << disturbance_[t] << "\n";
  }
}

void MySimulator::Simulate() {
  for (int i = time_; i < horizon_; ++i) {
    Step();
  }
}

void MySimulator::WriteToFile(const string &state_file,
                              const string &input_file) {
  cout << "writing state and input sequence to files: " << state_file << " and "
       << input_file << "\n";

  ofstream state;
  ofstream input;

  state.open(state_file);
  input.open(input_file);

  for (int i = 0; i < horizon_; ++i) {
    state << state_[i] << ";\n";
    input << input_[i] << ";\n";
  }

  state << state_[horizon_] << ";\n";
  state.close();
  input.close();
}
