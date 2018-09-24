#ifndef MY_SIMULATOR_H
#define MY_SIMULATOR_H

#include "../Eigen-3.3/Eigen/Core"
#include "System.h"
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace Eigen;
using namespace std;

class MySimulator {
private:
  int horizon_;
  int time_;
  bool verbose_ = true;

  System *system_;
  vector<VectorXd> state_;
  vector<VectorXd> input_;
  vector<VectorXd> disturbance_;

public:
  MySimulator(int horizon, const VectorXd &x0, System *system)
      : horizon_{horizon}, time_{0}, system_{system} {
    state_.reserve(horizon_ + 1);
    input_.reserve(horizon_);
    disturbance_.reserve(horizon_);
    state_.push_back(x0);
  }
  void Step();
  void Simulate();
  void WriteToFile(const string &state_file, const string &input_file);
  ~MySimulator() { cout << "Calling MySimulator destructor\n"; }
};

#endif /* MY_SIMULATOR_H */
