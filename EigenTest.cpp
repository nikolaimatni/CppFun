#include "Eigen-3.3/Eigen/Core"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::Vector3d;
using namespace std;

int main()
{
  MatrixXd m(2,2);
  m(0,0) = 3;
  m(1,0) = 2.5;
  m(0,1) = -1;
  m(1,1) = m(1,0) + m(0,1);
  cout << m << "\n";

  Vector3d v {1,-3.3, 9};
  cout << v << "\n";
}
