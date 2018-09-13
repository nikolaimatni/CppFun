#include "Circle.h"
#include "Polygon.h"

#include <iostream>
#include <vector>

using Eigen::Vector2d;

using namespace std;

int main()
{

  Vector2d me {1,1};

  me(1);
  
  Circle circle {Vector2d {0,0}, 1};

  vector<Vector2d> simplexVertices {Vector2d {0,0}, Vector2d {0,1}, Vector2d {1,0}};

  Polygon simplex {simplexVertices};
  
  cout << "We made a circle with " << circle << "\n";
  cout << "I'm sitting at (" << me.transpose() << ")\n";
  cout << (circle.collision(me) ? "I'm in the cirlce!" : "I'm outside the circle!") << "\n";
  cout << simplex << "\n";
  
  return 0;
}
