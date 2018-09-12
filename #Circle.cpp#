#include "Circle.h"
#include <iostream>

using namespace std;
using Eigen::Vector2d;

bool Circle::collision(const Vector2d& position)
{
  return (position-m_centroid).norm() <= m_radius;
}

ostream& operator<<(ostream &out, const Circle &circle)
{
  return out << "(center, radius) = ((" << circle.m_centroid.transpose() << "), " << circle.m_radius << ")";
}
