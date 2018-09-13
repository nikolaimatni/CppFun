#ifndef OBJECT_H
#define OBJECT_H

#include "../Eigen-3.3/Eigen/Core"

using Eigen::Vector2d;

class Object
{
 protected:
  Vector2d m_centroid;
  
 public:
  Object(const Vector2d& centroid = Vector2d {0,0}) : m_centroid {centroid} {}
  
  virtual bool collision(const Vector2d &position) = 0;
  //virtual bool collision(const Object &object) = 0;
};

#endif
