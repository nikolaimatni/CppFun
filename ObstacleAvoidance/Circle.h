#ifndef CIRCLE_H
#define CIRCLE_H

#include "Object.h"

using namespace std;

class Circle : public Object
{
 private:
  double m_radius;

 public:
 Circle(const Vector2d& center = Vector2d {0,0}, double radius = 0.0) :
  Object {center}, m_radius {radius} {}

   virtual bool collision(const Vector2d& position);

   friend ostream& operator<<(ostream &out, const Circle &circle);
  
};

#endif
