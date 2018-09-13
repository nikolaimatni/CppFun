#ifndef POLYGON_H
#define POLYGON_H

#include "Object.h"
#include <vector>

using namespace std;

class Polygon : public Object
{
 private:
  vector<Vector2d> m_vertices;
  void setCentroid();

 public:
  
 Polygon(const vector<Vector2d> &vertices) :  m_vertices {vertices}
  {
    setCentroid();
  }

  Vector2d getCentroid() const { return m_centroid; }

   virtual bool collision(const Vector2d &position);

   friend ostream& operator<<(ostream& out, const Polygon &poly);
};

#endif
