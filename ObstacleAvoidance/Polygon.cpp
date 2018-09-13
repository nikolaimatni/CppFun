#include "Polygon.h"
#include <iostream>

using namespace std;
using Eigen::Vector2d;

bool Polygon::collision(const Vector2d& position)
{
  bool temp = false;
  return temp;
}


void Polygon::setCentroid()
{
  double Cx {0}, Cy {0}, A {0};

  for ( int n= 0; n < m_vertices.size()-1;  ++n){
    Cx += (m_vertices[n](0) + m_vertices[n+1](0))*(m_vertices[n](0)*m_vertices[n+1](1) - m_vertices[n+1](0)*m_vertices[n](1));
    Cy += (m_vertices[n](1) + m_vertices[n+1](1))*(m_vertices[n](0)*m_vertices[n+1](1) - m_vertices[n+1](0)*m_vertices[n](1));
    A += .5*(m_vertices[n](0)*m_vertices[n+1](1) - m_vertices[n+1](0)*m_vertices[n](1));
  }

  Cx /= (6*A);
  Cy /= (6*A);
 
  m_centroid(0) = Cx;
  m_centroid(1) = Cy;
  
}

ostream& operator<<(ostream &out, const Polygon &poly)
{
  return out << "I have a polygon with centroid (" << poly.m_centroid.transpose() << ")";
}
