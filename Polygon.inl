bool Polygon::collision(const Vector2d& position)
{
  return false;
}

void Polygon::setCentroid()
{
  double Cx {0}, Cy {0}, A {0};
  for ( Vector2d& vec : m_vertices){
    Cx += 
    }
  return 
}

ostream& operator<<(ostream &out, const Polygon &poly)
{
  return out << "I have a polygon at  = (" << poly.m_centroid.transpose() << ")";
}
