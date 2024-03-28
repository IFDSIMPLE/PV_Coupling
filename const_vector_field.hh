#ifndef CONST_VECTOR_FIELD_HH
#define CONST_VECTOR_FIELD_HH
#include<Eigen/Dense>
#include "vector.hh"

class Const_vector_field
{
  public:

  Eigen::VectorXd const_vector_field_values_x; 
  Eigen::VectorXd const_vector_field_values_y;
  Eigen::VectorXd const_vector_field_values_z;

  Eigen::VectorXd const_vector_old_field_values_x; 
  Eigen::VectorXd const_vector_old_field_values_y;
  Eigen::VectorXd const_vector_old_field_values_z;

  Const_vector_field(int);

  void set_const_vector_field_values(Vector obj);
  void set_old_vector_const_field_values();

};
#endif