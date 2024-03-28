#ifndef CONST_SCALAR_FIELD_HH
#define CONST_SCALAR_FIELD_HH

#include<Eigen/Dense>

class Const_scalar_field 
{
   public:
   
   Eigen::VectorXd const_field_values;
   Eigen::VectorXd const_old_field_values;

   Const_scalar_field(int);

   void set_const_field_values(double);
   void set_old_const_field_values();
};


#endif