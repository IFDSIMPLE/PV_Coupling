#include "const_scalar_field.hh"
#include<Eigen/Dense>

Const_scalar_field::Const_scalar_field(int number_of_cells)
{
    const_field_values.resize(number_of_cells);
    const_old_field_values.resize(number_of_cells);

    const_field_values.setZero();
    const_old_field_values.setZero();
}


void Const_scalar_field::set_const_field_values(double a)
{
   for(int i = 0; i<const_field_values.size(); i++)
   {
         const_field_values(i) = a ;
   }

}

void Const_scalar_field::set_old_const_field_values()
{
  const_old_field_values = const_field_values;
}