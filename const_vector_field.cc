#include "const_vector_field.hh"


Const_vector_field::Const_vector_field(int n_cells)
{
    const_vector_field_values_x.resize(n_cells);
    const_vector_field_values_y.resize(n_cells);
    const_vector_field_values_z.resize(n_cells);

    const_vector_old_field_values_x.resize(n_cells);
    const_vector_old_field_values_y.resize(n_cells);
    const_vector_old_field_values_y.resize(n_cells);

    const_vector_field_values_x.setZero();
    const_vector_field_values_y.setZero();
    const_vector_field_values_z.setZero();

    const_vector_old_field_values_x.setZero();
    const_vector_old_field_values_y.setZero();
    const_vector_old_field_values_z.setZero();
}

void Const_vector_field::set_const_vector_field_values(Vector obj)
{
   Vector const_values = obj;
 
    for(int i = 0; i<const_vector_field_values_x.rows(); i++)
    { 
        const_vector_field_values_x(i) = const_values[0] ;
        const_vector_field_values_y(i) = const_values[1] ;
        const_vector_field_values_z(i) = const_values[2] ;
    }
}

void Const_vector_field::set_old_vector_const_field_values()
{

        const_vector_old_field_values_x = const_vector_field_values_x;
        const_vector_old_field_values_y = const_vector_field_values_y;
        const_vector_old_field_values_z = const_vector_field_values_z;

}