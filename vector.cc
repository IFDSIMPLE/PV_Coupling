#include "vector.hh"
#include<cmath>

Vector::Vector()
{
    x_component = 0;
    y_component = 0;
    z_component = 0;
}

Vector::Vector(double a, double b, double c)
{
    x_component = a;
    y_component = b;
    z_component = c;
}

double Vector::dot_product(Vector &object)
{
    double dot_product_result = ((x_component*(object.x_component)) + (y_component*(object.y_component)) + (z_component*(object.z_component)));

    return dot_product_result;
}

Vector Vector::cross_product(Vector &object)
{
    double cross_product_x_component = ((y_component * (object.z_component)) - (z_component * (object.y_component)));
    double cross_product_y_component = -1*((x_component * (object.z_component)) - (z_component* (object.x_component)));
    double cross_product_z_component = ((x_component * (object.y_component)) - (y_component * (object.x_component)));

    Vector V(cross_product_x_component, cross_product_y_component, cross_product_z_component);

    return V;
}

double Vector::magnitude_of_vector()
{
    double magnitude  = sqrt(pow(x_component,2) + pow(y_component,2) + pow(z_component,2));

    return magnitude;
}

Vector Vector::operator+(Vector object)
{
    double x_component_sum = x_component + (object.x_component);
    double y_component_sum = y_component + (object.y_component);
    double z_component_sum = z_component + (object.z_component);

    Vector V(x_component_sum, y_component_sum, z_component_sum);

    return V;
}

Vector Vector::operator-(Vector object)
{

    double x_component_difference = x_component - (object.x_component);
    double y_component_difference = y_component - (object.y_component);
    double z_component_difference = z_component - (object.z_component);

    Vector V(x_component_difference, y_component_difference, z_component_difference);

    return V;
}

Vector Vector::operator*(double num)
{
    double x_component_multiplication  = x_component * num;
    double y_component_multiplication  = y_component * num;
    double z_component_multiplication  = z_component * num;

    Vector V(x_component_multiplication, y_component_multiplication, z_component_multiplication);

    return V;
}

Vector Vector::operator/(double num)
{
   double x_component_division = x_component / num;
   double y_component_division = y_component / num;
   double z_component_division = z_component / num;

   Vector V(x_component_division, y_component_division, z_component_division);

   return V;

}

Vector& Vector::operator=(const Vector &object)
{
    (*this)[0] = object[0];
    (*this)[1] = object[1];
    (*this)[2] = object[2];
    return (*this);
}


double Vector::operator[](int component) const
{ 
    if(component ==0)
        return x_component;

    else if(component ==1)
    return y_component;

    else
    return z_component;
}

double& Vector::operator[](int i)
{
    if(i == 0)
        return x_component;
    else if(i == 1)
        return y_component;
    else
        return z_component;
}
