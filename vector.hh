#ifndef VECTOR_HH
#define VECTOR_HH

class Vector
{
    private:

    double x_component;
    double y_component;
    double z_component;

    public:

    Vector();                          //Default Constructor
    Vector(double, double, double);    // Constructor

    double dot_product(Vector &object);
    Vector cross_product(Vector &object);
    double magnitude_of_vector();

    Vector& operator=(const Vector &object);
    Vector operator+(Vector object);
    Vector operator-(Vector object);
    Vector operator*(double num);
    Vector operator/(double num);
    double operator[](int component) const;
    double& operator[](int i);

};

#endif