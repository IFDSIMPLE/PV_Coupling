#ifndef BOUNDARY_HH
#define BOUNDARY_HH
#include<string>
#include<vector>

#include "face.hh"

class Boundary
{
    static int total_number_of_boundaries;

    private:
    std::string boundary_name;
    std::vector<Face*> list_of_boundary_faces;
    int boundary_index;

    bool is_scalar_dirichlet;
    bool is_scalar_neumann;

    bool is_vector_dirichlet;
    bool is_vector_neumann;

    double scalar_dirichlet_boundary_value;
    double scalar_neumann_boundary_value;

    std::vector<double> vector_dirichlet_boundary_value;
    std::vector<double> vector_neumann_boundary_value;

    public:
    Boundary(std::vector<Face*> list_of_boundary_faces, std::string name_of_boundary);

    void set_boundary_faces(); 

    void set_is_scalar_dirichlet();
    void set_is_scalar_neumann();
    void set_is_vector_dirichlet();
    void set_is_vector_neumann();
    void set_scalar_dirichlet_boundary_value(double);
    void set_scalar_neumann_boundary_value(double);
    void set_vector_dirichlet_boundary_value(std::vector<double>);
    void set_vector_neumann_boundary_value(std::vector<double>);

    bool get_is_scalar_dirichlet();
    bool get_is_scalar_neumann();
    bool get_is_vector_dirichlet();
    bool get_is_vector_neumann();
    double get_scalar_dirichlet_boundary_value();
    double get_scalar_neumann_boundary_value();
    std::vector<double> get_vector_dirichlet_boundary_value();
    std::vector<double> get_vector_neumann_boundary_value();

    int get_boundary_index() const;
    std::vector<Face*> get_list_of_faces() const;
    std::string return_boundary_name() const;

    static std::vector<Boundary> generate_list_of_all_boundary_faces(std::string filename, std::vector<Face> &list_of_all_faces);
};

#endif