#ifndef MESH_HH
#define MESH_HH

 #include "cell.hh"
 #include "boundary.hh"

class Mesh
{
    private:
    std::vector<Point> list_of_all_points;
    std::vector<Face> list_of_all_faces;
    std::vector<Cell> list_of_all_cells;
    std::vector<Boundary> list_of_all_boundaries;

    public:
    Mesh(std::string, std::string, std::string, std::string);
 
    void cell_neighbour_assignment();
    void face_interpolation_factor_delta_assignment();
    void face_delta_assignment();

    void set_scalar_boundary_conditions(std::vector<std::string>, std::vector<double>);
    void set_vector_boundary_conditions(std::vector<std::string>, std::vector<std::vector<double>>);
    
    std::vector<Point> return_list_of_all_points() const;
    std::vector<Face> return_list_of_all_faces() const;
    std::vector<Cell> return_list_of_all_cells() const;
    std::vector<Boundary> return_list_of_all_boundaries() const;

};

#endif