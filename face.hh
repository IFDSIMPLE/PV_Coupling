#ifndef FACE_HH
#define FACE_HH

#include "point.hh"
#include<vector>

class Face
{
    static int number_of_faces;
    

    private:

    std::vector<Point> list_of_points_for_this_face;
    double face_area;
    Vector face_normal;
    Vector face_centre;

    double face_interpolation_factor;
    double face_delta_value;
    double face_flux;
    
    int face_index;
    int owner_cell_index;
    int neighbour_cell_index;
    int boundary_index;

    bool is_boundary_face;
    bool is_face_flux_set;
   
    public:

    Face(std::vector<Point> list_of_points_for_this_face);
    
    void calculate_face_area();
    void calculate_face_normal();
    void calculate_face_centre();

    void set_owner_cell_index(int);
    void set_neighbour_cell_index(int);
    void set_boundary_face();
    void set_boundary_index(int a);
    void set_face_interpolation_factor(double);
    void set_face_delta(double);
    void set_face_flux(double);
    
    double get_face_area() const;
    Vector get_face_normal() const;
    Vector get_face_centre() const;
    int get_face_index() const;
    int get_owner_index() const;
    int get_neighbour_index() const;  
    bool get_is_boundary_face() const;
    bool get_is_face_flux_set() const;
    double get_face_interpolation_factor() const;
    double get_face_delta() const;
    int get_boundary_index() const;
    double get_face_flux() const;
    

    double operator[](int);

    std::vector<Point> get_list_of_points () const;

    static std::vector<Face> generate_list_of_all_faces(std::string file_name, std::vector<Point> list_of_all_points);
};

#endif