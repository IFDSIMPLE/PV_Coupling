#ifndef CELL_HH
#define CELL_HH

#include "face.hh"
#include<vector>

class Cell
{
    static int total_number_of_cells;

    private:

    std::vector<Face*> list_of_all_faces_for_this_cell; 
    std::vector<int> list_of_neighbouring_cells;
    double cell_volume;
    Vector cell_centre;
    int cell_index;

    double sum_of_fluxes;

    public:

    Cell (std::vector<Face*> list_of_all_faces_for_this_cell); 
    void calculate_cell_volume();
    void calculate_cell_centre();

    void set_owner_neighbour_index_for_faces();
    void set_sum_of_fluxes(double);
    void change_sum_of_fluxes(double);
    void assign_cell_neighbours();

    double get_cell_volume() const;
    Vector get_cell_centre() const;
    int get_cell_index() const;
    double get_sum_of_fluxes_through_cell() const;
    std::vector<int> get_list_of_neighbouring_cells() const;

    std::vector<Face*> get_list_of_faces() const;

    static std::vector<Cell> generate_list_of_all_cells(std::string filename, std::vector<Face> &list_of_all_faces);
};

#endif