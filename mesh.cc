#include "mesh.hh"

Mesh::Mesh(std::string point_file, std::string face_file, std::string cell_file, std::string boundary_file)
{
    list_of_all_points = Point::obtain_list_of_all_points(point_file);


    list_of_all_faces = Face::generate_list_of_all_faces(face_file, list_of_all_points);
   

    list_of_all_cells = Cell::generate_list_of_all_cells(cell_file, list_of_all_faces);
    

    list_of_all_boundaries = Boundary::generate_list_of_all_boundary_faces(boundary_file, list_of_all_faces);
 

   cell_neighbour_assignment();
   face_interpolation_factor_delta_assignment();
 

}

void Mesh::cell_neighbour_assignment()
{
      for(int i = 0; i < list_of_all_cells.size(); i++)
   {
      list_of_all_cells[i].assign_cell_neighbours();
   }
}

void Mesh::face_interpolation_factor_delta_assignment() //define in mesh class, use list of faces
{
    for(int i=0; i<list_of_all_faces.size(); i++)
    {
      Face *f = &list_of_all_faces[i];

       if(f->get_is_boundary_face() == true )
       {
          int owner_cell_index = f->get_owner_index(); 

          Cell c_o = list_of_all_cells[owner_cell_index];
       
          Vector face_centre = f->get_face_centre();
          Vector owner_cell_centre = c_o.get_cell_centre();

          Vector face_to_owner_cell = owner_cell_centre - face_centre;

          double face_to_owner_magnitude = face_to_owner_cell.magnitude_of_vector();

          double IF = 1.0;
          double DELTA = 1.0/face_to_owner_magnitude;

           f->set_face_interpolation_factor(IF);
           f->set_face_delta(DELTA);
       }

       else
       {
           int owner_cell_index = f->get_owner_index(); 
          int neighbour_cell_index = f->get_neighbour_index();

          Cell c_o = list_of_all_cells[owner_cell_index];
          Cell c_n = list_of_all_cells[neighbour_cell_index];                                    //Do we need pointers here

          Vector neighbour_cell_centre = c_n.get_cell_centre();
          Vector face_centre = f->get_face_centre();
          Vector owner_cell_centre = c_o.get_cell_centre();

          Vector face_to_neighbour_cell = neighbour_cell_centre - face_centre;
          Vector owner_cell_to_neighbour_cell = neighbour_cell_centre - owner_cell_centre;

          double face_to_neighbour_magnitude = face_to_neighbour_cell.magnitude_of_vector();
          double cell_to_neighbour_magnitude = owner_cell_to_neighbour_cell.magnitude_of_vector();

          double IF = face_to_neighbour_magnitude/cell_to_neighbour_magnitude;
          double DELTA = 1.0/cell_to_neighbour_magnitude;

           f->set_face_interpolation_factor(IF);
           f->set_face_delta(DELTA);
       }
    }
}


void Mesh::set_scalar_boundary_conditions(std::vector<std::string>boundary_type, std::vector<double>boundary_value)
{
   for(int i = 0 ;i<list_of_all_boundaries.size(); i++)
   {
      Boundary *b = &list_of_all_boundaries[i];

      if(boundary_type[i] == "dirichlet")
      {
          b->set_is_scalar_dirichlet();
          b->set_scalar_dirichlet_boundary_value(boundary_value[i]);
      }

      else
      {
         b->set_is_scalar_neumann();
         b->set_scalar_neumann_boundary_value(boundary_value[i]);
      }
   }
}

void Mesh::set_vector_boundary_conditions(std::vector<std::string>boundary_type, std::vector<std::vector<double>> boundary_value)
{
   for(int i = 0 ;i<list_of_all_boundaries.size(); i++)
   {
      Boundary *b = &list_of_all_boundaries[i];

      if(boundary_type[i] == "dirichlet")
      {
          b->set_is_vector_dirichlet();
          b->set_vector_dirichlet_boundary_value(boundary_value[i]);
      }

      else
      {
         b->set_is_vector_neumann();
         b->set_vector_neumann_boundary_value(boundary_value[i]);
      }
   }
}


std::vector<Point> Mesh::return_list_of_all_points() const
{
   return list_of_all_points;
}

std::vector<Face> Mesh::return_list_of_all_faces() const
{
   return list_of_all_faces;
}

std::vector<Cell> Mesh::return_list_of_all_cells() const
{
    return list_of_all_cells;
}

std::vector<Boundary> Mesh::return_list_of_all_boundaries() const
{
    return list_of_all_boundaries;
}