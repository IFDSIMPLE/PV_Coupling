#include "cell.hh"
#include "vector.hh"
#include<iostream>
#include<cmath>
#include<string>
#include<sstream>
#include<fstream>

int Cell::total_number_of_cells = 0;

Cell::Cell(std::vector<Face*> list_of_all_faces_for_this_cell)
{
    this -> list_of_all_faces_for_this_cell = list_of_all_faces_for_this_cell;

    cell_index = total_number_of_cells;

    total_number_of_cells ++;

    calculate_cell_volume(); 
    calculate_cell_centre();
    set_owner_neighbour_index_for_faces();

    sum_of_fluxes = 0.0;
}


void Cell::calculate_cell_volume()
{
   double cell_volume_temp = 0;  
   Vector cell_naive_centroid;
   double volume_c = 0;


   for(int i=0; i<list_of_all_faces_for_this_cell.size(); i++)
   {
       Vector face_centre_each_face =  list_of_all_faces_for_this_cell[i]->get_face_centre();

       cell_naive_centroid = cell_naive_centroid + face_centre_each_face;
   }

        cell_naive_centroid = cell_naive_centroid/list_of_all_faces_for_this_cell.size(); 
    
   for(int i = 0; i<list_of_all_faces_for_this_cell.size(); i++)
   {
      Face *f = list_of_all_faces_for_this_cell[i];

      std::vector<Point> list_of_points_for_this_face = f->get_list_of_points();

      for(int j=0; j<list_of_points_for_this_face.size(); j++)
      {
        Point p1 = list_of_points_for_this_face[j];
        Point p2 = list_of_points_for_this_face[(j + 1) % list_of_points_for_this_face.size()];
     
        Vector v1 = p1.convert_point_to_vector();
        Vector v2 = p2.convert_point_to_vector();

        Vector face_centre_this_face = f->get_face_centre();  

        Vector v1_to_face_centroid = v1 - face_centre_this_face;
        Vector v2_to_face_centroid = v2 - face_centre_this_face;

        Vector face_centroid_to_naive_centroid = cell_naive_centroid - face_centre_this_face;
 
        Vector cross_product_result = v1_to_face_centroid.cross_product(v2_to_face_centroid);

        double dot_product_result  = fabs(cross_product_result.dot_product(face_centroid_to_naive_centroid))*(1.0/6.0);
     
        volume_c = volume_c + dot_product_result ;

      }    
      cell_volume = volume_c;
   }
}


void Cell::calculate_cell_centre()
{   
    Vector cell_naive_centroid;
    Vector aivi;
    Vector sigma_aivi;
    double cell_volume_temp;
    double cell_volume_small;

   for(int i=0; i<list_of_all_faces_for_this_cell.size(); i++)
   {
       Vector face_centre_each_face =  list_of_all_faces_for_this_cell[i]->get_face_centre();

       cell_naive_centroid = cell_naive_centroid + face_centre_each_face;
   }
       cell_naive_centroid = cell_naive_centroid/list_of_all_faces_for_this_cell.size(); 


   for(int i = 0; i<list_of_all_faces_for_this_cell.size(); i++)
   {
      Face *f = list_of_all_faces_for_this_cell[i];

      std::vector<Point> list_of_points_for_this_face = f->get_list_of_points();

      for(int j=0; j<list_of_points_for_this_face.size(); j++)
      {
        Vector face_centroid_a = f->get_face_centre();

        Point p1 = list_of_points_for_this_face[j];
        Point p2 = list_of_points_for_this_face[(j + 1) % list_of_points_for_this_face.size()];
     
        Vector v1 = p1.convert_point_to_vector();
        Vector v2 = p2.convert_point_to_vector();

        Vector small_cell_centroid = (v1 + v2 + face_centroid_a + cell_naive_centroid)/4 ;

        Vector face_centre_this_face = f->get_face_centre();

        Vector v1_to_face_centroid = v1 - face_centre_this_face;
        Vector v2_to_face_centroid = v2 - face_centre_this_face;
        Vector face_centroid_to_naive_centroid = cell_naive_centroid - face_centre_this_face;

        Vector cross_product_result = v1_to_face_centroid.cross_product(v2_to_face_centroid);
        double dot_product_result  = cross_product_result.dot_product(face_centroid_to_naive_centroid);

        cell_volume_small = fabs(dot_product_result)*(1.0/6.0);
     
        aivi = small_cell_centroid*cell_volume_small;     
        sigma_aivi = sigma_aivi + aivi;

      }    
   }
        cell_centre = sigma_aivi/cell_volume; //Cell volume calculated in erlier function used
}

void Cell::set_owner_neighbour_index_for_faces()
{
   for(int i=0; i<list_of_all_faces_for_this_cell.size(); i++)
   {
      Face *f = list_of_all_faces_for_this_cell[i];
      Vector this_face_centre = f->get_face_centre();
      Vector this_face_normal = f->get_face_normal();
      Vector this_cell_centre = get_cell_centre();

      Vector cell_centre_to_face_centre = this_face_centre - this_cell_centre;

      double dot_p = cell_centre_to_face_centre.dot_product(this_face_normal);

      if(dot_p >= 0)
      {
         f->set_owner_cell_index(cell_index);
      }

      else
      {
         f->set_neighbour_cell_index(cell_index);
      }
   }
} 

void Cell::assign_cell_neighbours()
{
    for(int i=0; i<list_of_all_faces_for_this_cell.size(); i++)
    {
        Face* f = list_of_all_faces_for_this_cell[i];
        if(f->get_is_boundary_face())
        {
            continue;
        }
        if(cell_index == list_of_all_faces_for_this_cell[i]->get_owner_index())
        {
            list_of_neighbouring_cells.push_back(list_of_all_faces_for_this_cell[i]->get_neighbour_index()) ;
        }

        else
        {
            list_of_neighbouring_cells.push_back(list_of_all_faces_for_this_cell[i]->get_owner_index()) ; 
        }
    }
}

std::vector<Cell> Cell::generate_list_of_all_cells(std::string file_name, std::vector<Face> &list_of_all_faces)
{
    std::string single_line;
    std::vector<std::string> vector_of_lines;
    std::vector<Cell> list_of_all_cells;
    int  k;

    std::ifstream file_object(file_name);

    if(file_object.is_open())
    {
        while(std::getline(file_object, single_line))
        {
            vector_of_lines.push_back(single_line);
        } 
    }

    file_object.close();
     
    for(int i = 2; i<vector_of_lines.size() - 1; i++)
    {    
        single_line = vector_of_lines[i];
        
        std::string temp_number_storage = "";
      
        for(int j = 0; j < single_line.size(); j++)
        {   
            if(single_line[j]=='(')
            {
                break;
            }

            else
            {
                temp_number_storage = temp_number_storage + single_line[j];
            }
        
            k = j; 
        }
       
       int number_of_faces_for_the_cell = stoi(temp_number_storage);
      
       std::string temp_string = "";
 
        for(int l= k+2; l<single_line.length() - 1; l++)
        {  
            temp_string = temp_string + single_line[l];
        }

         
        std::istringstream iss(temp_string);

        std::vector<Face*> list_of_faces;
        
        int x;

        for(int m = 0; m<number_of_faces_for_the_cell; m++)
        {    
             iss>>x;
             list_of_faces.push_back(&list_of_all_faces[x]); 
        }
        
        Cell c(list_of_faces);
        
        list_of_all_cells.push_back(c) ; 
        
    }
      return list_of_all_cells;
}

double Cell::get_cell_volume() const
{
   return cell_volume;
}

Vector Cell::get_cell_centre() const
{
    return cell_centre;
}

std::vector<Face*> Cell::get_list_of_faces() const
{
    return list_of_all_faces_for_this_cell;
}

std::vector<int> Cell::get_list_of_neighbouring_cells() const
{
    return list_of_neighbouring_cells;
}

int Cell::get_cell_index() const
{
   return cell_index;
}

double Cell::get_sum_of_fluxes_through_cell() const
{
   return sum_of_fluxes;
}

void Cell::set_sum_of_fluxes(double a)
{
    sum_of_fluxes = a;
}

void Cell::change_sum_of_fluxes(double a)
{
    sum_of_fluxes = sum_of_fluxes - a;
}