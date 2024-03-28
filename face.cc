#include "face.hh"
#include "vector.hh"
#include "point.hh"
#include<iostream>
#include<string>
#include<sstream>
#include<fstream>

int Face::number_of_faces = 0;

Face::Face(std::vector<Point> list_of_points_for_this_face)
{
   this -> list_of_points_for_this_face = list_of_points_for_this_face;

   boundary_index = -1;
   owner_cell_index = -1;
   neighbour_cell_index = -1;
   is_boundary_face = false; 
   is_face_flux_set = false;

   face_index  = number_of_faces;
   number_of_faces ++ ;

   calculate_face_area();
   calculate_face_normal();
   calculate_face_centre();
}

void Face::calculate_face_area()
{
    
    face_area = 0;

    Vector naive_centroid;
    Vector temporary_sum_vector;

    for(int i = 0; i<list_of_points_for_this_face.size(); i++)
    {
         Point p = list_of_points_for_this_face[i];
         Vector v = p.convert_point_to_vector();

         temporary_sum_vector  = temporary_sum_vector + v;
    }

    naive_centroid = temporary_sum_vector/list_of_points_for_this_face.size();

    for(int i = 0; i<list_of_points_for_this_face.size(); i++)
    {
        Point p1 = list_of_points_for_this_face[i];
        Point p2 = list_of_points_for_this_face[(i + 1) % list_of_points_for_this_face.size()];
     
        Vector v1 = p1.convert_point_to_vector();
        Vector v2 = p2.convert_point_to_vector();

        Vector naive_centre_to_v1 = naive_centroid - v1;
        Vector naive_centre_to_v2 = naive_centroid - v2;

        Vector cross_product = (naive_centre_to_v1).cross_product(naive_centre_to_v2);

        Vector temp_area_vector = cross_product*0.5;

        double temp_area_magnitude = temp_area_vector.magnitude_of_vector();

        face_area = face_area + temp_area_magnitude;

    }
}

void Face::calculate_face_normal()
{
    Vector naive_centroid;
    Vector temporary_sum_vector;

    for(int i = 0; i<list_of_points_for_this_face.size(); i++)
    {
         Point p = list_of_points_for_this_face[i];
         Vector v = p.convert_point_to_vector();

         temporary_sum_vector  = temporary_sum_vector + v;
    }

    naive_centroid = temporary_sum_vector/list_of_points_for_this_face.size();

    for(int i = 0; i<list_of_points_for_this_face.size(); i++)
    {
        Point p1 = list_of_points_for_this_face[i];
        Point p2 = list_of_points_for_this_face[(i + 1) % list_of_points_for_this_face.size()];
     
        Vector v1 = p1.convert_point_to_vector();
        Vector v2 = p2.convert_point_to_vector();

        Vector naive_centre_to_v1 = naive_centroid - v1;
        Vector naive_centre_to_v2 = naive_centroid - v2;

        Vector cross_product = (naive_centre_to_v1).cross_product(naive_centre_to_v2);  //why cant above calculations be included in the bracket (is it because of operator overloading)?
        Vector temp_area_vector = cross_product*0.5;

        face_normal = face_normal + temp_area_vector;
    }

}

void Face::calculate_face_centre()
{
    Vector naive_centroid;
    Vector temporary_sum_vector;
    Vector area_moment;
    Vector sum_of_area_moments;

    for(int i = 0; i<list_of_points_for_this_face.size(); i++)
    {
         Point p = list_of_points_for_this_face[i];
         Vector v = p.convert_point_to_vector();

         temporary_sum_vector  = temporary_sum_vector + v;
    }
    
    naive_centroid = temporary_sum_vector/list_of_points_for_this_face.size();

    for(int i = 0; i<list_of_points_for_this_face.size() ; i++)
    {
        Point p1 = list_of_points_for_this_face[i];
        Point p2 = list_of_points_for_this_face[(i + 1) % list_of_points_for_this_face.size()];

        Vector v1 = p1.convert_point_to_vector();
        Vector v2 = p2.convert_point_to_vector();

        Vector small_triangle_centroid = (v1 + v2 + naive_centroid)/3;

        Vector naive_centre_to_v1 = naive_centroid - v1;
        Vector naive_centre_to_v2 = naive_centroid - v2;

        Vector cross_product = (naive_centre_to_v1).cross_product(naive_centre_to_v2);
        Vector temp_area_vector = cross_product*0.5;

        double temp_area_magnitude = temp_area_vector.magnitude_of_vector();

        area_moment = small_triangle_centroid * temp_area_magnitude;

        sum_of_area_moments = sum_of_area_moments + area_moment;

    }
    
    face_centre = sum_of_area_moments/face_area;

}


std::vector<Face> Face::generate_list_of_all_faces(std::string file_name, std::vector<Point> list_of_all_points)
{   
    std::string single_line;
    std::vector<std::string> vector_of_lines;
    std::vector<Face> list_of_all_faces;
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
       
       int number_of_points_for_the_face = stoi(temp_number_storage);
      
        std::string temp_string = "";
 
        for(int l= k+2; l<single_line.length() - 1; l++)
        {  
            temp_string = temp_string + single_line[l];
        }

         
        std::istringstream iss(temp_string);

        std::vector<int> point_number;

        std::vector<Point> list_of_points;
        
        int x;

        for(int m = 0; m<number_of_points_for_the_face; m++)
        {    
             iss>>x;
             list_of_points.push_back(list_of_all_points[x]); 
        }
        
        Face f(list_of_points);
        
        list_of_all_faces.push_back(f) ; 
        
    }
      return list_of_all_faces;
}

double Face::get_face_area() const
{
    return face_area;
}

Vector Face::get_face_normal() const
{
   return face_normal;
}

Vector Face::get_face_centre() const
{
    return face_centre;
}

int Face::get_face_index() const
{
    return face_index;
}

std::vector<Point> Face::get_list_of_points() const
{
    return list_of_points_for_this_face;
}

void Face::set_owner_cell_index(int a) 
{
    owner_cell_index = a;
}

void Face::set_neighbour_cell_index(int a) 
{
    neighbour_cell_index = a;
}

int Face::get_owner_index() const
{
    return owner_cell_index;
}

int Face::get_neighbour_index() const
{
    return neighbour_cell_index;
}

void Face::set_boundary_face()
{
    is_boundary_face = true;
}


void Face::set_boundary_index(int a)
{
    boundary_index = a;
}

bool Face::get_is_boundary_face() const
{
    return is_boundary_face;
}

bool Face::get_is_face_flux_set() const
{
    return is_face_flux_set;
}

void Face::set_face_interpolation_factor(double a)
{
    face_interpolation_factor = a;
}

double Face::get_face_interpolation_factor() const
{
    return face_interpolation_factor;
}

void Face::set_face_delta(double a)
{
    face_delta_value = a;
}

double Face::get_face_delta() const
{
    return face_delta_value;
}

int Face::get_boundary_index() const
{
    return boundary_index;
}

void Face::set_face_flux(double a)
{
    is_face_flux_set = true;
    face_flux = a;
}

double Face::get_face_flux() const
{
    return face_flux;
}