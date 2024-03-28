#include "boundary.hh"

#include<vector>
#include<cmath>
#include<sstream>
#include<fstream>
#include<string>
#include<cstring>
#include<iostream>

int Boundary::total_number_of_boundaries = 0;

Boundary::Boundary(std::vector<Face*> list_of_boundary_faces, std::string name_of_boundary)
{
    
  this->list_of_boundary_faces = list_of_boundary_faces;
  this->boundary_name = name_of_boundary;
  boundary_index = total_number_of_boundaries;
  total_number_of_boundaries++ ;

  is_scalar_dirichlet = false;
  is_scalar_neumann = false;
  is_vector_dirichlet = false;
  is_vector_neumann = false;

  set_boundary_faces();  
}

void Boundary::set_boundary_faces()
{  
    for(int i=0; i<list_of_boundary_faces.size(); i++)
    {   
        list_of_boundary_faces[i] -> set_boundary_face();
        list_of_boundary_faces[i] -> set_boundary_index(boundary_index);
    }
}

std::vector<Boundary> Boundary:: generate_list_of_all_boundary_faces(std::string filename, std::vector<Face> &list_of_all_faces)
{
    
    std::string single_line, single_line_2, single_line_3;
    std::vector<std::string> vector_of_lines;
    std::vector<Boundary> list_of_all_boundaries;

    std::ifstream file_object(filename);

    if(file_object.is_open())
    {
        while(std::getline(file_object, single_line))
        {
            vector_of_lines.push_back(single_line);
        }
    }

    file_object.close();
   
    for(int i =2 ;i<vector_of_lines.size() - 1; i = i + 5)
    {   
        int l = i + 1;
        int j = i + 3;
        std::string temp_string = ""; 
        single_line = vector_of_lines[i] ;
        single_line_2 = vector_of_lines[j];
        single_line_3 = vector_of_lines[l];
        
        int n = stoi(single_line_3);

        std::vector<Face*> list_of_boundary_faces;

        for(int k=0; k<single_line_2.size(); k++)
        {
            temp_string = temp_string + single_line_2[k];
        }
         
        int x;
        std::istringstream iss(temp_string);

        for(int k=0; k<n; k++)
        {
          iss>>x;
          list_of_boundary_faces.push_back(&list_of_all_faces[x]);
        }
       
        Boundary b(list_of_boundary_faces, single_line);

        list_of_all_boundaries.push_back(b);

    }
    return list_of_all_boundaries;
}

int Boundary::get_boundary_index() const
{
    return boundary_index;
}

std::vector<Face*> Boundary:: get_list_of_faces() const
{
    return list_of_boundary_faces;
}

std::string Boundary::return_boundary_name() const
{
    return boundary_name;
}

void Boundary::set_is_scalar_dirichlet()
{
    is_scalar_dirichlet = true;
}

void Boundary::set_is_scalar_neumann()
{
    is_scalar_neumann = true;
}

void Boundary::set_is_vector_dirichlet()
{
    is_vector_dirichlet = true;
}

void Boundary::set_is_vector_neumann()
{
    is_vector_neumann = true;
}

void Boundary::set_scalar_dirichlet_boundary_value(double a)
{
   scalar_dirichlet_boundary_value = a;
}

void Boundary::set_scalar_neumann_boundary_value(double a)
{
   scalar_neumann_boundary_value = a;
}

void Boundary::set_vector_dirichlet_boundary_value(std::vector<double> a)
{
   vector_dirichlet_boundary_value = a;
}

void Boundary::set_vector_neumann_boundary_value(std::vector<double> a)
{
   vector_neumann_boundary_value = a;
}

bool Boundary::get_is_scalar_dirichlet()
{
    return is_scalar_dirichlet;
}

bool Boundary::get_is_scalar_neumann()
{
    return is_scalar_neumann;
}

bool Boundary::get_is_vector_dirichlet()
{
    return is_vector_dirichlet;
}

bool Boundary::get_is_vector_neumann()
{
    return is_vector_neumann;
}

double Boundary::get_scalar_dirichlet_boundary_value()
{
   return scalar_dirichlet_boundary_value;
}

double Boundary::get_scalar_neumann_boundary_value()
{
   return scalar_neumann_boundary_value;
}

std::vector<double> Boundary::get_vector_dirichlet_boundary_value()
{
   return vector_dirichlet_boundary_value;
}

std::vector<double> Boundary::get_vector_neumann_boundary_value()
{
   return vector_neumann_boundary_value;
}