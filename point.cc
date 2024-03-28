#include "point.hh"
#include<iostream>
#include<fstream>
#include<string>
#include<sstream>

int Point::number_of_points = 0;

Point::Point(double x, double y, double z)
{
    x_component_distance = x;
    y_component_distance = y;
    z_component_distance = z;

    point_index = number_of_points;

    number_of_points++;
}

Vector Point::convert_point_to_vector()
{
    Vector V(x_component_distance, y_component_distance, z_component_distance);

    return V;
}


int Point::return_point_index() const
{
  return point_index;
}

int Point::return_number_of_points() const
{
  return number_of_points;
}


std::vector<Point> Point::obtain_list_of_all_points(std::string file_name)
{  
   
   double x_coordinate, y_coordinate, z_coordinate;

   std::vector<Point> list_of_all_points ;
   std::vector<std::string> vector_of_lines;

   std::string single_line = "";

  
   std::ifstream file_object(file_name);

   if(file_object.is_open())
    {
            while(std::getline(file_object, single_line))
           
                {
                    vector_of_lines.push_back(single_line);
                }    
      
        file_object.close();


        for(int i = 2; i<vector_of_lines.size()-1 ; i++)
        {
            single_line = vector_of_lines[i];

            std::string temp_string = "";

            for(int j = 1; j<single_line.size() ; j++)
            {
                 temp_string = temp_string + single_line[j];
            }

            std::istringstream iss(temp_string);
             

            iss>>x_coordinate>>y_coordinate>>z_coordinate;

            Point p(x_coordinate,y_coordinate,z_coordinate);

            list_of_all_points.push_back(p);
        }
    }

    return list_of_all_points;

}


double Point::operator[](int component)
{

       if(component == 0)
       {
          return x_component_distance;
       }


       else if(component == 1)
       {
          return y_component_distance;
       }

       else
       {
          return z_component_distance;
       }

}

int Point::get_point_index() const
{
    return point_index;
}