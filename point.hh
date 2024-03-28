#ifndef POINT_HH
#define POINT_HH

#include<vector>
#include<string>

#include "vector.hh"

class Point
{
   static int number_of_points;

   private:

   double x_component_distance;
   double y_component_distance;
   double z_component_distance;
   
   int point_index;

   public:
   
   Point(double, double, double);
   Vector convert_point_to_vector();

   int return_point_index() const;
   int return_number_of_points() const;

   int get_point_index() const;

   static std::vector<Point> obtain_list_of_all_points(std::string file_name);

   double operator[](int);

   static std::vector<Point> list_of_all_points;

};

#endif