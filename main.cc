#include<iostream>
#include<vector>
#include<fstream>
#include<Eigen/Dense>

#include "const_vector_field.hh"
#include "vector.hh"
#include "face.hh"
#include "cell.hh"
#include "boundary.hh"
#include "mesh.hh"

#include "const_scalar_field.hh"
#include "scalar_field.hh"
#include "vector_field.hh"

int main (void)
{ 
  int rows = 40, columns = 40, total_cells = rows*columns; 
  double delta_time = 0.05, nu_val = 0.001, alpha_u = 0.7;

  double x_length = 0.1, y_length = 0.1;
  double delta_x = x_length/columns;
  double delta_y = y_length/rows;

  std::vector<double> temp_ap_coeffs;
  std::vector<double> pressure_corrections_for_flux;
  std::vector<double> pressure_corrections_for_velocity_x(total_cells);
  std::vector<double> pressure_corrections_for_velocity_y(total_cells);
  std::vector<double> cell_fluxes_sum;
  std::ofstream div_u("div_u_cell.txt");

  std::vector<double> x_distance (columns);
  std::vector<double> y_distance (rows);
  
  for(int i=0; i<columns; i++)
  {
     x_distance[i] = (0.5 + i)*delta_x;
  }

    for(int j=0; j<rows; j++)
  {
     y_distance[j] = (0.5 + j)*delta_y;
  }

  Const_scalar_field nu(total_cells);
  nu.set_const_field_values(nu_val);

  Vector_field velocity(total_cells);
  std::vector<double> const_vel = {0.0, 0.0, 0.0};

  Vector const_velocity(const_vel[0], const_vel[1], const_vel[2]);
  velocity.set_const_vector_field_values(const_velocity);

  Scalar_field pressure(total_cells);
  pressure.set_const_field_values(0.0);
   pressure.set_old_const_field_values();

   //std::cout << "Reached here!" << std::endl;

  Mesh m("../points.txt","../faces.txt","../cells.txt","../boundary.txt");

  

  std::vector<Cell> list_of_cells = m.return_list_of_all_cells();
  std::vector<Face> list_of_faces = m.return_list_of_all_faces();

  std::vector<std::string> boundary_types_scalar = {"neumann","neumann","neumann"} ;
  std::vector<double> boundary_values_scalar = {0.0, 0.0, 0.0};

  std::vector<std::string> boundary_types_vector = {"dirichlet","dirichlet","dirichlet"};
  std::vector<std::vector<double>> boundary_values_vector = {{1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

  m.set_scalar_boundary_conditions(boundary_types_scalar, boundary_values_scalar);
  m.set_vector_boundary_conditions(boundary_types_vector, boundary_values_vector);

  std::vector<Boundary> list_of_boundaries = m.return_list_of_all_boundaries();

  std::ofstream sum_initial_residual_ux("sum_residual_initial_ux.txt");
  std::ofstream sum_initial_residual_uy("sum_residual_initial_uy.txt");

  std::ofstream sum_final_residual_ux("sum_residual_final_ux.txt");
  std::ofstream sum_final_residual_uy("sum_residual_final_uy.txt");

  std::ofstream sum_initial_pressure_residual("sum_initial_pressure_residual.txt");
  std::ofstream sum_final_pressure_residual("sum_final_pressure_residual.txt");

  std::ofstream pressure_convergence_initial_h("pressure_convergence.txt");
  std::ofstream x_velocity_convergence_initial("x_velocity_convergence.txt");
  std::ofstream y_velocity_convergence_initial("y_velocity_convergence.txt");
   std::ofstream x_velocity_convergence_final("x_velocity_convergence_final.txt");
  std::ofstream y_velocity_convergence_final("y_velocity_convergence_final.txt");


   int iteration_no = 0;

   while(iteration_no < 200)
   {
     std::cout<<"Iteration Number : "<<iteration_no<<std::endl;
     iteration_no = iteration_no +1; 
     
     velocity.compute_diffusion_matrix(list_of_cells, nu, list_of_boundaries);
     velocity.compute_convection_matrix(list_of_cells, list_of_boundaries, velocity);
     velocity.combine_a_matrices();
     velocity.combine_b_matrices();
     velocity.under_relaxation(list_of_cells, alpha_u);  

     velocity.calculate_initial_residuals(x_distance, y_distance, iteration_no, sum_initial_residual_ux, sum_initial_residual_uy);
     velocity.solve_matrices();
     
     velocity.plot_convergence_initial_x(x_velocity_convergence_initial, iteration_no);
     velocity.plot_convergence_initial_y(y_velocity_convergence_initial, iteration_no);

     velocity.set_old_vector_const_field_values();

     velocity.calculate_final_residuals(x_distance, y_distance, iteration_no, sum_final_residual_ux, sum_final_residual_uy);
     
     for(int i = 0; i < list_of_cells.size(); i++)
     {
      velocity.a_matrix_combined(i,i) *= alpha_u;
     }

     temp_ap_coeffs = velocity.store_ap_coefficients();

     velocity.set_face_and_cell_fluxes(list_of_cells, list_of_boundaries);
 
     pressure.compute_diffusion_matrix(list_of_cells, temp_ap_coeffs, list_of_boundaries, iteration_no);    
     pressure.compute_source_matrix(list_of_cells);  
     pressure.calculate_initial_residuals_p(x_distance, y_distance, iteration_no, sum_initial_pressure_residual);                            
     pressure.combine_and_solve_matrices(list_of_cells);
     pressure.plot_convergence_initial(pressure_convergence_initial_h, iteration_no);
     pressure.set_old_const_field_values();
     
     pressure.calculate_final_residuals_p(x_distance, y_distance, iteration_no, sum_final_pressure_residual);
     
     pressure.compute_flux_correction(list_of_cells, temp_ap_coeffs, list_of_boundaries);

  if(iteration_no == 199)
   {
     for(int i=0; i<y_distance.size(); i++)
     {
       for(int j=0; j<x_distance.size(); j++)
       {
          div_u<<x_distance[j]<<"\t\t"<<y_distance[i]<<"\t\t"<<list_of_cells[j + (i*x_distance.size())].get_sum_of_fluxes_through_cell()<<std::endl;
       }
     }
   }

      Const_vector_field grad_p = pressure.compute_gauss_gradient(list_of_cells);

      for(int i = 0; i < list_of_cells.size(); i++)
      {
         velocity.const_vector_field_values_x(i) -= temp_ap_coeffs[i] * grad_p.const_vector_field_values_x(i);
         velocity.const_vector_field_values_y(i) -= temp_ap_coeffs[i] * grad_p.const_vector_field_values_y(i);
      }
 
      pressure.pressure_under_relax();
  
   }

  velocity.output_vector_matrix_coefficients_to_file(total_cells);
  velocity.output_vector_field_to_file(x_distance, y_distance);
  
  pressure.output_scalar_matrix_coefficients_to_file(total_cells);
  pressure.output_scalar_field_to_file(x_distance, y_distance, list_of_cells);

  div_u.close();

  sum_initial_residual_ux.close();
  sum_initial_residual_uy.close();
  sum_final_residual_ux.close();
  sum_final_residual_uy.close();
  sum_initial_pressure_residual.close();
  sum_final_pressure_residual.close(); 
  pressure_convergence_initial_h.close();
  x_velocity_convergence_initial.close();
  y_velocity_convergence_initial.close();
  x_velocity_convergence_final.close();
  y_velocity_convergence_final.close();

  return 0;
}
