#include "scalar_field.hh"
#include "const_vector_field.hh"
#include "vector.hh"
#include<iostream>
#include<fstream>
#include <vector>

Scalar_field::Scalar_field(int no_of_cells):Const_scalar_field(no_of_cells)
{
    a_matrix_rate_of_change.resize(no_of_cells, no_of_cells);
    a_matrix_diffusion.resize(no_of_cells, no_of_cells);
    a_matrix_convection.resize(no_of_cells, no_of_cells);
    a_matrix_combined.resize(no_of_cells, no_of_cells);

    b_vector_rate_of_change.resize(no_of_cells);
    b_vector_source.resize(no_of_cells);
    b_vector_diffusion.resize(no_of_cells);
    b_vector_convection.resize(no_of_cells);
    b_vector_combined.resize(no_of_cells);

   a_matrix_rate_of_change.setZero();
   a_matrix_diffusion.setZero();
   a_matrix_convection.setZero();
   a_matrix_combined.setZero();

   b_vector_rate_of_change.setZero();
   b_vector_source.setZero();
   b_vector_diffusion.setZero();
   b_vector_convection.setZero();
   b_vector_combined.setZero();

}

void Scalar_field::reset_rate_of_change()
{
    a_matrix_rate_of_change.setZero();
    b_vector_rate_of_change.setZero();
}

void Scalar_field::reset_source()
{
    b_vector_source.setZero();
}

void Scalar_field::reset_diffusion()
{
    a_matrix_diffusion.setZero();
    b_vector_diffusion.setZero();
}

void Scalar_field::reset_convection()
{
    a_matrix_convection.setZero();
    b_vector_convection.setZero();
}

void Scalar_field::reset_combined()
{
    a_matrix_combined.setZero();
    b_vector_combined.setZero();
}

void Scalar_field::update_a_matrix_rate_of_change(double a, int index)
{
   a_matrix_rate_of_change(index, index) = a_matrix_rate_of_change(index, index) + a ;         
}

void Scalar_field::update_b_vector_rate_of_change(double a, int index)
{
    b_vector_rate_of_change(index) = b_vector_rate_of_change(index) + a;
}

void Scalar_field::update_b_vector_source(double a, int index)
{
    b_vector_source(index) = b_vector_source(index) + a;
}

void Scalar_field::update_a_matrix_diffusion(double a, int row_index, int column_index)
{ 
    a_matrix_diffusion(row_index, column_index) = a_matrix_diffusion(row_index, column_index) + a;                
}                                                                         

void Scalar_field::update_b_vector_diffusion(double a, int index)
{
    b_vector_diffusion(index) = a;           
}

void Scalar_field::update_a_matrix_convection(double a, int row_index, int column_index)
{ 
    a_matrix_convection(row_index, column_index) = a_matrix_convection(row_index, column_index) + a;                 
}                                                                         

void Scalar_field::update_b_vector_convection(double a, int index)
{
    b_vector_convection(index) = b_vector_convection(index) + a;           
}

void Scalar_field::combine_a_and_b_matrices()
{
    for(int i = 0; i<a_matrix_combined.rows(); i++)
    {
        for(int j=0; j<a_matrix_combined.cols(); j++)
        {
            a_matrix_combined(i,j) = a_matrix_rate_of_change(i,j) + a_matrix_diffusion(i,j) + a_matrix_convection(i,j);
        }
    }

    for(int k =0 ; k <b_vector_combined.size(); k++)
    {
        b_vector_combined(k) = b_vector_rate_of_change(k) + b_vector_source(k) + b_vector_diffusion(k) + b_vector_convection(k) ;
    }
}

  void Scalar_field::compute_source_matrix(std::vector<Cell> list_of_cells)                                                               
 {  
    reset_source(); 
    
    for(int i = 0; i<list_of_cells.size(); i++)
    {  
        int cell_index = list_of_cells[i].get_cell_index();

        double source_value = list_of_cells[i].get_sum_of_fluxes_through_cell();

        update_b_vector_source(source_value, cell_index);
    }
 } 

void Scalar_field::compute_diffusion_matrix(std::vector<Cell> list_of_cells, std::vector <double> temp_ap_coeffs, std::vector<Boundary> list_of_boundaries, int iteration_no)                                                               
 { 
      reset_diffusion();  
        
      for(int i = 0; i<list_of_cells.size(); i++)
      {     
            double a_diagonal, a_neighbour, b_source;
            
            int cell_index = list_of_cells[i].get_cell_index();

            std::vector<Face*> list_of_faces  =  list_of_cells[cell_index].get_list_of_faces(); 
           
            for(int k = 0; k<list_of_faces.size() ; k++)
            {
                Face *f =  list_of_faces[k];

                int owner_index = f->get_owner_index();
                int neighbour_index = f->get_neighbour_index();

                double interpolation_factor = f->get_face_interpolation_factor();
                double face_delta = f->get_face_delta();
               
                Vector face_normal = f->get_face_normal();

                if(cell_index == neighbour_index)
                {
                    face_normal = face_normal * (-1);
                }

                double face_normal_magnitude = face_normal.magnitude_of_vector();

                if(f-> get_is_boundary_face() == true)
                {
                     
                     int boundary_index = f->get_boundary_index();
                     Boundary *b = &list_of_boundaries[boundary_index];

                     if(b->return_boundary_name() == "empty")
                     {
                        continue;
                     }

                     else
                     {
                            double IF = f->get_face_interpolation_factor();

                            double ap_owner = temp_ap_coeffs[owner_index];
                            double ap_boundary_face = IF * ap_owner;

                        if(b->get_is_scalar_dirichlet() == true)
                        {   
                            double dirichlet_value = b->get_scalar_dirichlet_boundary_value();
                            a_diagonal = -1*ap_boundary_face*face_normal_magnitude*face_delta;

                            update_a_matrix_diffusion(a_diagonal, cell_index, cell_index);  

                            b_source = ap_boundary_face*face_normal_magnitude*face_delta*dirichlet_value;
                            update_b_vector_diffusion(b_source, cell_index) ;
                        }

                        if(b->get_is_scalar_neumann() == true)
                        {   
                            double neumann_gradient =  b->get_scalar_neumann_boundary_value();
                            b_source = -1*ap_owner * face_normal_magnitude * neumann_gradient;

                            update_b_vector_diffusion(b_source, cell_index);
                        }
                     }
                }

                else
                {
                   double ap_owner_cell = temp_ap_coeffs[owner_index];
                   double ap_neighbour_cell = temp_ap_coeffs[neighbour_index];

                   double ap_face = (interpolation_factor*ap_owner_cell) + ((1 - interpolation_factor)*(ap_neighbour_cell)) ;

                   a_diagonal = -1*ap_face * face_normal_magnitude * face_delta;
                   
                   update_a_matrix_diffusion(a_diagonal, cell_index, cell_index);

                   a_neighbour = -1 * a_diagonal;

                   if(owner_index == cell_index)
                   {
                    update_a_matrix_diffusion(a_neighbour, cell_index, neighbour_index);
                   }

                   else
                   {
                    update_a_matrix_diffusion(a_neighbour, cell_index, owner_index);
                   }
                   
                }

            }
      }
 } 
    

void Scalar_field::calculate_initial_residuals_p(std::vector<double> x_distance, std::vector<double> y_distance, int n, std::ofstream &sum_initial_residual_p)
{
     if(n == 199)
  { 
   std::ofstream file_initial_residuals_p("pressure_initial_residuals.txt");

   initial_residuals_p = b_vector_combined - a_matrix_combined*const_field_values;
   
   for(int i=0; i<y_distance.size(); i++)
   {
      for(int j=0; j<x_distance.size(); j++)
      {
         file_initial_residuals_p<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<const_field_values(j + (i*x_distance.size()))<<std::endl;
            
      }
   } 

   file_initial_residuals_p.close();
  } 

  initial_residuals_p = b_vector_combined - a_matrix_combined*const_field_values;

  double sum_r_p = 0.0;

  for(int i=0; i<initial_residuals_p.size(); i++)
  {
     sum_r_p = sum_r_p + initial_residuals_p(i);
  } 

  sum_initial_residual_p<<n<<"\t"<<sum_r_p<<std::endl;

}
        
void Scalar_field::combine_and_solve_matrices(std::vector<Cell> list_of_cells)
{
    combine_a_and_b_matrices();
    int size = list_of_cells.size();

    const_field_values = a_matrix_combined.fullPivHouseholderQr().solve(b_vector_combined);

}

void Scalar_field::pressure_under_relax()
{
    
    const_field_values = const_old_field_values + (0.7) * (const_field_values - const_old_field_values); 
} 

void Scalar_field::calculate_final_residuals_p(std::vector<double> x_distance, std::vector<double> y_distance, int n, std::ofstream &sum_final_residual_p)
{
  if(n == 199)
  { 
   std::ofstream file_final_residuals_p("pressure_final_residuals.txt");

   final_residuals_p = b_vector_combined - a_matrix_combined*const_field_values;
   
   for(int i=0; i<y_distance.size(); i++)
   {
      for(int j=0; j<x_distance.size(); j++)
      {
            file_final_residuals_p<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<const_field_values(j + (i*x_distance.size()))<<std::endl;     
      }
   } 

   file_final_residuals_p.close();
  }


  final_residuals_p = b_vector_combined - a_matrix_combined*const_field_values;

  double sum_r_p = 0.0;

  for(int i=0; i<final_residuals_p.size(); i++)
  {
     sum_r_p = sum_r_p + final_residuals_p(i);
  } 

  sum_final_residual_p<<n<<"\t"<<sum_r_p<<std::endl;
}


void Scalar_field::compute_flux_correction(std::vector<Cell> &list_of_cells ,std::vector<double> temp_ap_coeffs, std::vector<Boundary> list_of_boundaries)
{       
   double ap_face, grad_p_face, new_face_flux_value;

      for(int i = 0; i<list_of_cells.size(); i++)
      {           
            double sum_face_flux = 0.0;

            int cell_index = list_of_cells[i].get_cell_index();

            std::vector<Face*> list_of_faces  =  list_of_cells[cell_index].get_list_of_faces(); 
           
            for(int k = 0; k<list_of_faces.size() ; k++)
            {
                Face *f =  list_of_faces[k];

                int owner_index = f->get_owner_index();
                int neighbour_index = f->get_neighbour_index();

                double interpolation_factor = f->get_face_interpolation_factor();
                double face_delta = f->get_face_delta();
                double face_area = f->get_face_area();
               
                Vector face_normal = f->get_face_normal();

                if(cell_index == neighbour_index)
                {
                    face_normal = face_normal * (-1);
                }

                double face_normal_magnitude = face_normal.magnitude_of_vector();

                if(f-> get_is_boundary_face() == true)
                {
                     
                     int boundary_index = f->get_boundary_index();
                     Boundary *b = &list_of_boundaries[boundary_index];


                     if(b->return_boundary_name() == "empty")
                     {
                        continue;
                     }

                     else
                     {
                            double IF = f->get_face_interpolation_factor();
                            double ap_owner = temp_ap_coeffs[owner_index];

                            grad_p_face = 0;
                            new_face_flux_value = ap_face * grad_p_face;
                     }
                }

                else
                {
                   double ap_owner_cell = temp_ap_coeffs[owner_index];
                   double ap_neighbour_cell = temp_ap_coeffs[neighbour_index];

                   ap_face = (interpolation_factor*ap_owner_cell) + ((1 - interpolation_factor)*(ap_neighbour_cell)) ;
                   grad_p_face = ((const_field_values(neighbour_index) - const_field_values(owner_index))*face_delta)* face_area; 
                   new_face_flux_value =  ap_face*grad_p_face;
                   
                }

              
              f->set_face_flux(new_face_flux_value); 

              sum_face_flux = sum_face_flux + new_face_flux_value;

            }
            double old_sum = list_of_cells[i].get_sum_of_fluxes_through_cell();

            list_of_cells[i].set_sum_of_fluxes(old_sum - sum_face_flux );    
            
      }
}

Const_vector_field Scalar_field::compute_gauss_gradient(const std::vector<Cell> &listOfCells) const
{
    Const_vector_field grad_p(listOfCells.size());
    for(int i = 0; i < listOfCells.size(); i++)
    {
        Vector grad(0, 0, 0);
        std::vector<Face*> listOfFaces = listOfCells[i].get_list_of_faces();
        for(int j = 0; j < listOfFaces.size(); j++)
        {
            Vector sf = listOfFaces[j]->get_face_normal();
            if(listOfFaces[j]->get_is_boundary_face())
            {
                grad = grad + sf * const_field_values(i);
            }
            else
            {
                int N = listOfFaces[j]->get_neighbour_index();
                int P = listOfFaces[j]->get_owner_index();
                double fx = listOfFaces[j]->get_face_interpolation_factor();
                double p_f = fx * const_field_values(P) + (1 - fx) * const_field_values(N);
                if(i == P)
                {
                    grad = grad + sf * p_f;
                }
                else
                {
                    grad = grad - sf * p_f;
                }
            }
        }
        grad_p.const_vector_field_values_x(i) = grad[0] / listOfCells[i].get_cell_volume();
        grad_p.const_vector_field_values_y(i) = grad[1] / listOfCells[i].get_cell_volume();
        grad_p.const_vector_field_values_z(i) = grad[2] / listOfCells[i].get_cell_volume();
    }
    return grad_p;
}

             
void Scalar_field::output_scalar_matrix_coefficients_to_file(double total_cells)
{
     std::ofstream pressure_convection_coeffs_a("pressure_convection_coeffs_a.txt");
     std::ofstream pressure_convection_coeffs_b("pressure_convection_coeffs_b.txt");
     std::ofstream pressure_diffusion_coeffs_a("pressure_diffusion_coeffs_a.txt");
     std::ofstream pressure_diffusion_coeffs_b("pressure_diffusion_coeffs_b.txt");
     std::ofstream pressure_roc_coeffs_a("pressure_roc_coeffs_a.txt");
     std::ofstream pressure_roc_coeffs_b("pressure_roc_coeffs_b.txt");
     std::ofstream pressure_source_coeffs_b("pressure_source_coeffs_b.txt");
     std::ofstream pressure_combined_coeffs_a("pressure_combined_coeffs_a.txt");
     std::ofstream pressure_combined_coeffs_b("pressure_combined_coeffs_b.txt");
    
     for(int i = 0; i<total_cells; i++)
     {
        for(int j=0; j<total_cells; j++)
              {
                if(j < total_cells - 1)
                {
                 pressure_convection_coeffs_a<<a_matrix_convection(i,j)<<"\t";
                 pressure_diffusion_coeffs_a<<a_matrix_diffusion(i,j)<<"\t";
                 pressure_roc_coeffs_a<<a_matrix_rate_of_change(i,j)<<"\t";
                 pressure_combined_coeffs_a<<a_matrix_combined(i,j)<<"\t";

                }

                else
                {
                 pressure_convection_coeffs_a<<a_matrix_convection(i,j)<<"\n";
                 pressure_diffusion_coeffs_a<<a_matrix_diffusion(i,j)<<"\n";
                 pressure_roc_coeffs_a<<a_matrix_rate_of_change(i,j)<<"\n";
                 pressure_combined_coeffs_a<<a_matrix_combined(i,j)<<"\n";

                }
              }
     }

     for(int k=0; k<total_cells; k++)
     {
                pressure_convection_coeffs_b<<b_vector_convection(k)<<"\n";
                 pressure_diffusion_coeffs_b<<b_vector_diffusion(k)<<"\n";
                 pressure_roc_coeffs_b<<b_vector_rate_of_change(k)<<"\n";
                 pressure_source_coeffs_b<<b_vector_source(k)<<"\n";
                 pressure_combined_coeffs_b<<b_vector_combined(k)<<"\n";

     }

   pressure_convection_coeffs_a.close();
   pressure_convection_coeffs_b.close();
   pressure_diffusion_coeffs_a.close();
   pressure_diffusion_coeffs_b.close();
   pressure_roc_coeffs_a.close();
   pressure_roc_coeffs_b.close();
   pressure_combined_coeffs_a.close();
   pressure_combined_coeffs_b.close();
   pressure_source_coeffs_b.close();
}


void Scalar_field::output_scalar_field_to_file(std::vector<double> x_distance, std::vector<double> y_distance, std::vector<Cell> list_of_cells)
{
     std::ofstream scalar_field_profiles("pressure.txt");
     std::ofstream scalar_field_profiles_volume("Pressure_vol.txt");
     std::ofstream pressure_x_line("pressure_x_line.txt");
     std::ofstream pressure_y_line("pressure_y_line.txt");

        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                double cell_v = list_of_cells[ j + i*x_distance.size()].get_cell_volume();
                scalar_field_profiles_volume<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<const_field_values(j + i*x_distance.size())/cell_v<<std::endl;
                scalar_field_profiles<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<const_field_values(j + i*x_distance.size())<<std::endl;
              }
         }

        for(int j=0; j<x_distance.size(); j++)
         {
            pressure_x_line<<x_distance[j]<<"\t"<< const_field_values(((y_distance.size()/2) - 1.0)*x_distance.size() + j)<<std::endl;
         }

         for(int i=0; i<y_distance.size(); i++)
         {
            pressure_y_line<<y_distance[i]<<"\t"<< const_field_values((x_distance.size()/2) + i*x_distance.size())<<std::endl;
         }

        scalar_field_profiles.close();
        pressure_x_line.close();
        pressure_y_line.close();
}     
                
 void Scalar_field::plot_convergence_initial(std::ofstream &plot_convergence_initial_h, int iteration_no)
 {
    Eigen::VectorXd difference = (const_field_values - const_old_field_values);

    for(int i=0; i<const_field_values.size(); i++)
    {
        difference(i) = abs(difference(i));
    }

     double max = difference(0);
    
    for(int i=1; i<const_field_values.size(); i++)
    {   double next_term = difference(i);
        double last_term = difference(i-1);
        
        if(next_term > last_term)
        {
         max = difference(i);
        }
    }
     plot_convergence_initial_h<<iteration_no<<"\t"<<max<<std::endl;
 }               
