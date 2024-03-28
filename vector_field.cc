#include "vector_field.hh"
#include<iostream>
#include<fstream>

Vector_field::Vector_field(int no_of_cells):Const_vector_field(no_of_cells)
{
    a_matrix_rate_of_change.resize(no_of_cells, no_of_cells);
    a_matrix_diffusion.resize(no_of_cells, no_of_cells);
    a_matrix_convection.resize(no_of_cells, no_of_cells);
    a_matrix_combined.resize(no_of_cells, no_of_cells);

    b_vector_rate_of_change.resize(no_of_cells);
    b_vector_source_ux.resize(no_of_cells);
    b_vector_source_uy.resize(no_of_cells);
    b_vector_diffusion_ux.resize(no_of_cells);
    b_vector_diffusion_uy.resize(no_of_cells);
    b_vector_convection_ux.resize(no_of_cells);
    b_vector_convection_uy.resize(no_of_cells);
    b_vector_combined_ux.resize(no_of_cells);
    b_vector_combined_uy.resize(no_of_cells);

   a_matrix_rate_of_change.setZero();
   a_matrix_diffusion.setZero();
   a_matrix_convection.setZero();
   a_matrix_combined.setZero();

   b_vector_rate_of_change.setZero();
   b_vector_source_ux.setZero();
   b_vector_source_ux.setZero();
   b_vector_diffusion_ux.setZero();
   b_vector_diffusion_uy.setZero();
   b_vector_convection_ux.setZero();
   b_vector_convection_uy.setZero();
   b_vector_combined_ux.setZero();
   b_vector_combined_uy.setZero();
}

void Vector_field::reset_rate_of_change()
{
    a_matrix_rate_of_change.setZero();
    b_vector_rate_of_change.setZero();
}

void Vector_field::reset_source()
{
    b_vector_source_ux.setZero();
    b_vector_source_uy.setZero();
}

void Vector_field::reset_diffusion()
{
    a_matrix_diffusion.setZero();
    b_vector_diffusion_ux.setZero();
    b_vector_diffusion_uy.setZero();
}

void Vector_field::reset_convection()
{
    a_matrix_convection.setZero();
    b_vector_convection_ux.setZero();
    b_vector_convection_uy.setZero();

}

void Vector_field::reset_combined()
{
    a_matrix_combined.setZero();
    b_vector_combined_ux.setZero();
    b_vector_combined_uy.setZero();
}


void Vector_field::update_a_matrix_rate_of_change(double a, int index)
{
   a_matrix_rate_of_change(index, index) = a_matrix_rate_of_change(index, index) + a ;         
}

void Vector_field::update_b_vector_rate_of_change(double a, int index)
{
    b_vector_rate_of_change(index) = b_vector_rate_of_change(index) + a;
}

void Vector_field::update_b_vector_source_ux(double a, int index)
{
    b_vector_source_ux(index) = b_vector_source_ux(index) + a;
}

void Vector_field::update_b_vector_source_uy(double a, int index)
{
    b_vector_source_uy(index) = b_vector_source_uy(index) + a;
}

void Vector_field::update_a_matrix_diffusion(double a, int row_index, int column_index)
{ 
    a_matrix_diffusion(row_index, column_index) = a_matrix_diffusion(row_index, column_index) + a;                
}                                                                         

void Vector_field::update_b_vector_diffusion_ux(double a, int index)
{
    b_vector_diffusion_ux(index) = b_vector_diffusion_ux(index) + a;           
}

void Vector_field::update_b_vector_diffusion_uy(double a, int index)
{
    b_vector_diffusion_uy(index) = b_vector_diffusion_uy(index) + a;           
}

void Vector_field::update_a_matrix_convection(double a, int row_index, int column_index)
{ 
    a_matrix_convection(row_index, column_index) = a_matrix_convection(row_index, column_index) + a;                 
}                                                                         

void Vector_field::update_b_vector_convection_ux(double a, int index)
{
    b_vector_convection_ux(index) = b_vector_convection_ux(index) + a;           
}

void Vector_field::update_b_vector_convection_uy(double a, int index)
{
    b_vector_convection_uy(index) = b_vector_convection_uy(index) + a;           
}

void Vector_field::combine_a_matrices()
{
    for(int i = 0; i<a_matrix_combined.rows(); i++)
    {
        for(int j=0; j<a_matrix_combined.cols(); j++)
        {
            a_matrix_combined(i,j) = a_matrix_rate_of_change(i,j) + a_matrix_diffusion(i,j) + a_matrix_convection(i,j);
        }
    }

}

void Vector_field::combine_b_matrices()
{
    for(int k =0 ; k <b_vector_combined_ux.size(); k++)
    {
        b_vector_combined_ux(k) = b_vector_rate_of_change(k) + b_vector_source_ux(k) + b_vector_diffusion_ux(k) + b_vector_convection_ux(k) ;
    }

    for(int k =0 ; k <b_vector_combined_ux.size(); k++)
    {
        b_vector_combined_uy(k) = b_vector_rate_of_change(k) + b_vector_source_uy(k) + b_vector_diffusion_uy(k) + b_vector_convection_uy(k) ;
    }
}

 void Vector_field::compute_rate_of_change_matrix(std::vector<Cell> list_of_cells, double delta_time)                                                               
 {
      reset_rate_of_change();  

      for(int i = 0; i<list_of_cells.size(); i++)
      {  
            int cell_index = list_of_cells[i].get_cell_index();
            double cell_vol = list_of_cells[cell_index].get_cell_volume();
            double temperature_old = const_vector_field_values_x[cell_index];

            double rate_of_change_diagonal_contribution = cell_vol/delta_time;
            double rate_of_change_source_contribution = (cell_vol * temperature_old)/ delta_time;

             update_a_matrix_rate_of_change(rate_of_change_diagonal_contribution, cell_index);
             update_b_vector_rate_of_change(rate_of_change_source_contribution, cell_index);
      }

 } 


void Vector_field::under_relaxation(std::vector<Cell> list_of_cells, double alpha)                                                               
 {  
    double alpha_source_factor = (1 - alpha)/alpha ;

    for(int i=0; i<list_of_cells.size(); i++)
    {
        b_vector_combined_ux(i)  = b_vector_combined_ux(i) + (alpha_source_factor*a_matrix_combined(i,i) * const_vector_field_values_x(i)) ;
        b_vector_combined_uy(i)  = b_vector_combined_uy(i) + (alpha_source_factor*a_matrix_combined(i,i) * const_vector_field_values_y(i)) ; 
        a_matrix_combined(i,i) = a_matrix_combined(i,i) / alpha;
    }
  
 }  


void Vector_field::compute_diffusion_matrix(std::vector<Cell> list_of_cells, Const_scalar_field nu_vector, std::vector<Boundary> list_of_boundaries)                                                               
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
                            double nu_owner = nu_vector.const_field_values(owner_index);
                            double nu_boundary_face = IF * nu_owner;
                        


                        if(b->get_is_vector_dirichlet() == true)
                        {   
                            std::vector<double> dirichlet_values = b->get_vector_dirichlet_boundary_value();
                            
                            a_diagonal = (nu_boundary_face*face_normal_magnitude*face_delta);
                            update_a_matrix_diffusion(a_diagonal, cell_index, cell_index);  
                            
                            for(int q = 0; q<2; q++)
                         {   
                            double dirichlet_value = dirichlet_values[q];
                            b_source = nu_boundary_face*face_normal_magnitude*face_delta*dirichlet_value;

                            if(q == 0)
                            {
                                update_b_vector_diffusion_ux(b_source, cell_index) ;
                            }

                           if(q == 1)
                            {   
                                update_b_vector_diffusion_uy(b_source, cell_index) ;
                            }
                         }
                     
                        }

                        if(b->get_is_vector_neumann() == true)
                        {   
                            double neumann_gradient =  b->get_scalar_neumann_boundary_value();
                            b_source = nu_owner * face_normal_magnitude * neumann_gradient;

                            update_b_vector_diffusion_ux(b_source, cell_index);
                            update_b_vector_diffusion_uy(b_source, cell_index); 
                        }
                     }
                }

                else
                {

                   double nu_owner_cell = nu_vector.const_field_values(owner_index);
                   double nu_neighbour_cell = nu_vector.const_field_values(neighbour_index);

                   double nu_face = (interpolation_factor*nu_owner_cell) + ((1 - interpolation_factor)*(nu_neighbour_cell)) ;

                   a_diagonal = nu_face * face_normal_magnitude * face_delta;

                   update_a_matrix_diffusion(a_diagonal, cell_index, cell_index);

                   a_neighbour = -a_diagonal;
                   
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


 void Vector_field::compute_convection_matrix(std::vector<Cell> list_of_cells, std::vector<Boundary> list_of_boundaries, Vector_field velocity)                                                               
 { 
      reset_convection();  
        
      for(int i = 0; i<list_of_cells.size(); i++)
      {  
            int cell_index = list_of_cells[i].get_cell_index();

             double a_diagonal, a_neighbour, b_source;
             Vector face_velocity;

             std::vector<Face*> list_of_faces  =  list_of_cells[cell_index].get_list_of_faces() ; 
           
             for(int k = 0; k<list_of_faces.size() ; k++)
             {
                Face *f =  list_of_faces[k];
                int P = f->get_owner_index();
                int N = f->get_neighbour_index();
                Vector face_normal = f->get_face_normal();

                if(cell_index == N)
                {
                    face_normal = face_normal *(-1);
                }

                if(f->get_is_boundary_face())
                {
                    Boundary b = list_of_boundaries[f->get_boundary_index()];
                    std::vector<double> vel = b.get_vector_dirichlet_boundary_value();
                    face_velocity = Vector(vel[0], vel[1], vel[2]);
                }
                else
                { 
                    double fx = f->get_face_interpolation_factor();
                    double owner_vx = velocity.const_vector_field_values_x(P);
                    double owner_vy = velocity.const_vector_field_values_x(P);
                    double owner_vz = velocity.const_vector_field_values_z(P);
                    double neighbour_vx = velocity.const_vector_field_values_x(N);
                    double neighbour_vy = velocity.const_vector_field_values_y(N);
                    double neighbour_vz = velocity.const_vector_field_values_z(N);
                    Vector owner_velocity(owner_vx, owner_vy, owner_vz);
                    Vector neighbour_velocity(neighbour_vx, neighbour_vy, neighbour_vz);
                    face_velocity = owner_velocity * fx + neighbour_velocity * (1 - fx);
                }
                double velocity_flux = face_velocity.dot_product(face_normal);
                double face_normal_magnitude = face_normal.magnitude_of_vector();
                int IF;
                double face_delta = f->get_face_delta();
                if(velocity_flux >= 0)
                {
                    IF = 1;
                }

                else
                {
                    IF = 0;
                }
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
                       if(b->get_is_vector_dirichlet() == true)
                        {   
                            double b_source_x = -1*face_velocity[0]*velocity_flux;
                            double b_source_y = -1*face_velocity[1]*velocity_flux;
                            update_b_vector_convection_ux(b_source_x, cell_index) ;
                            update_b_vector_convection_uy(b_source_y, cell_index) ;
                        }

                       if(b->get_is_vector_neumann() == true)
                        {
                            a_diagonal = velocity_flux;
                            update_a_matrix_convection(a_diagonal, cell_index, cell_index);
                            
                           for(int q = 0; q<2; q++)
                           {
                            std::vector<double> neumann_gradients = b->get_vector_neumann_boundary_value();
                            double neumann_gradient =  neumann_gradients[q];

                            b_source = -1*(velocity_flux/face_delta) * neumann_gradient;

                            if(q == 1)
                            {
                              update_b_vector_convection_ux(b_source, cell_index) ;
                            }

                            if(q == 2)
                            {
                              update_b_vector_convection_uy(b_source, cell_index) ;
                            }

                           }  
                        }
                     }
                }

                else
                {
                   a_diagonal = IF*velocity_flux;

                   update_a_matrix_convection(a_diagonal, cell_index, cell_index);

                   a_neighbour = (1 - IF)*velocity_flux;

                   if(P == cell_index)
                   {
                    update_a_matrix_convection(a_neighbour, cell_index, N);
                   }

                   else
                   {
                    update_a_matrix_convection(a_neighbour, cell_index, P);
                   }
                }

             } 
      }
    }  


void Vector_field::calculate_initial_residuals(std::vector<double> x_distance, std::vector<double> y_distance, int n, std::ofstream &sum_initial_residual_ux, 
                                               std::ofstream &sum_initial_residual_uy)
{ 
   if(n == 199)
  { 
   std::ofstream file_initial_residuals_ux("initial_residuals_ux.txt");
   std::ofstream file_initial_residuals_uy("initial_residuals_uy.txt");

   initial_residuals_ux = b_vector_combined_ux - a_matrix_combined*const_vector_field_values_x;
   initial_residuals_uy = b_vector_combined_uy - a_matrix_combined*const_vector_field_values_y;

   for(int i=0; i<y_distance.size(); i++)
   {
      for(int j=0; j<x_distance.size(); j++)
      {
            file_initial_residuals_ux<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<initial_residuals_ux(j + (i*x_distance.size()))<<std::endl;
            file_initial_residuals_uy<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<initial_residuals_uy(j + (i*x_distance.size()))<<std::endl;
      }

   } 

   file_initial_residuals_ux.close();
   file_initial_residuals_uy.close();
  } 
 
  initial_residuals_ux = b_vector_combined_ux - a_matrix_combined*const_vector_field_values_x;
  initial_residuals_uy = b_vector_combined_uy - a_matrix_combined*const_vector_field_values_y;

  double sum_r_ux = 0.0, sum_r_uy = 0.0;

  for(int i=0; i<initial_residuals_ux.size(); i++)
  {
     sum_r_ux = sum_r_ux + initial_residuals_ux(i);
     sum_r_uy = sum_r_uy + initial_residuals_uy(i);
  } 

  sum_initial_residual_ux<<n<<"\t"<<sum_r_ux<<std::endl;
  sum_initial_residual_uy<<n<<"\t"<<sum_r_ux<<std::endl;
 
}    


 void Vector_field::solve_matrices()
{
    const_vector_field_values_x = a_matrix_combined.lu().solve(b_vector_combined_ux);
    const_vector_field_values_y = a_matrix_combined.lu().solve(b_vector_combined_uy);
}  


void Vector_field::calculate_final_residuals(std::vector<double> x_distance, std::vector<double> y_distance, int n, std::ofstream &sum_final_residual_ux, std::ofstream &sum_final_residual_uy)
{
   if(n==199)
 {  
   std::ofstream file_final_residuals_ux("final_residuals_ux.txt");
   std::ofstream file_final_residuals_uy("final_residuals_uy.txt");
   
   final_residuals_ux = b_vector_combined_ux - a_matrix_combined*const_vector_field_values_x;
   final_residuals_uy = b_vector_combined_uy - a_matrix_combined*const_vector_field_values_y;

      for(int i=0; i<y_distance.size(); i++)
   {
      for(int j=0; j<x_distance.size(); j++)
      {
            file_final_residuals_ux<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<final_residuals_ux(j + (i*x_distance.size()))<<std::endl;
            file_final_residuals_uy<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<final_residuals_uy(j + (i*x_distance.size()))<<std::endl;
      }

   } 

   file_final_residuals_ux.close();
   file_final_residuals_uy.close();
 }

  final_residuals_ux = b_vector_combined_ux - a_matrix_combined*const_vector_field_values_x;
  final_residuals_uy = b_vector_combined_uy - a_matrix_combined*const_vector_field_values_y;

  double sum_r_ux = 0.0, sum_r_uy = 0.0;

  for(int i=0; i<initial_residuals_ux.size(); i++)
  {
     sum_r_ux = sum_r_ux + final_residuals_ux(i);
     sum_r_uy = sum_r_uy + final_residuals_uy(i);
  } 

  sum_final_residual_ux<<n<<"\t"<<sum_r_ux<<std::endl;
  sum_final_residual_uy<<n<<"\t"<<sum_r_ux<<std::endl;


}   

void Vector_field::set_face_and_cell_fluxes(std::vector<Cell> &list_of_cells, std::vector<Boundary> list_of_boundaries)
{
    std::vector<double> cell_fluxes;

    for(int i=0; i<list_of_cells.size(); i++)
    {
        double sum_of_face_fluxes_through_cell = 0.0; 
        int cell_index = list_of_cells[i].get_cell_index();

        std::vector<Face*> list_of_faces = list_of_cells[i].get_list_of_faces();

        for(int j=0; j<list_of_faces.size(); j++)
        { 
            double face_flux_value, duplicate_face_flux_value;
            int owner_index = list_of_faces[j] -> get_owner_index();
            int neighbour_index = list_of_faces[j] -> get_neighbour_index();
            double IF = list_of_faces[j] -> get_face_interpolation_factor();
            double x_velocity_neighbour = 0, y_velocity_neighbour = 0, z_velocity_neighbour = 0;

            Vector face_normal = list_of_faces[j] -> get_face_normal();

            if(cell_index != owner_index)
            {
                face_normal = face_normal * -1;  
            }

            if(list_of_faces[j] -> get_is_boundary_face() != true)
            { 
                double x_velocity_owner = const_vector_field_values_x(owner_index);
                double y_velocity_owner = const_vector_field_values_y(owner_index);
                double z_velocity_owner = const_vector_field_values_z(owner_index);

                x_velocity_neighbour = const_vector_field_values_x(neighbour_index);
                y_velocity_neighbour = const_vector_field_values_y(neighbour_index);
                z_velocity_neighbour = const_vector_field_values_z(neighbour_index);
            

                double x_face_velocity = (IF*x_velocity_owner) + ((1 - IF)*x_velocity_neighbour);
                double y_face_velocity = (IF*y_velocity_owner) + ((1 - IF)*y_velocity_neighbour);     
                double z_face_velocity = (IF*z_velocity_owner) + ((1 - IF)*z_velocity_neighbour); 
                
                Vector face_velocity(x_face_velocity, y_face_velocity, z_face_velocity);
                face_flux_value = face_normal.dot_product(face_velocity);
                duplicate_face_flux_value = face_flux_value;
            }  


            if(list_of_faces[j] -> get_is_boundary_face() == true)
            { 
                int b_index = list_of_faces[j] -> get_boundary_index();
                Boundary b = list_of_boundaries[b_index];

                std::vector<double> bound_values = b.get_vector_dirichlet_boundary_value();

                Vector V_n(bound_values[0], bound_values[1], bound_values[2]);

                face_flux_value = face_normal.dot_product(V_n);
                duplicate_face_flux_value = face_flux_value;
            }  
              
            if(cell_index != owner_index)
            {
                duplicate_face_flux_value = duplicate_face_flux_value * -1;
            }  
            list_of_faces[j]->set_face_flux(duplicate_face_flux_value);
   
            sum_of_face_fluxes_through_cell += face_flux_value;
        }
        
        list_of_cells[i].set_sum_of_fluxes(sum_of_face_fluxes_through_cell);

   } 
}

std::vector<double> Vector_field::store_ap_coefficients()
{
    std::vector<double> diag_coeffs(b_vector_combined_ux.size());

    for(int i=0; i<b_vector_combined_ux.size(); i++)
    {
        diag_coeffs[i] = 1.0/a_matrix_combined(i,i);
    }

    return diag_coeffs;
}


  void Vector_field::correct_cell_centre_velocities(std::vector<Cell> &list_of_cells, std::vector<double> px_correction, std::vector<double> py_correction, std::vector<double> ap_coefficients)
  {
     for(int i=0; i<list_of_cells.size(); i++)
    {
        int cell_index = list_of_cells[i].get_cell_index();

        const_vector_field_values_x(i) = const_vector_field_values_x(i) - (ap_coefficients[i]* px_correction[i]);
        const_vector_field_values_y(i) = const_vector_field_values_y(i) - (ap_coefficients[i]* py_correction[i]);

    }

  }

void Vector_field::output_vector_matrix_coefficients_to_file(double total_cells)
{
     std::ofstream convection_coeffs_a("vel_convection_coeffs_a.txt");
     std::ofstream convection_coeffs_b_x("vel_convection_coeffs_b_x.txt");
     std::ofstream convection_coeffs_b_y("vel_convection_coeffs_b_y.txt");
     std::ofstream diffusion_coeffs_a("vel_diffusion_coeffs_a.txt");
     std::ofstream diffusion_coeffs_b_x("vel_diffusion_coeffs_b_x.txt");
     std::ofstream diffusion_coeffs_b_y("vel_diffusion_coeffs_b_y.txt");
     std::ofstream roc_coeffs_a("vel_roc_coeffs_a.txt");
     std::ofstream roc_coeffs_b("vel_roc_coeffs_b.txt");
     std::ofstream source_coeffs_b_x("vel_source_coeffs_b_x.txt");
     std::ofstream source_coeffs_b_y("vel_source_coeffs_b_y.txt");
     std::ofstream combined_coeffs_a("vel_combined_coeffs_a.txt");
     std::ofstream combined_coeffs_b_x("vel_combined_coeffs_b_x.txt");
     std::ofstream combined_coeffs_b_y("vel_combined_coeffs_b_y.txt");
    
     for(int i = 0; i<total_cells; i++)
     {
        for(int j=0; j<total_cells; j++)
              {
                if(j < total_cells - 1)
                {
                 convection_coeffs_a<<a_matrix_convection(i,j)<<"\t";
                 diffusion_coeffs_a<<a_matrix_diffusion(i,j)<<"\t";
                 roc_coeffs_a<<a_matrix_rate_of_change(i,j)<<"\t";
                 combined_coeffs_a<<a_matrix_combined(i,j)<<"\t";

                }

                else
                {
                 convection_coeffs_a<<a_matrix_convection(i,j)<<"\n";
                 diffusion_coeffs_a<<a_matrix_diffusion(i,j)<<"\n";
                 roc_coeffs_a<<a_matrix_rate_of_change(i,j)<<"\n";
                 combined_coeffs_a<<a_matrix_combined(i,j)<<"\n";

                }
              }
     }

          for(int k=0; k < total_cells; k++)
          {
                 convection_coeffs_b_x<<b_vector_convection_ux(k)<<"\n";
                 convection_coeffs_b_y<<b_vector_convection_uy(k)<<"\n";

                 diffusion_coeffs_b_x<<b_vector_diffusion_ux(k)<<"\n";
                 diffusion_coeffs_b_y<<b_vector_diffusion_uy(k)<<"\n";

                 roc_coeffs_b<<b_vector_rate_of_change(k)<<"\n";

                 source_coeffs_b_x<<b_vector_source_ux(k)<<"\n";
                 source_coeffs_b_y<<b_vector_source_ux(k)<<"\n";

                 combined_coeffs_b_x<<b_vector_combined_ux(k)<<"\n";
                 combined_coeffs_b_y<<b_vector_combined_uy(k)<<"\n";
          }    
     

   convection_coeffs_a.close();
   convection_coeffs_b_x.close();
   convection_coeffs_b_y.close();
   diffusion_coeffs_a.close();
   diffusion_coeffs_b_x.close();
   diffusion_coeffs_b_y.close();
   roc_coeffs_a.close();
   roc_coeffs_b.close();
   combined_coeffs_a.close();
   combined_coeffs_b_x.close();
   combined_coeffs_b_y.close();
   source_coeffs_b_x.close();
   source_coeffs_b_y.close();
}


void Vector_field::output_vector_field_to_file(std::vector<double> x_distance, std::vector<double> y_distance)
{
     std::ofstream vector_field_profiles_x("velocity_field_x.txt");

        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {  
                vector_field_profiles_x<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<const_vector_field_values_x(j + i*x_distance.size())<<std::endl;
              }
         }

        std::ofstream vector_field_profiles_y("velocity_field_y.txt");

        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                vector_field_profiles_y<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<const_vector_field_values_y(j + i*x_distance.size())<<std::endl;
              }
         }

        std::ofstream vector_field_profiles_xy("velocity_magnitude.txt");

        for(int i=0; i<y_distance.size(); i++)
         {
              for(int j=0; j<x_distance.size(); j++)
              {
                 double vx = const_vector_field_values_x(j + i*x_distance.size());
                 double vy = const_vector_field_values_y(j + i*x_distance.size());
                 double v_mag = sqrt(vx*vx +  vy*vy);

                vector_field_profiles_xy<<x_distance[j]<<"\t"<<y_distance[i]<<"\t"<<v_mag<<std::endl;
              }
         }


        std::ofstream vx_horizontal("vx_horizontal.txt");

        for(int j=0; j<x_distance.size(); j++)
         {
            vx_horizontal<<x_distance[j]<<"\t"<< const_vector_field_values_x(((y_distance.size()/2) - 1.0)*x_distance.size() + j)<<std::endl;
         }



        std::ofstream vx_vertical("vx_vertical.txt");
        for(int i=0; i<y_distance.size(); i++)
         {
            vx_vertical<<const_vector_field_values_x((x_distance.size()/2) + i*x_distance.size())<<"\t"<<y_distance[i]<<std::endl;
         }


        std::ofstream vy_horizontal("vy_horizontal.txt");

        for(int j=0; j<x_distance.size(); j++)
         {
            vy_horizontal<<x_distance[j] <<"\t" <<const_vector_field_values_y(((y_distance.size()/2) - 1.0)*x_distance.size() + j)<<std::endl;
         }

        std::ofstream vy_vertical("vy_vertical.txt") ;

        for(int i=0; i<y_distance.size(); i++)
         {
            vy_vertical<<y_distance[i]<<"\t"<< const_vector_field_values_y((x_distance.size()/2) + i*x_distance.size())<<std::endl;
         }

 
        vector_field_profiles_x.close();
        vector_field_profiles_y.close();
        vector_field_profiles_xy.close();
        vx_horizontal.close();
        vx_vertical.close();
        vy_horizontal.close();
        vy_vertical.close();
} 

void Vector_field::plot_convergence_initial_x(std::ofstream &x_velocity_convergence_initial, int iteration_no)
{   int index;
    Eigen::VectorXd diff = (const_vector_field_values_x - const_vector_old_field_values_x);

    for(int i=0; i<const_vector_field_values_x.size(); i++)
    {
        diff(i) = abs(diff(i));
    }

     double max = diff(0);
    
    for(int i=1; i<const_vector_field_values_x.size(); i++)
    {   double next_term = diff(i);
        double last_term = diff(i-1);
        
        if(next_term > last_term)
        {
         max = diff(i);
         index = i;
        }
    }

     x_velocity_convergence_initial<<iteration_no<<"\t"<<max<<"\t"<<std::endl;
}

void Vector_field::plot_convergence_initial_y(std::ofstream &y_velocity_convergence_initial, int iteration_no)
{
    Eigen::VectorXd difference = (const_vector_field_values_y - const_vector_old_field_values_y);

    for(int i=0; i<const_vector_field_values_y.size(); i++)
    {
        difference(i) = abs(difference(i));
    }

     double max = difference(0);
    
    for(int i=1; i<const_vector_field_values_y.size(); i++)
    {   double next_term = difference(i);
        double last_term = difference(i-1);
        
        if(next_term > last_term)
        {
         max = difference(i);
        }
    }

     y_velocity_convergence_initial<<iteration_no<<"\t"<<max<<std::endl;

}


void Vector_field::plot_convergence_final_x(std::ofstream &x_velocity_convergence_final, int iteration_no)
{
       Eigen::VectorXd difference = (const_vector_field_values_x - const_vector_old_field_values_x);

    for(int i=0; i<const_vector_field_values_x.size(); i++)
    {
        difference(i) = abs(difference(i));
    }

     double max = difference(0);
    
    for(int i=1; i<const_vector_field_values_x.size(); i++)
    {   double next_term = difference(i);
        double last_term = difference(i-1);
        
        if(next_term > last_term)
        {
         max = difference(i);
        }
    }

     x_velocity_convergence_final<<iteration_no<<"\t"<<max<<std::endl;
}


void Vector_field::plot_convergence_final_y(std::ofstream &y_velocity_convergence_final, int iteration_no)
{
      Eigen::VectorXd difference = (const_vector_field_values_y - const_vector_old_field_values_y);

    for(int i=0; i<const_vector_field_values_y.size(); i++)
    {
        difference(i) = abs(difference(i));
    }

     double max = difference(0);
    
    for(int i=1; i<const_vector_field_values_y.size(); i++)
    {   double next_term = difference(i);
        double last_term = difference(i-1);
        
        if(next_term > last_term)
        {
         max = difference(i);
        }
    }

     y_velocity_convergence_final<<iteration_no<<"\t"<<max<<std::endl;
}


    