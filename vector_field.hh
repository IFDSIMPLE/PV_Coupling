#ifndef VECTOR_FIELD_HH
#define VECTOR_FIELD_HH

#include<Eigen/Dense>
#include<vector>
#include "cell.hh"
#include "boundary.hh"
#include "const_scalar_field.hh"
#include "const_vector_field.hh"

class Vector_field : public Const_vector_field
{
  public:
  
  Eigen::MatrixXd a_matrix_rate_of_change;
  Eigen::VectorXd b_vector_rate_of_change;

  Eigen::VectorXd b_vector_source_ux;
  Eigen::VectorXd b_vector_source_uy;

  Eigen::MatrixXd a_matrix_diffusion;
  Eigen::VectorXd b_vector_diffusion_ux;
  Eigen::VectorXd b_vector_diffusion_uy;

  Eigen::MatrixXd a_matrix_convection;
  Eigen::VectorXd b_vector_convection_ux;
  Eigen::VectorXd b_vector_convection_uy;

  Eigen::MatrixXd a_matrix_combined;
  Eigen::VectorXd b_vector_combined_ux;
  Eigen::VectorXd b_vector_combined_uy;

  Eigen::VectorXd initial_residuals_ux;
  Eigen::VectorXd initial_residuals_uy;
  Eigen::VectorXd final_residuals_ux;
  Eigen::VectorXd final_residuals_uy;

  Vector_field(int a);

  void update_a_matrix_rate_of_change(double, int);
  void update_b_vector_rate_of_change(double, int);

  void update_b_vector_source_ux(double, int);
  void update_b_vector_source_uy(double, int);

  void update_a_matrix_diffusion(double, int, int);
  void update_b_vector_diffusion_ux(double, int);
  void update_b_vector_diffusion_uy(double, int);

  void update_a_matrix_convection(double, int, int);
  void update_b_vector_convection_ux(double, int);
  void update_b_vector_convection_uy(double, int);

  void combine_a_matrices();
  void combine_b_matrices();

  void reset_rate_of_change();
  void reset_source();
  void reset_diffusion();
  void reset_convection();
  void reset_combined();

  void compute_rate_of_change_matrix(std::vector<Cell>, double);
  void compute_source_matrix(std::vector<Cell>);
  void compute_diffusion_matrix(std::vector<Cell>, Const_scalar_field obj_1, std::vector<Boundary>);
  void compute_convection_matrix(std::vector<Cell>, std::vector<Boundary>, Vector_field obj_2);
  void under_relaxation(std::vector<Cell>, double);
  void solve_matrices();

  void calculate_initial_residuals(std::vector<double>, std::vector<double>, int, std::ofstream &, std::ofstream &);
  void calculate_final_residuals(std::vector<double>, std::vector<double>, int, std::ofstream &, std::ofstream &);

  std::vector<double> store_ap_coefficients();

  void output_vector_matrix_coefficients_to_file(double);
  void output_vector_field_to_file(std::vector<double>, std::vector<double>);

  void correct_cell_fluxes(std::vector<Cell>, std::vector<double>);
  void correct_cell_centre_velocities(std::vector<Cell>&, std::vector<double> , std::vector<double> , std::vector<double>);

  //std::vector<double> set_face_and_cell_fluxes(std::vector<Cell> &);
  void set_face_and_cell_fluxes(std::vector<Cell> &, std::vector<Boundary>);

  void plot_convergence_initial_x(std::ofstream &, int);
  void plot_convergence_initial_y(std::ofstream &, int);
  void plot_convergence_final_x(std::ofstream &, int);
  void plot_convergence_final_y(std::ofstream &, int);

};

#endif