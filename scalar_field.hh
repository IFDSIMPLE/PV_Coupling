#ifndef SCALAR_FIELD_HH
#define SCALAR_FIELD_HH

#include<Eigen/Dense>
#include<vector>
#include<cmath>

#include "const_scalar_field.hh"
#include "const_vector_field.hh"
#include "cell.hh"
#include "boundary.hh"

class Scalar_field : public Const_scalar_field
{
  public:
   
  Eigen::MatrixXd a_matrix_rate_of_change;
  Eigen::VectorXd b_vector_rate_of_change;

  Eigen::VectorXd b_vector_source;

  Eigen::MatrixXd a_matrix_diffusion;
  Eigen::VectorXd b_vector_diffusion;

  Eigen::MatrixXd a_matrix_convection;
  Eigen::VectorXd b_vector_convection;

  Eigen::MatrixXd a_matrix_combined;
  Eigen::VectorXd b_vector_combined;

  Eigen::VectorXd initial_residuals_p;
  Eigen::VectorXd final_residuals_p;
  

  Scalar_field(int a);

  void update_a_matrix_rate_of_change(double, int);
  void update_b_vector_rate_of_change(double, int);

  void update_b_vector_source(double, int);

  void update_a_matrix_diffusion(double, int, int);
  void update_b_vector_diffusion(double, int);

  void update_a_matrix_convection(double, int, int);
  void update_b_vector_convection(double, int);

  void combine_a_and_b_matrices();

  void reset_rate_of_change();
  void reset_source();
  void reset_diffusion();
  void reset_convection();
  void reset_combined();

  void compute_rate_of_change_matrix(std::vector<Cell> , double);
  void compute_source_matrix(std::vector<Cell>);
  void compute_diffusion_matrix(std::vector<Cell>, std::vector<double>, std::vector<Boundary>, int);
  void compute_convection_matrix(std::vector<Cell>, std::vector<Boundary>, Const_vector_field obj_2);
  void combine_and_solve_matrices(std::vector<Cell>);

  void pressure_under_relax();

  Const_vector_field compute_gauss_gradient(const std::vector<Cell> &listOfCells) const;

  void calculate_initial_residuals_p(std::vector<double>, std::vector<double>, int, std::ofstream &);
  void calculate_final_residuals_p(std::vector<double>, std::vector<double>, int,  std::ofstream &);

  void compute_flux_correction(std::vector<Cell>& , std::vector<double>, std::vector<Boundary>);
  void compute_velocity_correction_terms(std::vector<Cell>&, std::vector<double>&, std::vector<double>&, std::vector<Boundary>);

  void output_scalar_matrix_coefficients_to_file(double);
  void output_scalar_field_to_file(std::vector<double>, std::vector<double>, std::vector<Cell>);
  void plot_convergence_initial(std::ofstream &, int);
};

#endif