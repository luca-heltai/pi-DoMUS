/*! \addtogroup equations
 * @{
 */

/**
 * This interface solves ALE Navier Stokes Equation:
 *
 */

#ifndef _pidomus_ALE_navier_stokes_h_
#define _pidomus_ALE_navier_stokes_h_

#include "pde_system_interface.h"
#include "../boundary_values.h"

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
 #include <deal.II/grid/grid_tools.h>

#include <deal2lkit/parsed_mapped_functions.h>
#include <deal2lkit/parsed_preconditioner/amg.h>
#include <deal2lkit/parsed_preconditioner/jacobi.h>

////////////////////////////////////////////////////////////////////////////////
/// ALE Navier Stokes interface:

template <int dim, int spacedim=dim, typename LAC=LATrilinos>
class ALENavierStokes
  :
  public PDESystemInterface<dim,spacedim,ALENavierStokes<dim,spacedim,LAC>, LAC>
{

public:
  ~ALENavierStokes () {};
  ALENavierStokes ();

  void declare_parameters (ParameterHandler &prm);
  void parse_parameters_call_back ();

  template <typename EnergyType, typename ResidualType>
  void
  energies_and_residuals(
    const typename DoFHandler<dim,spacedim>::active_cell_iterator &cell,
    FEValuesCache<dim,spacedim> &scratch,
    std::vector<EnergyType> &energies,
    std::vector<std::vector<ResidualType>> &residuals,
    bool compute_only_system_terms) const;

  void
  compute_system_operators(
    const std::vector<shared_ptr<LATrilinos::BlockMatrix>>,
    LinearOperator<LATrilinos::VectorType> &,
    LinearOperator<LATrilinos::VectorType> &,
    LinearOperator<LATrilinos::VectorType> &) const;
  
  virtual void connect_to_signals() const
  {
    auto &signals = this->get_signals();

    signals.postprocess_newly_created_triangulation.connect(
      [&,this](Triangulation<dim> &tria)
    {
      int index = (dim==2)? 1:0;
      Tensor<1, dim> shift_vec;
      shift_vec[index] = -1.318;
      GridTools::shift(shift_vec, tria);
    }
    );

    signals.begin_make_grid_fe.connect(
      [&,this]()
    {
      adaptive_preconditioners_on = this->get_adaptive_preconditioners();
      max_iterations_adaptive = this->get_max_iterations_adaptive();
    }
    );

    signals.end_solve_jacobian_system.connect(
      [&,this]()
    {
      double time = this->get_current_time();
      if (time > 0.005)
        iterations_last_step = this->get_solver_control()->last_step(); 
      else 
        iterations_last_step = 0;     
    }
    );

    signals.update_constraint_matrices.connect(
      [&,this](std::vector<std::shared_ptr<dealii::ConstraintMatrix> > &constraints, ConstraintMatrix &constraints_dot) 
    {
      auto &dof=this->get_dof_handler();
      auto &fe =this->get_fe();

      FEValuesExtractors::Vector displacements (0);
      FEValuesExtractors::Vector velocities (dim);
      ComponentMask displacement_mask = fe.component_mask (displacements);
      ComponentMask velocity_mask = fe.component_mask (velocities);

      // in 3D:
      // displacement_mask = [1 1 1 0 0 0 0]
      // velocity_mask     = [0 0 0 1 1 1 0]
      
      double timestep = this->get_current_time();
      double dt = this->get_timestep();

      //std::cout << "timestep " << timestep << std::endl;
      //std::cout << "dt " << dt << std::endl;
      // if timestep == nan
      if (timestep != timestep)
      {
        timestep = 0;
      }
      if (dt != dt) {
          dt = 1;
      }


      // dirichlet BC for d 
      if (dim==2)
      {
        // bottom face
        VectorTools::interpolate_boundary_values (dof,
                                                  2,
                                                  BoundaryValues<dim>(2, timestep, dt, false, 2),
                                                  *constraints[0],
                                                  displacement_mask);
        // left hull
        VectorTools::interpolate_boundary_values (dof,
                                                  0,
                                                  BoundaryValues<dim>(0, timestep, dt, true, 4),
                                                  *constraints[0],
                                                  displacement_mask);
        // right hull
        VectorTools::interpolate_boundary_values (dof,
                                                  1,
                                                  BoundaryValues<dim>(1, timestep, dt, true, 4),
                                                  *constraints[0],
                                                  displacement_mask);
        // top face
        VectorTools::interpolate_boundary_values (dof,
                                                  3,
                                                  BoundaryValues<dim>(3),
                                                  *constraints[0],
                                                  displacement_mask);
      }
      else
      {
        // bottom face
        VectorTools::interpolate_boundary_values (dof,
                                                  1,
                                                  BoundaryValues<dim>(1, timestep, dt, false, 2),
                                                  *constraints[0],
                                                  displacement_mask);
        // hull
        VectorTools::interpolate_boundary_values (dof,
                                                  0,
                                                  BoundaryValues<dim>(0, timestep, dt, true, 4),
                                                  *constraints[0],
                                                  displacement_mask);
        // top face
        VectorTools::interpolate_boundary_values (dof,
                                                  2,
                                                  BoundaryValues<dim>(2),
                                                  *constraints[0],
                                                  displacement_mask);
      }

      int n_faces = (dim==2)? 3 : 2;
      Point<3> degrees (2, 4, 4); 
      Point<3> bools (false, true, true);

      if(timestep < 0.005)  // 0.005 is the time of one heart interval
      {
        // time derivatives of dirichlet BC for d
        for (int i = 0; i < n_faces; ++i)
        { 
          VectorTools::interpolate_boundary_values (dof,
                                                    i,
                                                    ZeroFunction<dim>(2*dim+1),
                                                    constraints_dot,
                                                    displacement_mask);
        }
      }
      else
      {
        // time derivatives of dirichlet BC for d
        for (int j = 0; j < n_faces; ++j)
        {
          VectorTools::interpolate_boundary_values (dof,
                                                    j,
                                                    BoundaryValues<dim>(j, timestep, dt, bools(n_faces-1-j), degrees(n_faces-1-j), true),
                                                    constraints_dot,
                                                    displacement_mask);
        }
        // BC for u
        for (int j = 0; j < n_faces; ++j)
        {
          VectorTools::interpolate_boundary_values (dof,
                                                    j,
                                                    BoundaryValues<dim>(j, timestep, dt, bools(n_faces-1-j), degrees(n_faces-1-j), true),
                                                    *constraints[0],
                                                    velocity_mask);
        }
      }
    }
    );
  } 
  

  //void
  //set_matrix_couplings(std::vector<std::string> &couplings) const;

private:

// Physical parameter
  double nu;
  double rho;
  mutable unsigned int iterations_last_step;
  mutable unsigned int max_iterations_adaptive;
  mutable bool adaptive_preconditioners_on;

  bool Mp_use_inverse_operator;
  bool AMG_u_use_inverse_operator;
  bool AMG_d_use_inverse_operator;

  /**
   * AMG preconditioner for velocity-velocity matrix.
   */
  mutable ParsedAMGPreconditioner AMG_u;

  /**
   * AMG preconditioner for the pressure stifness matrix.
   */
  mutable ParsedAMGPreconditioner AMG_d;

  /**
   * Jacobi preconditioner for the pressure mass matrix.
   */
  mutable ParsedJacobiPreconditioner jac_M;
};


template <int dim, int spacedim, typename LAC>
ALENavierStokes<dim,spacedim, LAC>::
ALENavierStokes()
  :
  PDESystemInterface<dim,spacedim,ALENavierStokes<dim,spacedim,LAC>, LAC>(
    "ALE Navier Stokes Interface",
    dim+dim+1,
    2,
    "FESystem[FE_Q(2)^d-FE_Q(2)^d-FE_Q(1)]",
    (dim==2)? "d,d,u,u,p" : "d,d,d,u,u,u,p",
    "1,1,0"),
  AMG_u("AMG for u"),
  AMG_d("AMG for d"),
  jac_M("Jacobi for M")
{
  this->init();
}

template <int dim, int spacedim, typename LAC>
void ALENavierStokes<dim,spacedim,LAC>::
declare_parameters (ParameterHandler &prm)
{
  PDESystemInterface<dim,spacedim, ALENavierStokes<dim,spacedim,LAC>,LAC>::
  declare_parameters(prm);

  this->add_parameter(prm, &nu,
                      "nu [Pa s]", "1.0",
                      Patterns::Double(0.0),
                      "Viscosity");

  this->add_parameter(prm, &rho,
                      "rho [kg m^-d]", "1.0",
                      Patterns::Double(0.0),
                      "Density");

  this->add_parameter(prm, &Mp_use_inverse_operator,
                      "Invert Mp using inverse operator", "false",
                      Patterns::Bool(),
                      "Invert Mp usign inverse operator");

  this->add_parameter(prm, &AMG_d_use_inverse_operator,
                      "AMG d - use inverse operator", "false",
                      Patterns::Bool(),
                      "Enable the use of inverse operator for AMG d");

  this->add_parameter(prm, &AMG_u_use_inverse_operator,
                      "AMG u - use inverse operator", "false",
                      Patterns::Bool(),
                      "Enable the use of inverse operator for AMG u");
}

template <int dim, int spacedim, typename LAC>
void ALENavierStokes<dim,spacedim,LAC>::
parse_parameters_call_back ()
{}

//template <int dim, int spacedim, typename LAC>
//void ALENavierStokes<dim,spacedim,LAC>::
//set_matrix_couplings(std::vector<std::string> &couplings) const
//{
//  couplings[0] = "1,0,0; 0,1,1; 0,1,0"; // TODO: Select only not null entries
//  couplings[1] = "0,0,0; 0,0,0; 0,0,1";
//}

template <int dim, int spacedim, typename LAC>
template <typename EnergyType, typename ResidualType>
void
ALENavierStokes<dim,spacedim,LAC>::
energies_and_residuals(const typename DoFHandler<dim,spacedim>::active_cell_iterator &cell,
                       FEValuesCache<dim,spacedim> &fe_cache,
                       std::vector<EnergyType> &,
                       std::vector<std::vector<ResidualType> > &residual,
                       bool compute_only_system_terms) const
{
  const FEValuesExtractors::Vector displacement(0);
  const FEValuesExtractors::Vector velocity(dim);
  const FEValuesExtractors::Scalar pressure(2*dim);

  ResidualType et = 0;
  double dummy = 0.0;

  double h = cell->diameter();

  for (unsigned int face=0; face < GeometryInfo<dim>::faces_per_cell; ++face)
  {
    if (cell->face(face)->at_boundary())
    {
      this->reinit(et, cell, face, fe_cache);
      
      // Displacement:
      auto &d_dot_ = fe_cache.get_values( "solution_dot", "d_dot", displacement, et);

      // Velocity:
      auto &u_ = fe_cache.get_values("solution", "grad_u", velocity, et);

      // Pressure:
      // auto &grad_p_ = fe_cache.get_gradients("solution", "p", pressure, et);

      auto &fev = fe_cache.get_current_fe_values();
      auto &q_points = fe_cache.get_quadrature_points();
      auto &JxW = fe_cache.get_JxW_values();

      for (unsigned int q=0; q<q_points.size(); ++q)
      {
        // Displacement:
        const Tensor<1, dim, ResidualType> &d_dot = d_dot_[q];

        // Velocity:
        const Tensor<1, dim, ResidualType> &u = u_[q];

        for (unsigned int i=0; i<residual[0].size(); ++i)
        {
          // Test functions:
          // auto d_test = fev[displacement].value(i,q);
          auto u_test = fev[velocity].value(i,q);

          residual[0][i] += (1./h)*( (u - d_dot) * u_test )*JxW[q];
        }
      } // end loop over quadrature points
      break;
    }
  } // end loop over faces

  this->reinit (et, cell, fe_cache);

  // displacement:
  auto &grad_ds = fe_cache.get_gradients( "solution", "grad_d", displacement, et);
  auto &Fs = fe_cache.get_deformation_gradients( "solution", "Fd", displacement, et);
  auto &ds_dot = fe_cache.get_values( "solution_dot", "d_dot", displacement, et);

  // velocity:
  auto &us = fe_cache.get_values( "solution", "u", velocity, et);
  auto &grad_us = fe_cache.get_gradients( "solution", "grad_u", velocity, et);
  auto &div_us = fe_cache.get_divergences( "solution", "div_u", velocity, et);
  auto &sym_grad_us = fe_cache.get_symmetric_gradients( "solution", "u", velocity, et);
  auto &us_dot = fe_cache.get_values( "solution_dot", "u_dot", velocity, et);

  // Previous time step solution:
  auto &u_olds = fe_cache.get_values("explicit_solution","u",velocity,dummy);  

  // pressure:
  auto &ps = fe_cache.get_values( "solution", "p", pressure, et);

  // Jacobian:
  auto &JxW = fe_cache.get_JxW_values();

  const unsigned int n_quad_points = us.size();

  auto &fev = fe_cache.get_current_fe_values();

  for (unsigned int quad=0; quad<n_quad_points; ++quad)
  {
    // velocity:
    const ResidualType &div_u = div_us[quad];
    const Tensor<1, dim, ResidualType> &u_dot = us_dot[quad];
    const Tensor<2, dim, ResidualType> &grad_u = grad_us[quad];
    const Tensor <2, dim, ResidualType> &sym_grad_u = sym_grad_us[quad];
    
    // displacement
    const Tensor<1, dim, ResidualType> &d_dot = ds_dot[quad];
    const Tensor<2, dim, ResidualType> &grad_d = grad_ds[quad];

    // deformation gradient
    const Tensor <2, dim, ResidualType> &F = Fs[quad];
    ResidualType J = determinant(F);
    const Tensor <2, dim, ResidualType> &F_inv = invert(F);
    const Tensor <2, dim, ResidualType> &Ft_inv = transpose(F_inv);

    // Previous time step solution:
    const Tensor<1, dim, ResidualType> &u_old = u_olds[quad];

    // pressure:
    const ResidualType &p = ps[quad];

    // jacobian of ALE transformation
    auto J_ale = J; 
    // auto div_u_ale = (J_ale * (F_inv * u) );
    
    // pressure * identity matrix
    Tensor <2, dim, ResidualType> p_Id;
    for (unsigned int i = 0; i<dim; ++i)
      p_Id[i][i] = p;

    ResidualType my_rho = rho;
    const Tensor <2, dim, ResidualType> sigma = -p_Id + my_rho*(nu*sym_grad_u*F_inv + (Ft_inv * transpose(nu*sym_grad_u) ) );

    for (unsigned int i=0; i<residual[0].size(); ++i)
    {
      // test functions:

      // velocity:
      auto u_test = fev[velocity].value(i,quad);
      auto grad_u_test = fev[velocity].gradient(i,quad);
      auto div_u_test = fev[velocity].divergence(i,quad);

      // displacement:
      auto grad_d_test = fev[displacement].gradient(i,quad);

      // pressure:
      auto p_test = fev[pressure].value(i,quad);

      residual[1][i] += ( (1./nu)*p*p_test )*JxW[quad];

      residual[0][i] +=
        (
          // time derivative term
          rho*scalar_product( u_dot * J_ale , u_test )

          + rho*scalar_product( grad_u * ( F_inv * ( u_old - d_dot ) ) * J_ale , u_test )

          + scalar_product( J_ale * sigma * Ft_inv, grad_u_test )
          // divergence free constriant
          - div_u * p_test
          // pressure term
          - p * div_u_test
          // Impose harmony of d and u=d_dot
          + scalar_product( grad_d , grad_d_test )
        )*JxW[quad];

        }
    }

  (void)compute_only_system_terms;
}

template <int dim, int spacedim, typename LAC>
void
ALENavierStokes<dim,spacedim,LAC>::compute_system_operators(
  const std::vector<shared_ptr<LATrilinos::BlockMatrix>> matrices,
  LinearOperator<LATrilinos::VectorType> &system_op,
  LinearOperator<LATrilinos::VectorType> &prec_op,
  LinearOperator<LATrilinos::VectorType> &) const
{
  typedef LATrilinos::VectorType::BlockType BVEC;
  typedef LATrilinos::VectorType VEC;

// Preconditioners:
  const DoFHandler<dim,spacedim> &dh = this->get_dof_handler();
  const ParsedFiniteElement<dim,spacedim> fe = this->pfe;
  
// TODO 1: DONE find a way to check if local refinement is turned off in "compute_system_operators()" 

// TODO 2: DONE check convergence time of last step
// get iterations of last step, if iterations > 15: refine

  static int counter = 0;

  if (adaptive_preconditioners_on == true)
  {
    double time = this->get_current_time();
    // only refine at the beginning and if convergence is slow
    if (time <= 0.006 || iterations_last_step > max_iterations_adaptive)
    {
      AMG_d.initialize_preconditioner<dim, spacedim>( matrices[0]->block(0,0), fe, dh);
      AMG_u.initialize_preconditioner<dim, spacedim>( matrices[0]->block(1,1), fe, dh);
      jac_M.initialize_preconditioner<>(matrices[1]->block(2,2));
      counter++;
      std::cout << "gettig reinitialized!" << std::endl;
    }
  } 
  else
  {
      AMG_d.initialize_preconditioner<dim, spacedim>( matrices[0]->block(0,0), fe, dh);
      AMG_u.initialize_preconditioner<dim, spacedim>( matrices[0]->block(1,1), fe, dh);
      jac_M.initialize_preconditioner<>(matrices[1]->block(2,2));
      //std::cout << "gettig reinitialized always!" << std::endl;
  }
////////////////////////////////////////////////////////////////////////////
// SYSTEM MATRIX:

  std::array<std::array<LinearOperator< BVEC >, 3>, 3> S;
  for (unsigned int i = 0; i<3; ++i)
    for (unsigned int j = 0; j<3; ++j)
      S[i][j] = linear_operator< BVEC >(matrices[0]->block(i,j) );
  system_op = BlockLinearOperator< VEC >(S);

////////////////////////////////////////////////////////////////////////////
// PRECONDITIONER MATRIX:

  std::array<std::array<LinearOperator<TrilinosWrappers::MPI::Vector>, 3>, 3> P;
  for (unsigned int i = 0; i<3; ++i)
    for (unsigned int j = 0; j<3; ++j)
      P[i][j] = linear_operator<BVEC>( matrices[0]->block(i,j) );

  static ReductionControl solver_control_pre(5000, 1e-8);
  static SolverCG<BVEC> solver_CG(solver_control_pre);
  static SolverGMRES<BVEC> solver_GMRES(solver_control_pre);

  for (unsigned int i = 0; i<3; ++i)
    for (unsigned int j = 0; j<3; ++j)
      if (i!=j)
        P[i][j] = null_operator< TrilinosWrappers::MPI::Vector >(P[i][j]);

  auto A = linear_operator< BVEC>(matrices[0]->block(1,1) );
  auto B = linear_operator< BVEC>(matrices[0]->block(2,1) );
  auto Bt = transpose_operator< >(B);

  LinearOperator<BVEC> A_inv;
  if (AMG_u_use_inverse_operator)
    {
      A_inv = inverse_operator( S[1][1],
                                solver_GMRES,
                                AMG_u);
    }
  else
    {
      A_inv = linear_operator<BVEC>(matrices[0]->block(1,1),
                                    AMG_u);
    }

  auto Mp = linear_operator< TrilinosWrappers::MPI::Vector >( matrices[1]->block(2,2) );

  LinearOperator<BVEC> Mp_inv;
  if (Mp_use_inverse_operator)
    {
      Mp_inv = inverse_operator(Mp,
                                solver_GMRES,
                                jac_M);
    }
  else
    {
      Mp_inv = linear_operator<BVEC>(matrices[1]->block(2,2),
                                     jac_M);
    }

  auto Schur_inv = nu * Mp_inv;

  if (AMG_d_use_inverse_operator)
    {
      P[0][0] = inverse_operator( S[0][0],
                                  solver_CG,
                                  AMG_d);
    }
  else
    {
      P[0][0] = linear_operator<BVEC>(matrices[0]->block(0,0),
                                      AMG_d);
    }

  P[1][1] = A_inv;
  P[1][2] = A_inv * Bt * Schur_inv;
  P[2][1] = null_operator(B);
  P[2][2] = -1 * Schur_inv;


  prec_op = BlockLinearOperator< VEC >(P);
}

#endif

/*! @} */
