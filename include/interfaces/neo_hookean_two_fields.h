#ifndef _neo_hookean_two_fields_interface_h_
#define _neo_hookean_two_fields_interface_h_

#include "conservative_interface.h"
#include <deal2lkit/parsed_function.h>


#include <deal.II/fe/fe_values.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_precondition.h>

#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/block_linear_operator.h>
#include <deal.II/lac/packaged_operation.h>

#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>


template <int dim, int spacedim>
class NeoHookeanTwoFieldsInterface : public ConservativeInterface<dim,spacedim,dim+1, NeoHookeanTwoFieldsInterface<dim,spacedim> >
{
public:

  typedef FEValuesCache<dim,spacedim> Scratch;
  typedef Assembly::CopyData::piDoMUSPreconditioner<dim,dim> CopyPreconditioner;
  typedef Assembly::CopyData::piDoMUSSystem<dim,dim> CopySystem;
  typedef TrilinosWrappers::MPI::BlockVector VEC;

  /* specific and useful functions for this very problem */
  NeoHookeanTwoFieldsInterface();

  virtual void declare_parameters (ParameterHandler &prm);
  virtual void parse_parameters_call_back ();

  /* these functions MUST have the follwowing names
   *  because they are called by the ConservativeInterface class
   */
  template<typename Number>
  void preconditioner_energy(const typename DoFHandler<dim,spacedim>::active_cell_iterator &,
                             Scratch &,
                             CopyPreconditioner &,
                             Number &energy) const;

  template<typename Number>
  void system_energy(const typename DoFHandler<dim,spacedim>::active_cell_iterator &,
                     Scratch &,
                     CopySystem &,
                     Number &energy) const;

  virtual void compute_system_operators(const DoFHandler<dim,spacedim> &,
                                        const TrilinosWrappers::BlockSparseMatrix &,
                                        const TrilinosWrappers::BlockSparseMatrix &,
                                        BlockLinearOperator<VEC> &,
                                        BlockLinearOperator<VEC> &) const;

private:
  double E;
  double nu;
  double mu;
  double lambda;

  mutable shared_ptr<TrilinosWrappers::PreconditionAMG>    Amg_preconditioner;
  mutable shared_ptr<TrilinosWrappers::PreconditionJacobi> Mp_preconditioner;
  mutable shared_ptr<TrilinosWrappers::PreconditionJacobi> T_preconditioner;

};

template <int dim, int spacedim>
NeoHookeanTwoFieldsInterface<dim,spacedim>::NeoHookeanTwoFieldsInterface() :
  ConservativeInterface<dim,spacedim,dim+1,NeoHookeanTwoFieldsInterface<dim,spacedim> >("NeoHookean Interface",
      "FESystem[FE_Q(2)^d-FE_DGP(1)]",
      "u,u,u,p", "1,1; 1,0", "1,0; 0,1",
      "1,0")
{};



template <int dim, int spacedim>
template<typename Number>
void NeoHookeanTwoFieldsInterface<dim,spacedim>::preconditioner_energy(const typename DoFHandler<dim,spacedim>::active_cell_iterator &cell,
    Scratch &fe_cache,
    CopyPreconditioner &data,
    Number &energy) const
{
  Number alpha = this->alpha;
  this->reinit(alpha, cell, fe_cache);

  auto &JxW = fe_cache.get_JxW_values();

  const FEValuesExtractors::Vector displacement(0);
  const FEValuesExtractors::Scalar pressure(dim);
  auto &ps = fe_cache.get_values("solution","p", pressure, alpha);
  auto &grad_us = fe_cache.get_gradients("solution","grad_u",displacement, alpha);

  const unsigned int n_q_points = ps.size();

  energy = 0;
  for (unsigned int q=0; q<n_q_points; ++q)
    {
      const Number &p = ps[q];
      const Tensor <2, dim, Number> &grad_u = grad_us[q];

      energy += (scalar_product(grad_u,grad_u) +
                 0.5*p*p)*JxW[q];
    }
}

template <int dim, int spacedim>
template<typename Number>
void NeoHookeanTwoFieldsInterface<dim,spacedim>::system_energy(const typename DoFHandler<dim,spacedim>::active_cell_iterator &cell,
    Scratch &fe_cache,
    CopySystem &data,
    Number &energy) const
{
  Number alpha = this->alpha;
  this->reinit(alpha, cell, fe_cache);

  const FEValuesExtractors::Vector displacement(0);
  auto &us = fe_cache.get_values("solution", "u", displacement, alpha);
  auto &us_dot = fe_cache.get_values("solution_dot", "u_dot", displacement, alpha);
  auto &Fs = fe_cache.get_deformation_gradients("solution", "Fu", displacement, alpha);

  const FEValuesExtractors::Scalar pressure(dim);
  auto &ps = fe_cache.get_values("solution","p", pressure, alpha);

  auto &JxW = fe_cache.get_JxW_values();

  const unsigned int n_q_points = ps.size();

  energy = 0;
  for (unsigned int q=0; q<n_q_points; ++q)
    {
      Tensor <1, dim, double> B;

      const Tensor <1, dim, Number> &u = us[q];
      const Tensor <1, dim, Number> &u_dot = us_dot[q];
      const Number &p = ps[q];
      const Tensor <2, dim, Number> &F = Fs[q];
      const Tensor<2, dim, Number> C = transpose(F)*F;

      Number Ic = trace(C);
      Number J = determinant(F);

      Number psi = (mu/2.)*(Ic-dim) +p*(J-1.);
      energy += (u*u_dot + psi)*JxW[q];
    }
}


template <int dim, int spacedim>
void NeoHookeanTwoFieldsInterface<dim,spacedim>::declare_parameters (ParameterHandler &prm)
{
  ConservativeInterface<dim,spacedim,dim+1, NeoHookeanTwoFieldsInterface<dim,spacedim> >::declare_parameters(prm);
  this->add_parameter(prm, &E, "Young's modulus", "10.0", Patterns::Double(0.0));
  this->add_parameter(prm, &nu, "Poisson's ratio", "0.3", Patterns::Double(0.0));
}

template <int dim, int spacedim>
void NeoHookeanTwoFieldsInterface<dim,spacedim>::parse_parameters_call_back ()
{
  ConservativeInterface<dim,spacedim,dim+1, NeoHookeanTwoFieldsInterface<dim,spacedim> >::parse_parameters_call_back();
  mu = E/(2.0*(1.+nu));
  lambda = (E *nu)/((1.+nu)*(1.-2.*nu));
}

template <int dim, int spacedim>
void
NeoHookeanTwoFieldsInterface<dim,spacedim>::compute_system_operators(const DoFHandler<dim,spacedim> &dh,
    const TrilinosWrappers::BlockSparseMatrix &matrix,
    const TrilinosWrappers::BlockSparseMatrix &preconditioner_matrix,
    BlockLinearOperator<VEC> &system_op,
    BlockLinearOperator<VEC> &prec_op) const
{

  std::vector<std::vector<bool> > constant_modes;
  FEValuesExtractors::Vector displacement(0);
  DoFTools::extract_constant_modes (dh, dh.get_fe().component_mask(displacement),
                                    constant_modes);

  Mp_preconditioner.reset  (new TrilinosWrappers::PreconditionJacobi());
  Amg_preconditioner.reset (new TrilinosWrappers::PreconditionAMG());

  TrilinosWrappers::PreconditionAMG::AdditionalData Amg_data;
  Amg_data.constant_modes = constant_modes;
  Amg_data.elliptic = true;
  Amg_data.higher_order_elements = true;
  Amg_data.smoother_sweeps = 2;
  Amg_data.aggregation_threshold = 0.02;

  Mp_preconditioner->initialize (preconditioner_matrix.block(1,1));
  Amg_preconditioner->initialize (preconditioner_matrix.block(0,0),
                                  Amg_data);


  // SYSTEM MATRIX:
  auto A  = linear_operator< TrilinosWrappers::MPI::Vector >( matrix.block(0,0) );
  auto Bt = linear_operator< TrilinosWrappers::MPI::Vector >( matrix.block(0,1) );
  //  auto B =  transpose_operator(Bt);
  auto B     = linear_operator< TrilinosWrappers::MPI::Vector >( matrix.block(1,0) );
  auto ZeroP = 0*linear_operator< TrilinosWrappers::MPI::Vector >( matrix.block(1,1) );

  auto Mp    = linear_operator< TrilinosWrappers::MPI::Vector >( preconditioner_matrix.block(1,1) );

  static ReductionControl solver_control_pre(5000, 1e-8);
  static SolverCG<TrilinosWrappers::MPI::Vector> solver_CG(solver_control_pre);
  auto A_inv     = inverse_operator( A, solver_CG, *Amg_preconditioner);
  auto Schur_inv = inverse_operator( Mp, solver_CG, *Mp_preconditioner);

  auto P00 = A_inv;
  auto P01 = null_operator(Bt);
  auto P10 = Schur_inv * B * A_inv;
  auto P11 = -1 * Schur_inv;

  // ASSEMBLE THE PROBLEM:
  system_op  = block_operator<2, 2, VEC >({{
      {{ A, Bt }} ,
      {{ B, ZeroP }}
    }
  });


  //const auto S = linear_operator<VEC>(matrix);

  prec_op = block_operator<2, 2, VEC >({{
      {{ P00, P01 }} ,
      {{ P10, P11 }}
    }
  });
}


template class NeoHookeanTwoFieldsInterface <3,3>;

#endif
