#ifndef _NAVIER_STOKES_PROBLEM_
#define _NAVIER_STOKES_PROBLEM_

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>

#include <typeinfo>
#include <fstream>
#include <iostream>
#include <sstream>
#include <limits>
#include <locale>
#include <string>
#include <math.h>

#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/base/index_set.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>

#include "equation_data.h"
#include "linear_solver.h"
#include "assembly.h"
#include "solution.h"

#include "parsed_grid_generator.h"
#include "parsed_finite_element.h"
#include "error_handler.h"
#include "utilities.h"
#include "parsed_function.h"
#include "parsed_data_out.h"
#include "parameter_acceptor.h"

using namespace dealii;

template <int dim>
class NavierStokes : public ParameterAcceptor
{
	public:

		enum RefinementMode
		{
			global_refinement=0, 
			adaptive_refinement=1
		};

		NavierStokes (const RefinementMode refinement_mode);
    //    ~NavierStokes ();
		
		virtual void declare_parameters(ParameterHandler &prm);
		
		void run ();

	private:
		void make_grid_fe();
		void setup_dofs ();
		void assemble_navier_stokes_preconditioner ();
		void build_navier_stokes_preconditioner ();
		void assemble_navier_stokes_system ();
		void solve ();
		void output_results ();
		//void refine_mesh (const unsigned int max_grid_level);
		void refine_mesh ();
		void process_solution ();

		double       end_time;
		unsigned int initial_global_refinement;
		//unsigned int initial_global_refinement;
		unsigned int initial_adaptive_refinement;
		bool         generate_graphical_output;
		unsigned int graphical_output_interval;
		unsigned int adaptive_refinement_interval;
		double       stabilization_alpha;
		double       stabilization_c_R;
		double       stabilization_beta;
		unsigned int stokes_velocity_degree;
		bool         use_locally_conservative_discretization;

	private:
		ConditionalOStream                        pcout;

		shared_ptr<parallel::distributed::Triangulation<dim> > triangulation;

		double                                    global_Omega_diameter;

		const MappingQ<dim>                       mapping;

		shared_ptr<FiniteElement<dim,dim> >       navier_stokes_fe;

		shared_ptr<DoFHandler<dim> >           	  navier_stokes_dof_handler;

		ConstraintMatrix                          navier_stokes_constraints;

		TrilinosWrappers::BlockSparseMatrix       navier_stokes_matrix;
		TrilinosWrappers::BlockSparseMatrix       navier_stokes_preconditioner_matrix;

		TrilinosWrappers::MPI::BlockVector        navier_stokes_solution;
		TrilinosWrappers::MPI::BlockVector        old_navier_stokes_solution;
		TrilinosWrappers::MPI::BlockVector        navier_stokes_rhs;

		double                                    time_step;
		double                                    old_time_step;
		unsigned int                              timestep_number;

		std_cxx11::shared_ptr<TrilinosWrappers::PreconditionAMG>    Amg_preconditioner;
		std_cxx11::shared_ptr<TrilinosWrappers::PreconditionJacobi> Mp_preconditioner;
		std_cxx11::shared_ptr<TrilinosWrappers::PreconditionJacobi> T_preconditioner;

		bool                                      rebuild_navier_stokes_matrix;
		bool                                      rebuild_navier_stokes_preconditioner;

		TimerOutput                               computing_timer;

		void setup_navier_stokes_matrix (	const std::vector<IndexSet> &navier_stokes_partitioning,
											const std::vector<IndexSet> &navier_stokes_relevant_partitioning);
		void setup_navier_stokes_preconditioner (	const std::vector<IndexSet> &navier_stokes_partitioning,
													const std::vector<IndexSet> &navier_stokes_relevant_partitioning);


		void
		local_assemble_navier_stokes_preconditioner (const typename DoFHandler<dim>::active_cell_iterator &cell,
												Assembly::Scratch::NavierStokesPreconditioner<dim> &scratch,
												Assembly::CopyData::NavierStokesPreconditioner<dim> &data);

		void
		copy_local_to_global_navier_stokes_preconditioner (const Assembly::CopyData::NavierStokesPreconditioner<dim> &data);


		void
		local_assemble_navier_stokes_system (const typename DoFHandler<dim>::active_cell_iterator &cell,
										Assembly::Scratch::NavierStokesSystem<dim>  &scratch,
										Assembly::CopyData::NavierStokesSystem<dim> &data);

		void
		copy_local_to_global_navier_stokes_system (const Assembly::CopyData::NavierStokesSystem<dim> &data);

		class Postprocessor;

		const RefinementMode                    refinement_mode;

		ErrorHandler<1>                         eh;

		ParsedGridGenerator<dim,dim>            pgg;

		ParsedFiniteElement<dim,dim>            fe_builder;

		ParsedFunction<dim, dim+1>              boundary_conditions;

		ParsedFunction<dim, dim+1>              right_hand_side;

		ParsedDataOut<dim, dim+1> data_out;
};

template <int dim>
class NavierStokes<dim>::Postprocessor : public DataPostprocessor<dim>
{
public:
	Postprocessor (const unsigned int partition,
								const double       minimal_pressure);

	virtual
	void
	compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
																		std::vector<Vector<double> >                    &computed_quantities) const;

	virtual std::vector<std::string> get_names () const;

	virtual
	std::vector<DataComponentInterpretation::DataComponentInterpretation>
	get_data_component_interpretation () const;

	virtual UpdateFlags get_needed_update_flags () const;

private:
	const unsigned int partition;
	const double       minimal_pressure;
};

#endif