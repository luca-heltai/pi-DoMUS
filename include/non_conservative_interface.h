#include "interface.h"

using namespace deal2lkit;

template<int dim, int spacedim, int n_components, class Implementation,  typename LAC=LATrilinos>
class NonConservativeInterface : public Interface<dim,spacedim,n_components,LAC>
{

  typedef FEValuesCache<dim,spacedim> Scratch;
  typedef Assembly::CopyData::piDoMUSPreconditioner<dim,spacedim> CopyPreconditioner;
  typedef Assembly::CopyData::piDoMUSSystem<dim,spacedim> CopySystem;
public:

  virtual ~NonConservativeInterface() {};

  NonConservativeInterface(const std::string &name="",
                           const std::string &default_fe="FE_Q(1)",
                           const std::string &default_component_names="u",
                           const std::string &default_coupling="",
                           const std::string &default_preconditioner_coupling="",
                           const std::string &default_differential_components="") :
    Interface<dim,spacedim,n_components>(name, default_fe, default_component_names,
                                         default_coupling, default_preconditioner_coupling,
                                         default_differential_components) {};

  virtual void declare_parameters(ParameterHandler &prm)
  {
    Interface<dim,spacedim,n_components>::declare_parameters(prm);
  }
  virtual void parse_parameters_call_back()
  {
    Interface<dim,spacedim,n_components>::parse_parameters_call_back();
  }

  virtual void get_system_residual (const typename DoFHandler<dim,spacedim>::active_cell_iterator &cell,
                                    Scratch &scratch,
                                    CopySystem &data,
                                    std::vector<double> &local_residual) const
  {
    static_cast<const Implementation *>(this)->system_residual(cell, scratch, data, local_residual);
  }

  virtual void get_system_residual (const typename DoFHandler<dim,spacedim>::active_cell_iterator &cell,
                                    Scratch &scratch,
                                    CopySystem &data,
                                    std::vector<Sdouble> &local_residual) const
  {
    static_cast<const Implementation *>(this)->system_residual(cell, scratch, data, local_residual);
  }

  virtual void get_preconditioner_residual (const typename DoFHandler<dim,spacedim>::active_cell_iterator &cell,
                                            Scratch &scratch,
                                            CopyPreconditioner &data,
                                            std::vector<Sdouble> &local_residual) const
  {
    static_cast<const Implementation *>(this)->preconditioner_residual(cell, scratch, data, local_residual);
  }
};
