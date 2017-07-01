#include <pidomus.h>
#include "heat_interface.h"


using namespace dealii;
int main (int argc, char *argv[])
{

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,
                                                      numbers::invalid_unsigned_int);


  const int dim = 2;
  const int spacedim = 2;

  PoissonProblem<dim,spacedim,LADealII> p;
  piDoMUS<dim,spacedim,LADealII> solver ("pidomus",p);
  ParameterAcceptor::initialize("heat.prm", "used_parameters.prm");

  solver.run ();

  return 0;
}