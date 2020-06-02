#include "CellAngioModifier.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include <cmath>
#include <stdlib.h>     /* srand, rand */
#include "Debug.hpp"
#include <iostream>
#include <complex>      // std::complex, std::polar
#include "CellEpi.hpp"
#include "CellEndo.hpp"
#include "CellStalk.hpp"
#include "CellPeriph.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"

template<unsigned DIM>
CellAngioModifier<DIM>::CellAngioModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
CellAngioModifier<DIM>::~CellAngioModifier()
{
}

template<unsigned DIM>
void CellAngioModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
   VertexBasedCellPopulation<DIM>* p_cell_population = dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) ;
   for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
       cell_iter != rCellPopulation.End();
       ++cell_iter)
  {
    double n_neighbour_endo = 0 ;
    // bool is_bnd = 0 ;

    // MARK; --> PASSED
    if (cell_iter->template HasCellProperty<CellEpi>() )
      {

        // Get the location index corresponding to this cell
        // unsigned index = p_cell_population->GetLocationIndexUsingCell(*cell_iter);

        // Get the set of neighbouring location indices
        std::set<unsigned> neighbour_indices = p_cell_population->GetNeighbouringLocationIndices(*cell_iter);

        // If this cell has any neighbours (as defined by mesh/population/interaction distance)...
        if (!neighbour_indices.empty() )
        {

            for (std::set<unsigned>::iterator neighbour_iter = neighbour_indices.begin();
                 neighbour_iter != neighbour_indices.end();
                 ++neighbour_iter)
            {
                // Determine whether this neighbour is labelled
                CellPtr p_neighbour_cell = p_cell_population->GetCellUsingLocationIndex(*neighbour_iter);
                bool neighbour_is_endo = p_neighbour_cell->template HasCellProperty<CellStalk>();

                if ( neighbour_is_endo == 1)
                {
                    // Here both cell_iter and p_neighbour_cell are labelled, so set type_of_link to 2
                    n_neighbour_endo++;
                }

            }
        }



      }

      if (n_neighbour_endo > 0 )
      {
        cell_iter->template RemoveCellProperty<CellEpi>();
        cell_iter->AddCellProperty(CellPropertyRegistry::Instance()->Get<CellEndo>());
        cell_iter->AddCellProperty(CellPropertyRegistry::Instance()->Get<CellStalk>());
      }

    }
}

template<unsigned DIM>
void CellAngioModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{

}

template<unsigned DIM>
void CellAngioModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class CellAngioModifier<1>;
template class CellAngioModifier<2>;
template class CellAngioModifier<3>;
// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CellAngioModifier)
