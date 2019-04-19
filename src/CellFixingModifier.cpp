#include "CellFixingModifier.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include <cmath>
#include <stdlib.h>     /* srand, rand */
#include "Debug.hpp"
#include <iostream>
#include <complex>      // std::complex, std::polar
#include "CellLabel.hpp"
#include "CellEndo.hpp"
#include "CellStalk.hpp"
#include "CellTip.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"

template<unsigned DIM>
CellFixingModifier<DIM>::CellFixingModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
CellFixingModifier<DIM>::~CellFixingModifier()
{
}

template<unsigned DIM>
void CellFixingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
   VertexBasedCellPopulation<DIM>* p_cell_population = dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) ;
   for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
       cell_iter != rCellPopulation.End();
       ++cell_iter)
  {
    // MARK; --> PASSED
    if (cell_iter->template HasCellProperty<CellEndo>() ) // need to be changed to CellStalk
      {
        VertexElement<DIM, DIM>* p_element = p_cell_population->GetElementCorrespondingToCell(*cell_iter);

        // Iterate over nodes owned by this VertexElement
        unsigned num_nodes_in_vertex_element = p_element->GetNumNodes();
        for (unsigned local_index=0; local_index<num_nodes_in_vertex_element; local_index++)
             {

              unsigned node_index = p_element->GetNodeGlobalIndex(local_index);
              std::set<unsigned> elements_containing_node = p_cell_population->GetNode(node_index)->rGetContainingElementIndices();

              // std::cout << elements_containing_node.size() << std::endl ;
              double i = 0 ;

              for (std::set<unsigned>::iterator element_index = elements_containing_node.begin();
                   element_index != elements_containing_node.end();
                   ++element_index)
              {

              CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(*element_index);
              if (p_cell->HasCellProperty<CellTip>())
              {
                  // do nothing so the tip cell can freely move
              }
              else
              {
                  i++ ;
              }

              }

              //if (i > 0.0 && p_cell_population->GetNode(node_index)->IsBoundaryNode())
              if (i > 0.0 )
              {
                  p_cell_population->GetNode(node_index)->ClearAppliedForce();
                  std::cout << "The node " << node_index << " needs to be fixed" << std::endl ;
              }

             }
      }
    }
}

template<unsigned DIM>
void CellFixingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{

}

template<unsigned DIM>
void CellFixingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class CellFixingModifier<1>;
template class CellFixingModifier<2>;
template class CellFixingModifier<3>;
// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CellFixingModifier)
