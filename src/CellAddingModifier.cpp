#include "CellAddingModifier.hpp"
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
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"

template<unsigned DIM>
CellAddingModifier<DIM>::CellAddingModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
CellAddingModifier<DIM>::~CellAddingModifier()
{
}

template<unsigned DIM>
void CellAddingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
  bool add_cell = false ;

   for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
       cell_iter != rCellPopulation.End();
       ++cell_iter)
  {
    // MARK; --> PASSED
    if (cell_iter->template HasCellProperty<CellLabel>())
      {
        double cell_age = cell_iter->GetAge();
        // MARK; --> PASSED
        if (cell_age == 5.0)
        {
            add_cell = true ;
             // MARK; --> PASSED
        }
      }
    }
    // Define the condition under which you want to add a cell (change as required)
    // bool some_condition_is_met = true;

    if (add_cell)
    {
        // Specify the location of the new cell (change as required)
        c_vector<double, DIM> new_cell_position = zero_vector<double>(DIM);
        new_cell_position(0) = 3.0;
        new_cell_position(1) = 3.0;

        // MARK; --> PASSED ;

        // Initialise a new node object
        Node<DIM>* p_new_node = new Node<DIM>(rCellPopulation.GetNumNodes(), new_cell_position, false);
        p_new_node->ClearAppliedForce();

        // (Here you would need to specify or copy node attributes, if you use any)

        // Add the node object to the cell population
        unsigned new_node_index = static_cast<NodesOnlyMesh<DIM>*>(&(rCellPopulation.rGetMesh()))->AddNode(p_new_node);
        // Initialise the new cell object (change as required)
        std::vector<CellPtr> cells;
        UniformG1GenerationalCellCycleModel* p_cycle_model = new UniformG1GenerationalCellCycleModel();
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_type);
        p_cycle_model->SetDimension(DIM);

        CellPtr p_new_cell(new Cell(p_state, p_cycle_model));
        p_new_cell->SetCellProliferativeType(p_type);
        p_new_cell->SetBirthTime(0.0);
        p_new_cell->InitialiseCellCycleModel();

        // Update cells vector (not quite sure you can do this step...)
        rCellPopulation.rGetCells().push_back(p_new_cell);

        // MARK; --> PASSED ;

        // Update mappings between cells and location indices (I hope this bit is correct...)
        rCellPopulation.SetCellUsingLocationIndex(new_node_index, p_new_cell);

        MARK; 

      }
}

template<unsigned DIM>
void CellAddingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{

}

template<unsigned DIM>
void CellAddingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class CellAddingModifier<1>;
template class CellAddingModifier<2>;
template class CellAddingModifier<3>;
// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CellAddingModifier)
