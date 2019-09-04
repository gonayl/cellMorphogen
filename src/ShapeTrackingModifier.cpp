#include "ShapeTrackingModifier.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"

template<unsigned DIM>
ShapeTrackingModifier<DIM>::ShapeTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
ShapeTrackingModifier<DIM>::~ShapeTrackingModifier()
{
}

template<unsigned DIM>
void ShapeTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void ShapeTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void ShapeTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();
    /**
     * This hack is needed because in the case of a MeshBasedCellPopulation in which
     * multiple cell divisions have occurred over one time step, the Voronoi tessellation
     * (while existing) is out-of-date. Thus, if we did not regenerate the Voronoi
     * tessellation here, an assertion may trip as we try to access a Voronoi element
     * whose index exceeds the number of elements in the out-of-date tessellation.
     *
     * \todo work out how to properly fix this (#1986)
     */
    if (bool(dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation)))
    {
        static_cast<MeshBasedCellPopulation<DIM>*>(&(rCellPopulation))->CreateVoronoiTessellation();
    }

    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == NULL)
    {
        EXCEPTION("Perimeter Tracking Modifier is to be used with a VertexBasedCellPopulation only");
    }

    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);

    // Iterate over cell population
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        // Get the perimeter of this cell

        VertexElement<DIM, DIM>* p_element = p_cell_population->GetElementCorrespondingToCell(*cell_iter);
        unsigned p_index = p_element->GetIndex();
        double cell_perimeter = p_cell_population->rGetMesh().GetElongationShapeFactorOfElement(p_index);
        // Store the cell's perimeter in CellData
        cell_iter->GetCellData()->SetItem("perimeter", cell_perimeter);
    }
}

template<unsigned DIM>
void ShapeTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class ShapeTrackingModifier<1>;
template class ShapeTrackingModifier<2>;
template class ShapeTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ShapeTrackingModifier)
