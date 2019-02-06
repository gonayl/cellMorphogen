#include "PositionWeightTrackingModifier.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CellEndo.hpp"
#include "CellEpi.hpp"
#include "CellLumen.hpp"
#include "CellPolar.hpp"


template<unsigned DIM>
PositionWeightTrackingModifier<DIM>::PositionWeightTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
PositionWeightTrackingModifier<DIM>::~PositionWeightTrackingModifier()
{
}

template<unsigned DIM>
void PositionWeightTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void PositionWeightTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void PositionWeightTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
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
        EXCEPTION("Polar Tracking Modifier is to be used with a VertexBasedCellPopulation only");
    }

    // Iterate over cell population

    std::vector<double> max_dists;

    for (typename AbstractCellPopulation<DIM, DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {

      double xmoy = cell_iter->GetCellData()->GetItem("mass_center_x");
      double ymoy = cell_iter->GetCellData()->GetItem("mass_center_y");
      c_vector<double, DIM> cell_position = rCellPopulation.GetLocationOfCellCentre(*cell_iter) ;
      double xpos = cell_position[0] ;
      double ypos = cell_position[1] ;
      double dist = sqrt((xpos - xmoy)*(xpos - xmoy) + (ypos - ymoy)*(ypos - ymoy)) ;
      max_dists.push_back(dist);

      cell_iter->GetCellData()->SetItem("distancetocenter", dist);

    }

    double max_dist = *max_element(max_dists.begin(),max_dists.end());

    for (typename AbstractCellPopulation<DIM, DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {

      cell_iter->GetCellData()->SetItem("maxdistancetocenter", max_dist);

    }

}

template<unsigned DIM>
void PositionWeightTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class PositionWeightTrackingModifier<1>;
template class PositionWeightTrackingModifier<2>;
template class PositionWeightTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PositionWeightTrackingModifier)
