#include "PositionTrackingModifier.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "NodeBasedCellPopulation.hpp"
#include <cmath>
#include <stdlib.h>     /* srand, rand */
#include "Debug.hpp"
#include <iostream>
#include <complex>      // std::complex, std::polar
#include "CellLabel.hpp"
#include "CellEndo.hpp"


template<unsigned DIM>
PositionTrackingModifier<DIM>::PositionTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
PositionTrackingModifier<DIM>::~PositionTrackingModifier()
{
}

template<unsigned DIM>
void PositionTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void PositionTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
    MARK;
}

template<unsigned DIM>
void PositionTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
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

  /* /  if (bool(dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation)))
    {
        static_cast<MeshBasedCellPopulation<DIM>*>(&(rCellPopulation))->CreateVoronoiTessellation();
        MARK;
    } */



    //NodeBasedCellPopulation<DIM>* p_cell_population = static_cast<NodeBasedCellPopulation<DIM>*>(&(rCellPopulation));

    // Iterate over cell population
  /*
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
      if (cell_iter->template HasCellProperty<CellEndo>())
        {

        c_vector<double, DIM> cell_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter) ;
        double xposini = cell_iter->GetCellData()->GetItem("xposini") ;
        double yposini = cell_iter->GetCellData()->GetItem("yposini") ;
        double xpos = cell_location[0] ;
        double ypos = cell_location[1] ;

        cell_iter->GetCellData()->SetItem("xpos", xpos);
        cell_iter->GetCellData()->SetItem("ypos", ypos);
        }
    }
    */

    /*for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
      if (cell_iter->template HasCellProperty<CellEndo>())
        {

        c_vector<double, DIM> cell_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter) ;
        double xposini = cell_iter->GetCellData()->GetItem("xposini") ;
        double yposini = cell_iter->GetCellData()->GetItem("yposini") ;
        double xpos = cell_location[0] ;
        double ypos = cell_location[1] ;
        double xdiff = cell_iter->GetCellData()->GetItem("xdiff") ;
        double ydiff = cell_iter->GetCellData()->GetItem("ydiff") ;

        double xdiff2 = std::abs(cell_location[0] - xpos) + xdiff ;
        double ydiff2 = std::abs(cell_location[1] - ypos) + ydiff ;
        double dist = sqrt((xpos - xposini)*(xpos - xposini) + (ypos - yposini)*(ypos - yposini)) ;

        cell_iter->GetCellData()->SetItem("dist", dist);
        cell_iter->GetCellData()->SetItem("xpos", xpos);
        cell_iter->GetCellData()->SetItem("ypos", ypos);
        }
    } */
}

template<unsigned DIM>
void PositionTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class PositionTrackingModifier<1>;
template class PositionTrackingModifier<2>;
template class PositionTrackingModifier<3>;
// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PositionTrackingModifier)
