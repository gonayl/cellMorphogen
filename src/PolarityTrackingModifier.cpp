#include "PolarityTrackingModifier.hpp"
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
PolarityTrackingModifier<DIM>::PolarityTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
PolarityTrackingModifier<DIM>::~PolarityTrackingModifier()
{
}

template<unsigned DIM>
void PolarityTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void PolarityTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
    MARK;
}

template<unsigned DIM>
void PolarityTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
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



    // NodeBasedCellPopulation<DIM>* p_cell_population = static_cast<NodeBasedCellPopulation<DIM>*>(&(rCellPopulation));

    // Iterate over cell population
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
      if (cell_iter->template HasCellProperty<CellLabel>())
        {
        double magnitude = cell_iter->GetCellData()->GetItem("magnitude");
        double phase = cell_iter->GetCellData()->GetItem("phase") ;
        double phase2 = phase + 0.1;
        std::complex<double> polar (magnitude,phase2);
        double polabs = std::abs(polar) ;

        cell_iter->GetCellData()->SetItem("magnitude", magnitude);
        cell_iter->GetCellData()->SetItem("phase", phase2);
        cell_iter->GetCellData()->SetItem("polarity", polabs);
        }
    }
}

template<unsigned DIM>
void PolarityTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class PolarityTrackingModifier<1>;
template class PolarityTrackingModifier<2>;
template class PolarityTrackingModifier<3>;
// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PolarityTrackingModifier)
