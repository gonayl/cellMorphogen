#include "AdhesionCoefModifier.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "NodeBasedCellPopulation.hpp"
#include <cmath>
#include "Debug.hpp"
#include <iostream>
#include "CellLabel.hpp"
#include "CellEndo.hpp"
#include "DifferentialAdhesionForce.hpp"


template<unsigned DIM>
AdhesionCoefModifier<DIM>::AdhesionCoefModifier()

    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
AdhesionCoefModifier<DIM>::~AdhesionCoefModifier()
{
}

template<unsigned DIM>
void AdhesionCoefModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    MAKE_PTR(DifferentialAdhesionForce<2>, p_force);
    p_force->SetCoreCoreAdhesionEnergyParameter(5.0 + 5.0*(SimulationTime::Instance()->GetTime()/10.0)) ;
    MARK;
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void AdhesionCoefModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */

    UpdateCellData(rCellPopulation);
    MARK;
}


template<unsigned DIM>
void AdhesionCoefModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();


}

template<unsigned DIM>
void AdhesionCoefModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class AdhesionCoefModifier<1>;
template class AdhesionCoefModifier<2>;
template class AdhesionCoefModifier<3>;
// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(AdhesionCoefModifier)
