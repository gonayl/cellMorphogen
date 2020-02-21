#include "SimuInfoModifier.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include <stdlib.h>
#include <math.h>

#include "SimulationParameters.hpp"



template<unsigned DIM>
double SimuInfoModifier<DIM>::nbrPart = 100;
template<unsigned DIM>
double SimuInfoModifier<DIM>::nbrPartInscrite = 1;



template<unsigned DIM>
SimuInfoModifier<DIM>::SimuInfoModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
SimuInfoModifier<DIM>::~SimuInfoModifier()
{
}

template<unsigned DIM>
void SimuInfoModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void SimuInfoModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void SimuInfoModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();




    double tempsParcouru = SimulationTime::Instance()->GetTime()/SimulationParameters::TIME_OF_SIMULATION;
    int pourcent = tempsParcouru * 100;

    double cible = nbrPartInscrite/nbrPart;


    if(tempsParcouru > cible){
      std::cout << pourcent << "% des pas de temps effectués" << std::endl;
      nbrPartInscrite = nbrPartInscrite + 1;
    }



}



template<unsigned DIM>
void SimuInfoModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class SimuInfoModifier<1>;
template class SimuInfoModifier<2>;
template class SimuInfoModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SimuInfoModifier)
