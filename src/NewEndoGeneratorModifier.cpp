#include "NewEndoGeneratorModifier.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CellEndo.hpp"
#include "CellStalk.hpp"
#include "CellTip.hpp"
#include "CellEpi.hpp"
#include "CellLumen.hpp"
#include "CellPolar.hpp"
#include "CellPeriph.hpp"
#include "CellCore.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "CounterSingleton.hpp"

#include <stdlib.h>
#include <numeric>
#include <iostream>
#include <cmath>
using namespace std ;


template<unsigned DIM>
NewEndoGeneratorModifier<DIM>::NewEndoGeneratorModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
NewEndoGeneratorModifier<DIM>::~NewEndoGeneratorModifier()
{
}

template<unsigned DIM>
void NewEndoGeneratorModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void NewEndoGeneratorModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void NewEndoGeneratorModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
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
        EXCEPTION("Border Tracking Modifier is to be used with a VertexBasedCellPopulation only");
    }

    // MAKE_PTR(CellBoundary, p_border);

    double time_max = 1.0 ;
    double chance_2_endo ;
    double count = CounterSingleton::Instance()->GetCount() ;
    double time_elapsed = SimulationTime::Instance()->GetTime() - count ;

    // Iterate over cell population
    //VertexBasedCellPopulation<DIM>* p_cell_population = dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) ;

    for (typename AbstractCellPopulation<DIM, DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {

      //VertexElement<DIM,DIM>* p_element = p_cell_population->GetElementCorrespondingToCell(*cell_iter);

      //RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
      double rando100 = rand() % 101 ;

      chance_2_endo = (rando100/100) + (time_elapsed / time_max)  ;
      //cout << chance_2_endo << endl ;

      if (cell_iter->template HasCellProperty<CellPeriph>() && cell_iter->template HasCellProperty<CellEpi>() && chance_2_endo > 1 )
      {

        cell_iter->template RemoveCellProperty<CellEpi>();
        cell_iter->AddCellProperty(CellPropertyRegistry::Instance()->Get<CellEndo>());
        cell_iter->AddCellProperty(CellPropertyRegistry::Instance()->Get<CellStalk>());
        cell_iter->SetCellProliferativeType(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());

        CounterSingleton::Instance()->IncrementCounter();

        //double countbis = CounterSingleton::Instance()->GetCount() ;

        //cout << countbis << " new endo cells" << endl;

      }


    }

}

template<unsigned DIM>
void NewEndoGeneratorModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class NewEndoGeneratorModifier<1>;
template class NewEndoGeneratorModifier<2>;
template class NewEndoGeneratorModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NewEndoGeneratorModifier)
