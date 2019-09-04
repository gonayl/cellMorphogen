#include "LabelTrackingModifier.hpp"
#include "CellEndo.hpp"
#include "CellLumen.hpp"
#include "CellEpi.hpp"
#include "CellLabel.hpp"
#include "CellPolar.hpp"
#include "CellPeriph.hpp"
#include "RandomNumberGenerator.hpp"
#include "AbstractSimpleGenerationalCellCycleModel.hpp"
#include <stdlib.h>

template<unsigned DIM>
LabelTrackingModifier<DIM>::LabelTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
LabelTrackingModifier<DIM>::~LabelTrackingModifier()
{
}

template<unsigned DIM>
void LabelTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void LabelTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void LabelTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
      // boost::shared_ptr<AbstractCellProperty> p_label = cell_iter->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<CellLabel>();
      // mpCell->RemoveCellProperty<CellLabel>();

      double proba_lumen = 0.2 ;
      double age = cell_iter->GetAge() ;
      double dt = SimulationTime::Instance()->GetTimeStep();
      //double n_endo = cell_iter->GetCellData()->GetItem("nendoneighbours");
      double n_polar = cell_iter->GetCellData()->GetItem("npolarneighbours");

      bool is_epi = cell_iter-> template HasCellProperty<CellEpi>() ;
      bool is_polar = cell_iter-> template HasCellProperty<CellPolar>() ;
      bool is_bnd = cell_iter-> template HasCellProperty<CellPeriph>() ;
      unsigned gen = static_cast<AbstractSimpleGenerationalCellCycleModel*>(cell_iter->GetCellCycleModel())->GetGeneration();


      // std::cout << "epi ? " << " " << is_epi << " " << "n endo? " << " " << n_endo << std::endl;

      if (is_epi ==1 && n_polar != 0 && is_polar == 0)
      {
        // std::cout << "Should divide into a lumen ! " << std::endl;
        proba_lumen = 0.99 ;
      }
      else
      {
        // std::cout << "Should stay an epithelial cell !" << std::endl;
        proba_lumen = 0.01 ;
      }

       if (age < dt )
      {
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        if (p_gen->ranf() < proba_lumen && gen > 0 && is_bnd == 0)
        {
          // std::cout << "Should divide into a lumen ! " << std::endl;
          cell_iter->AddCellProperty(rCellPopulation.GetCellPropertyRegistry()->template Get<CellLumen>());
          cell_iter->GetCellData()->SetItem("target area", 0.48);
        }
      }
    }
}

template<unsigned DIM>
void LabelTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class LabelTrackingModifier<1>;
template class LabelTrackingModifier<2>;
template class LabelTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(LabelTrackingModifier)
