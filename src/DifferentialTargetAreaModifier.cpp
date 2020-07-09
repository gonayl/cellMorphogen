#include "DifferentialTargetAreaModifier.hpp"
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
#include "CellEpi.hpp"
#include "CellCore.hpp"
#include "CellPeriph.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"

template<unsigned DIM>
DifferentialTargetAreaModifier<DIM>::DifferentialTargetAreaModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
DifferentialTargetAreaModifier<DIM>::~DifferentialTargetAreaModifier()
{
}

template<unsigned DIM>
void DifferentialTargetAreaModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{

  //std::cout << "Importing size data from txt" << std::endl ;
  ifstream inFile ;
  double x ;
  //inFile.open("testoutput/SimpleConditionInit/test_label_simple.txt") ;
  inFile.open("testoutput/cell_ini_areas.txt") ;
  std::vector<double> area_input;
  if(!inFile)
  {
    cout << "Unable to open file" << endl ;
  }
  while(inFile >> x)
  {
    area_input.push_back(x) ;
  }
  inFile.close() ;

   //VertexBasedCellPopulation<DIM>* p_cell_population = dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) ;
   for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
       cell_iter != rCellPopulation.End();
       ++cell_iter)
  {
    double cell_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter) ;

    // MARK; --> PASSED
    if (cell_iter->template HasCellProperty<CellEpi>() && cell_iter->template HasCellProperty<CellCore>() )
      {
        double target = 0.5 ;
        AbstractSimpleGenerationalCellCycleModel* p_model = static_cast<AbstractSimpleGenerationalCellCycleModel*>(cell_iter->GetCellCycleModel());
        double cell_area = area_input[cell_index];
        double cell_age = cell_iter->GetAge();
        double max_age = p_model->GetG1Duration() / 0.38 ;
      //  double target_area = cell_area * ((max_age - cell_age)/max_age) + (cell_age/max_age * 0.3) ;
        double target_area = sqrt(cell_area*(1 - (cell_age/max_age))) + (cell_age/max_age * target) ;
      //  std::cout << cell_index << " ; " << max_age << " ; " << cell_area << " ; " << cell_age << " ; " << target_area << std::endl;
        cell_iter->GetCellData()->SetItem("target area", target_area);
      }
    else if (cell_iter->template HasCellProperty<CellEpi>() && cell_iter->template HasCellProperty<CellPeriph>() )
      {
        double target = 0.5 ;
        AbstractSimpleGenerationalCellCycleModel* p_model = static_cast<AbstractSimpleGenerationalCellCycleModel*>(cell_iter->GetCellCycleModel());
        double cell_area = area_input[cell_index];
        double cell_age = cell_iter->GetAge();
        double max_age = p_model->GetG1Duration() / 0.38 ;
      //  double target_area = cell_area * ((max_age - cell_age)/max_age) + (cell_age/max_age * 0.3) ;
        double target_area = 0.8*sqrt(cell_area*(1 - (cell_age/max_age))) + (cell_age/max_age * target) ;
      //  std::cout << cell_index << " ; " << max_age << " ; " << cell_area << " ; " << cell_age << " ; " << target_area << std::endl;
        cell_iter->GetCellData()->SetItem("target area", target_area);
      }
      else if (cell_iter->template HasCellProperty<CellEndo>())
        {
          double target_area = 0.55 ;
          cell_iter->GetCellData()->SetItem("target area", target_area);
        }


    }
}

template<unsigned DIM>
void DifferentialTargetAreaModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{

}

template<unsigned DIM>
void DifferentialTargetAreaModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class DifferentialTargetAreaModifier<1>;
template class DifferentialTargetAreaModifier<2>;
template class DifferentialTargetAreaModifier<3>;
// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DifferentialTargetAreaModifier)
