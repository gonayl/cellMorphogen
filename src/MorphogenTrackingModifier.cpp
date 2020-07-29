#include "MorphogenTrackingModifier.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "NodeBasedCellPopulation.hpp"
#include <cmath>
#include "Debug.hpp"
#include <iostream>
#include "CellLabel.hpp"
#include "CellEndo.hpp"
#include "CellTip.hpp"
#include "CellEpi.hpp"


template<unsigned DIM>
MorphogenTrackingModifier<DIM>::MorphogenTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
MorphogenTrackingModifier<DIM>::~MorphogenTrackingModifier()
{
}

template<unsigned DIM>
void MorphogenTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void MorphogenTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
    MARK;
}

template<unsigned DIM>
void MorphogenTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
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
      if (cell_iter->template HasCellProperty<CellEndo>())
      {
           cell_iter->GetCellData()->SetItem("morphogen",0.0);
      }
      if (cell_iter->template HasCellProperty<CellTip>())
      {
        std::vector<double> moy_x_morphogen_grad;
        std::vector<double> moy_y_morphogen_grad;

        double morphogen_grad_x = 0.0 ;
        double morphogen_grad_y = 0.0 ;

        // Get the set of neighbouring location indices
        std::set<unsigned> neighbour_indices = rCellPopulation.GetNeighbouringLocationIndices(*cell_iter);

        // If this cell has any neighbours (as defined by mesh/population/interaction distance)
        if (!neighbour_indices.empty() )
        {

            for (std::set<unsigned>::iterator neighbour_iter = neighbour_indices.begin();
                 neighbour_iter != neighbour_indices.end();
                 ++neighbour_iter)
            {
              CellPtr p_neighbour_cell = rCellPopulation.GetCellUsingLocationIndex(*neighbour_iter);
              bool neighbour_is_epi = p_neighbour_cell->template HasCellProperty<CellEpi>();

              if ( neighbour_is_epi == 1)
              {
                morphogen_grad_x = p_neighbour_cell->GetCellData()->GetItem("morphogen_grad_x");
                morphogen_grad_y = p_neighbour_cell->GetCellData()->GetItem("morphogen_grad_y");
                moy_x_morphogen_grad.push_back(morphogen_grad_x) ;
                moy_y_morphogen_grad.push_back(morphogen_grad_y) ;

              }
            }

        }



          double x = std::accumulate(moy_x_morphogen_grad.begin(), moy_x_morphogen_grad.end(), 0.0)/moy_x_morphogen_grad.size();
          double y = std::accumulate(moy_y_morphogen_grad.begin(), moy_y_morphogen_grad.end(), 0.0)/moy_y_morphogen_grad.size();

          cell_iter->GetCellData()->SetItem("morphogen_grad_x_moy",x);
          cell_iter->GetCellData()->SetItem("morphogen_grad_y_moy",y);
      }

    }
}

template<unsigned DIM>
void MorphogenTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class MorphogenTrackingModifier<1>;
template class MorphogenTrackingModifier<2>;
template class MorphogenTrackingModifier<3>;
// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MorphogenTrackingModifier)
