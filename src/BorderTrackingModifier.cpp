#include "BorderTrackingModifier.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CellEndo.hpp"
#include "CellEpi.hpp"
#include "CellLumen.hpp"
#include "CellPolar.hpp"
#include "CellPeriph.hpp"
#include "CellCore.hpp"
#include <stdlib.h>


template<unsigned DIM>
BorderTrackingModifier<DIM>::BorderTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
BorderTrackingModifier<DIM>::~BorderTrackingModifier()
{
}

template<unsigned DIM>
void BorderTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void BorderTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void BorderTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
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

    double num_cells = 0 ;
    // Iterate over cell population
    VertexBasedCellPopulation<DIM>* p_cell_population = dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) ;
    for (typename AbstractCellPopulation<DIM, DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
      double n_boundary_nodes = 0 ;
      double n_periph_nodes = 0 ;
      double n_endo_neighbours = 0;
      VertexElement<DIM,DIM>* p_element = p_cell_population->GetElementCorrespondingToCell(*cell_iter);
      unsigned num_nodes_in_element = p_element->GetNumNodes();
      for (unsigned local_index=0; local_index<num_nodes_in_element; local_index++)
      {
          unsigned node_index = p_element->GetNodeGlobalIndex(local_index);

          if (rCellPopulation.GetNode(node_index)->IsBoundaryNode())
          {

              n_boundary_nodes++ ;

          }
      }

      std::set<unsigned> neighbour_indices = p_cell_population->GetNeighbouringLocationIndices(*cell_iter);
      // If this cell has any neighbours (as defined by mesh/population/interaction distance)...
      if (!neighbour_indices.empty() )
      {

          for (std::set<unsigned>::iterator neighbour_iter = neighbour_indices.begin();
               neighbour_iter != neighbour_indices.end();
               ++neighbour_iter)
          {

            CellPtr pnCell = p_cell_population->GetCellUsingLocationIndex(*neighbour_iter);

            bool neighbour_is_endo = pnCell->template HasCellProperty<CellEndo>();
            if(neighbour_is_endo)
            {
              n_endo_neighbours++ ;
            }
          }
      }

      if (n_boundary_nodes == 0 && n_endo_neighbours == 0)
      {
        cell_iter->template RemoveCellProperty<CellPeriph>();
        cell_iter->AddCellProperty(CellPropertyRegistry::Instance()->Get<CellCore>());
        for (unsigned local_index=0; local_index<num_nodes_in_element; local_index++)
        {
            unsigned node_index = p_element->GetNodeGlobalIndex(local_index);
            rCellPopulation.GetNode(node_index)->SetAsPeriphNode(false);

        }
      }
      else if (n_boundary_nodes > 0)
      {
        cell_iter->template RemoveCellProperty<CellCore>();
        cell_iter->AddCellProperty(CellPropertyRegistry::Instance()->Get<CellPeriph>());

        for (unsigned local_index=0; local_index<num_nodes_in_element; local_index++)
        {
            unsigned node_index = p_element->GetNodeGlobalIndex(local_index);
            rCellPopulation.GetNode(node_index)->SetAsPeriphNode(true);

        }

      }

      else if (n_endo_neighbours > 0)
      {
        cell_iter->template RemoveCellProperty<CellCore>();
        cell_iter->AddCellProperty(CellPropertyRegistry::Instance()->Get<CellPeriph>());
        for (unsigned local_index=0; local_index<num_nodes_in_element; local_index++)
        {
            unsigned node_index = p_element->GetNodeGlobalIndex(local_index);
            rCellPopulation.GetNode(node_index)->SetAsPeriphNode(true);

        }
      }

      else

      {
        EXCEPTION("Problem with boundary nodes") ;
      }

      cell_iter->GetCellData()->SetItem("nboundarynodes", n_boundary_nodes);
      cell_iter->GetCellData()->SetItem("nperiphnodes", n_periph_nodes);

      num_cells++;
    }


    double num_cells2 = 0 ;
    // Iterate over cell population

    for (typename AbstractCellPopulation<DIM, DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
      double n_periph_nodes = 0 ;
      VertexElement<DIM,DIM>* p_element = p_cell_population->GetElementCorrespondingToCell(*cell_iter);
      unsigned num_nodes_in_element = p_element->GetNumNodes();

      if (cell_iter->template HasCellProperty<CellPeriph>())
      {
        for (unsigned local_index=0; local_index<num_nodes_in_element; local_index++)
        {
            unsigned node_index = p_element->GetNodeGlobalIndex(local_index);
            rCellPopulation.GetNode(node_index)->SetAsPeriphNode(true);
            n_periph_nodes++ ;
        }
      }
      else if (cell_iter->template HasCellProperty<CellCore>())
      {
        for (unsigned local_index=0; local_index<num_nodes_in_element; local_index++)
        {
            unsigned node_index = p_element->GetNodeGlobalIndex(local_index);
            rCellPopulation.GetNode(node_index)->SetAsPeriphNode(false);
        }

      }

      else

      {
        EXCEPTION("Problem with boundary nodes") ;
      }

      cell_iter->GetCellData()->SetItem("nperiphnodes", n_periph_nodes);

      num_cells2++;
    }

}

template<unsigned DIM>
void BorderTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class BorderTrackingModifier<1>;
template class BorderTrackingModifier<2>;
template class BorderTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(BorderTrackingModifier)
