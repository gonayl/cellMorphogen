#include "RepulsionForce.hpp"
#include "CellPeriph.hpp"
#include "VertexBasedCellPopulation.hpp"
#include <numeric>
#include <iostream>
#include <cmath>

template<unsigned DIM>
RepulsionForce<DIM>::RepulsionForce(double k)
    : AbstractForce<DIM>(),
      mTreshold(k)
{
assert(mTreshold > 0.0);
}

template<unsigned DIM>
RepulsionForce<DIM>::~RepulsionForce()
{
}

template<unsigned DIM>
void RepulsionForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Throw an exception message if not using a VertexBasedCellPopulation
    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == NULL)
    {
        EXCEPTION("RepulsionForce is to be used with a VertexBasedCellPopulation only");
    }

    // Helper variable that is a static cast of the cell population
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);

    unsigned num_nodes = p_cell_population->GetNumNodes();

    double treshold_reached;
    double neighbouring_cell;
    double same_cell;

    // Iterate over cells in the cell population
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)

    	{
        if (cell_iter->template HasCellProperty<CellPeriph>()) // only for peripherical cells (not core cells)
        {

	         VertexElement<DIM, DIM>* p_element = p_cell_population->GetElementCorrespondingToCell(*cell_iter); // Get the element corresponding to the cell
           double current_cell_location_index = p_cell_population->GetLocationIndexUsingCell(*cell_iter); // Get the location index of the CELL (see after)

           // Iterate over nodes owned by this VertexElement
		       unsigned num_nodes_in_vertex_element = p_element->GetNumNodes();
		       for (unsigned local_index=0; local_index<num_nodes_in_vertex_element; local_index++)
       		 {
             treshold_reached = 0 ;
             neighbouring_cell = 0 ;
             same_cell = 0 ;

             unsigned current_node_global_index = p_element->GetNodeGlobalIndex(local_index);
             Node<DIM>* p_current_node = p_cell_population->GetNode(current_node_global_index );
             c_vector<double, DIM> node_a_location = p_current_node->rGetLocation();
             std::set<unsigned> neighbours_of_current_cell = p_cell_population->GetNeighbouringLocationIndices(*cell_iter); // Given a cell, returns the set of location indices corresponding to neighbouring cells.


             if (p_current_node->IsBoundaryNode()) // only for boundary node within a peripherical cell
             {
               for (unsigned node_index=0; node_index<num_nodes; node_index++)
               {
                 Node<DIM>* p_node_b = p_cell_population->GetNode(node_index);
                 c_vector<double, DIM> node_b_location = p_node_b->rGetLocation();
                 c_vector<double, DIM> dist_nodes = p_cell_population->rGetMesh().GetVectorFromAtoB(node_a_location, node_b_location); // distance between current node (node A) and node B

                 std::set<unsigned> elements_containing_node_b = p_node_b->rGetContainingElementIndices(); // Return a set of indices of elements containing this node as a vertex.

                 // Iterate over cells containg node b
                 for (std::set<unsigned>::iterator containing_cells_iter = elements_containing_node_b .begin();
                 containing_cells_iter != elements_containing_node_b .end();
                 ++containing_cells_iter)
                 {
                   CellPtr p_containing_cell = p_cell_population->GetCellUsingLocationIndex(*containing_cells_iter);
                   double cell_b_location_index = p_cell_population->GetLocationIndexUsingCell(p_containing_cell);

                   if (neighbours_of_current_cell.find(cell_b_location_index) != neighbours_of_current_cell.end() && cell_b_location_index == current_cell_location_index)
                   {
                    neighbouring_cell++ ;// increment if one if the cell containing b is also a neighbour of the curent cell
                   }

                   if (cell_b_location_index == current_cell_location_index)
                   {
                    same_cell++ ;// increment if one if the cell containing b is also the curent cell
                   }

                 }

                 if (norm_2(dist_nodes) < mTreshold && neighbouring_cell == 0 && same_cell == 0 )
                 {
                  treshold_reached++ ;
                 }
               }

               if (treshold_reached > 0 )
               {
                 // Calculate the force applied on current node
                 c_vector<double, DIM> force = p_current_node->rGetAppliedForce();
                 c_vector<double, DIM> negative_force = -1 * force ;

                 // Add the force contribution to each node
                 p_current_node->AddAppliedForceContribution(negative_force);
               }

             }
          }
        }
    }
}

template<unsigned DIM>
void RepulsionForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class RepulsionForce<1>;
template class RepulsionForce<2>;
template class RepulsionForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RepulsionForce)
