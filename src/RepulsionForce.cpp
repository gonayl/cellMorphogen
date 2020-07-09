#include "RepulsionForce.hpp"
#include "CellEpi.hpp"
#include "CellPeriph.hpp"
#include "CellStalk.hpp"
#include "CellTip.hpp"
#include "CellVessel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include <numeric>
#include <iostream>
#include <cmath>

using namespace std ;

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

    /*
    //double treshold_reached;
    double neighbouring_cell;
    double same_cell;
    //double count = 0 ;


    // Iterate over cells in the cell population
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
    	{

        //double cell_location_index = p_cell_population->GetLocationIndexUsingCell(*cell_iter);
        if ( cell_iter->template HasCellProperty<CellPeriph>()) // only for peripherical cells (not core cells)
        {

	         VertexElement<DIM, DIM>* p_element = p_cell_population->GetElementCorrespondingToCell(*cell_iter); // Get the element corresponding to the cell
           double current_cell_location_index = p_cell_population->GetLocationIndexUsingCell(*cell_iter); // Get the location index of the CELL (see after)
          // cout << "Working on cell : " << current_cell_location_index << endl ;

           std::set<unsigned> neighbours_of_current_cell = p_cell_population->GetNeighbouringLocationIndices(*cell_iter); // Given a cell, returns the set of location indices corresponding to neighbouring cells.
           //cout << "that has " << neighbours_of_current_cell.size() << " neighbouring cells" << endl ;

           // Iterate over nodes owned by this VertexElement
		       unsigned num_nodes_in_vertex_element = p_element->GetNumNodes();
		       for (unsigned local_index=0; local_index<num_nodes_in_vertex_element; local_index++)
       		 {
             //treshold_reached = 0 ;
             neighbouring_cell = 0 ;
             same_cell = 0 ;

             unsigned current_node_global_index = p_element->GetNodeGlobalIndex(local_index);
             Node<DIM>* p_current_node = p_cell_population->GetNode(current_node_global_index );
             c_vector<double, DIM> node_a_location = p_current_node->rGetLocation();

             c_vector<double, DIM> repulsion_vector ; // initialise the "repulsion vector" (direction from node A to B)

             if (p_current_node->IsBoundaryNode()) // only for boundary node within a peripherical cell
             {
              // cout << "Working on bnd node : " << current_node_global_index << endl ;
               //std::cout << "node " << current_node_global_index << " is boundary" << std::endl ;

               for (unsigned node_index=0; node_index<num_nodes; node_index++) // loop over all nodes
               {
                Node<DIM>* p_node_b = p_cell_population->GetNode(node_index);
                if (p_node_b->IsBoundaryNode()) // only for boundary node within a peripherical cell
                {
                 c_vector<double, DIM> node_b_location = p_node_b->rGetLocation();
                 repulsion_vector = p_cell_population->rGetMesh().GetVectorFromAtoB(node_a_location, node_b_location); // distance between current node (node A) and node B

                if (norm_2(repulsion_vector) < mTreshold)
                {

                 std::set<unsigned> elements_containing_node_b = p_node_b->rGetContainingElementIndices(); // Return a set of indices of elements containing this node as a vertex.

                 // Iterate over cells containg node b
                 for (std::set<unsigned>::iterator containing_cells_iter = elements_containing_node_b .begin();
                 containing_cells_iter != elements_containing_node_b .end();
                 ++containing_cells_iter)
                 {
                   CellPtr p_containing_cell = p_cell_population->GetCellUsingLocationIndex(*containing_cells_iter);
                   double cell_b_location_index = p_cell_population->GetLocationIndexUsingCell(p_containing_cell);

                   if ((neighbours_of_current_cell.find(cell_b_location_index) != neighbours_of_current_cell.end()) == 1 )
                   {
                    neighbouring_cell++ ;// increment if one if the cell containing b is also a neighbour of the curent cell

                    //std::cout << " nodes " << current_node_global_index << " and " << node_index << " are too close but are on neighbouring cells" << std::endl ;
                   }

                   else if (cell_b_location_index == current_cell_location_index)
                   {
                    same_cell++ ;// increment if one if the cell containing b is also the curent cell

                    //std::cout << " nodes " << current_node_global_index << " and " << node_index << " are too close but are on the same cell" << std::endl ;
                   }

                 } // end iterating over cells containing node_b

                                  if (neighbouring_cell == 0 && same_cell == 0)
                                  {

                                 c_vector<double, DIM> force_current = p_current_node->rGetAppliedForce();

                                 //std::cout << "force applied before = " << norm_2(force_current) << std::endl ;

                                 c_vector<double, DIM> negative_force_current = - 1.01 * force_current ;

                                 //c_vector<double, DIM> repulsion_force_current =  repulsion_vector * 2.0 ;
                                //  c_vector<double, DIM> repulsion_force_b =  repulsion_vector * -2.0 ;

                                // std::cout << "repulsion force norm = " << norm_2(negative_force_current) << " applied on node " << current_node_global_index << std::endl ;

                                // std::cout << " nodes " << current_node_global_index << " from cell " << current_cell_location_index << " and " << node_index << " are too close " << std::endl ;

                                  p_current_node->ClearAppliedForce();
                                  p_node_b->ClearAppliedForce() ;

                                  //p_current_node->AddAppliedForceContribution(repulsion_force_current);
                                //  p_node_b->AddAppliedForceContribution(repulsion_force_b);

                                 c_vector<double, DIM> force_after = p_current_node->rGetAppliedForce();
                                 //std::cout << "force applied after = " << norm_2(force_after) << std::endl ;

                                 }

               } // end "if dist(current,b) < treshold"

               } // end "if node_b isbnb"

              } // end of iterating over all nodes

            } // end of "if current_node isbnb"
          }
        }
    }
    */

    // std::cout << "Starting repulsion force for bnd stalk cells" << std::endl ;

    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
        Node<DIM>* p_node = p_cell_population->GetNode(node_index);
        int stalk = 0 ;
        int tip = 0 ;
        int motile = 0 ;

        std::set<unsigned> elements_containing_node = p_cell_population->GetNode(node_index)->rGetContainingElementIndices();

        for (std::set<unsigned>::iterator element_index = elements_containing_node.begin();
             element_index != elements_containing_node.end();
             ++element_index)
        {
          CellPtr p_cell = p_cell_population->GetCellUsingLocationIndex(*element_index);
          if (p_cell->HasCellProperty<CellStalk>())
          {
            stalk++ ;
          }
          if (p_cell->HasCellProperty<CellTip>())
          {
            tip++ ;
          }
        }

        if (stalk > 0 && tip > 0 )
        {
          motile = 1 ;
        }

        /* else if (stalk == 0 && tip == 0 )
        {
          std::cout << "No tip or stalk cells" << std::endl ;
        } */

        for (std::set<unsigned>::iterator element_index = elements_containing_node.begin();
             element_index != elements_containing_node.end();
             ++element_index)
        {
          CellPtr p_cell = p_cell_population->GetCellUsingLocationIndex(*element_index);
          if (p_cell->HasCellProperty<CellStalk>() && p_node->IsBoundaryNode() && motile == 0)
          {
            /*
            double rand_mult = rand() % 16 + 95 ;
            rand_mult = rand_mult / 100 ;

            c_vector<double, DIM> force = p_node->rGetAppliedForce(); // get the force currently applied to the node to be fixed
            c_vector<double, DIM> new_force = - rand_mult * force ; // apply a force 90% lighter

            // Add the force contribution to each node
            p_node->AddAppliedForceContribution(new_force);
            */

            p_node->ClearAppliedForce();
          }
          else if (p_cell->HasCellProperty<CellVessel>())
          {
            /*
            double rand_mult = rand() % 16 + 95 ;
            rand_mult = rand_mult / 100 ;

            c_vector<double, DIM> force = p_node->rGetAppliedForce(); // get the force currently applied to the node to be fixed
            c_vector<double, DIM> new_force = - rand_mult * force ;

            // Add the force contribution to each node
            p_node->AddAppliedForceContribution(new_force); // add a opposite force contribution to dynamically fix the node
            */
            p_node->ClearAppliedForce() ;
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
