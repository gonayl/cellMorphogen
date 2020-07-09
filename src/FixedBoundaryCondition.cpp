#include "FixedBoundaryCondition.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CellStalk.hpp"
#include "CellTip.hpp"
#include "CellVessel.hpp"
#include <stdlib.h>
using namespace std ;

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
FixedBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::FixedBoundaryCondition(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>* pCellPopulation)
        : AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM>(pCellPopulation)
{

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void FixedBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::ImposeBoundaryCondition(const std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> >& rOldLocations)
{
    if (dynamic_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(this->mpCellPopulation)==nullptr)
    {
        EXCEPTION("FixedBoundaryCondition requires a subclass of AbstractOffLatticeCellPopulation.");
    }

    assert((dynamic_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(this->mpCellPopulation))
            || (SPACE_DIM==ELEMENT_DIM && (dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>*>(this->mpCellPopulation))) );

    if (SPACE_DIM != 1)
    {
        if (dynamic_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(this->mpCellPopulation))
        {
            // to be implemented
        }
        else
        {

            assert(SPACE_DIM == ELEMENT_DIM);
            assert(dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>*>(this->mpCellPopulation));
            // Iterate over all nodes and update their positions according to the boundary conditions
            //std::vector<double> pos_x;
            //std::vector<double> pos_y;

            unsigned num_nodes = this->mpCellPopulation->GetNumNodes();

            for (unsigned node_index=0; node_index<num_nodes; node_index++)
            {
                Node<SPACE_DIM>* p_node = this->mpCellPopulation->GetNode(node_index);
                c_vector<double, SPACE_DIM> node_location = p_node->rGetLocation();
                int stalk = 0 ;
                int tip = 0 ;
                int motile = 0 ;

                std::set<unsigned> elements_containing_node = this->mpCellPopulation->GetNode(node_index)->rGetContainingElementIndices();

                for (std::set<unsigned>::iterator element_index = elements_containing_node.begin();
                     element_index != elements_containing_node.end();
                     ++element_index)
                {
                  CellPtr p_cell = this->mpCellPopulation->GetCellUsingLocationIndex(*element_index);
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
                  CellPtr p_cell = this->mpCellPopulation->GetCellUsingLocationIndex(*element_index);
                  if (p_cell->HasCellProperty<CellStalk>() && p_node->IsBoundaryNode() && motile == 0)
                  {
                    c_vector<double, SPACE_DIM> old_node_location = rOldLocations.find(p_node)->second;

                    //cout << "node location : " << node_location[0] << " , " << node_location[1] << endl ;
                    //std::cout << node_index << " should be fixed at (" << old_node_location(0) << "," << old_node_location(1) << ")." << std::endl ;

                    // p_node->rGetModifiableLocation() = old_node_location; // fix the node at the current location

                    c_vector<double, SPACE_DIM> force = p_node->rGetAppliedForce(); // get the force currently applied to the node to be fixed
                    c_vector<double, SPACE_DIM> new_force = - 10 * force ; // apply a force 90% lighter

                    // Add the force contribution to each node
                    p_node->AddAppliedForceContribution(new_force);
                  }
                  else if (p_cell->HasCellProperty<CellVessel>())
                  {
                    c_vector<double, SPACE_DIM> old_node_location = rOldLocations.find(p_node)->second;

                    //cout << "node location : " << node_location[0] << " , " << node_location[1] << endl ;
                    //std::cout << node_index << " should be fixed at (" << old_node_location(0) << "," << old_node_location(1) << ")." << std::endl ;
                    // p_node->rGetModifiableLocation() = old_node_location;

                    c_vector<double, SPACE_DIM> force = p_node->rGetAppliedForce(); // get the force currently applied to the node to be fixed
                    c_vector<double, SPACE_DIM> new_force = - 10 * force ;

                    // Add the force contribution to each node
                    p_node->AddAppliedForceContribution(new_force); // add a opposite force contribution to dynamically fix the node
                  }


                }

            }

        }
    }
    else
    {
        // DIM == 1
        NEVER_REACHED;
        //FixedBoundaryCondition::ImposeBoundaryCondition is not implemented in 1D
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool FixedBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::VerifyBoundaryCondition()
{
    bool condition_satisfied = true;

    if (SPACE_DIM == 1)
    {
        EXCEPTION("FixedBoundaryCondition is not implemented in 1D");
    }
    else
    {

      unsigned num_nodes = this->mpCellPopulation->GetNumNodes();

      for (unsigned node_index=0; node_index<num_nodes; node_index++)
      {
          Node<SPACE_DIM>* p_node = this->mpCellPopulation->GetNode(node_index);
          c_vector<double, SPACE_DIM> node_location = p_node->rGetLocation();

          std::set<unsigned> elements_containing_node = this->mpCellPopulation->GetNode(node_index)->rGetContainingElementIndices();

          for (std::set<unsigned>::iterator element_index = elements_containing_node.begin();
               element_index != elements_containing_node.end();
               ++element_index)
          {
            CellPtr p_cell = this->mpCellPopulation->GetCellUsingLocationIndex(*element_index);

            if (p_cell->HasCellProperty<CellStalk>() && p_node->IsBoundaryNode() )
            {
              //Node<SPACE_DIM>* p_node = this->mpCellPopulation->GetNode(node_index);
              //c_vector<double, SPACE_DIM> node_location = p_node->rGetLocation();
              //c_vector<double, SPACE_DIM> old_node_location = rOldLocations.find(p_node)->second;
                /*
                if (node_location != old_node_location )
                {
                    // condition_satisfied = false ;
                    //break ; TO BE IMPLEMENTED

                    std::cout << "MARK;" << std::endl;
                }
                */
            }
          }
      }

  }

    return condition_satisfied;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void FixedBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{

    // Call method on direct parent class
    AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

// Explicit instantiation
template class FixedBoundaryCondition<1,1>;
template class FixedBoundaryCondition<1,2>;
template class FixedBoundaryCondition<2,2>;
template class FixedBoundaryCondition<1,3>;
template class FixedBoundaryCondition<2,3>;
template class FixedBoundaryCondition<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(FixedBoundaryCondition)
