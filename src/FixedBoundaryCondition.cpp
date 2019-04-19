#include "FixedBoundaryCondition.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CellLabel.hpp"

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
            /*
            if (SimulationTime::Instance()->GetTime() == 0)
            {
              std::vector<double> pos_x_ini ;
              std::vector<double> pos_y_ini ;
              for (unsigned node_index=0; node_index<num_nodes; node_index++)
              {
                Node<SPACE_DIM>* p_node = this->mpCellPopulation->GetNode(node_index);
                c_vector<double, SPACE_DIM> node_location = p_node->rGetLocation();
                pos_x_ini.push_back(node_location[0]) ;
                pos_y_ini.push_back(node_location[1]) ;
                }
            }
            */

            for (unsigned node_index=0; node_index<num_nodes; node_index++)
            {
                Node<SPACE_DIM>* p_node = this->mpCellPopulation->GetNode(node_index);
                c_vector<double, SPACE_DIM> node_location = p_node->rGetLocation();
                //pos_x.push_back(node_location[0]) ;
                //pos_y.push_back(node_location[1]) ;

                std::set<unsigned> elements_containing_node = this->mpCellPopulation->GetNode(node_index)->rGetContainingElementIndices();
                // double i = 0 ;

                for (std::set<unsigned>::iterator element_index = elements_containing_node.begin();
                     element_index != elements_containing_node.end();
                     ++element_index)
                {
                  CellPtr p_cell = this->mpCellPopulation->GetCellUsingLocationIndex(*element_index);
                  if (p_cell->HasCellProperty<CellLabel>())
                  {
                    c_vector<double, SPACE_DIM> new_pos;

                    new_pos(0) = node_location[0] ;
                    new_pos(1) = 0.0 ; // y-pos fixed
                    std::cout << node_index << " should be fixed at (" << new_pos(0) << "," << new_pos(1) << ")." << std::endl ;
                    p_node->rGetModifiableLocation() = new_pos;
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

            if (p_cell->HasCellProperty<CellLabel>() && p_node->IsBoundaryNode() && node_location[1] < 0.25)
            {
              if (node_location[1] != 0)
              {
                condition_satisfied = false ;
                break ;
              }
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
