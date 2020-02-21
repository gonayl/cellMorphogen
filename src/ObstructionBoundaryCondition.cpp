#include "ObstructionBoundaryCondition.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CellStalk.hpp"
#include "CellTip.hpp"
#include "CellVessel.hpp"
#include <stdlib.h>
using namespace std ;

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ObstructionBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::ObstructionBoundaryCondition(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>* pCellPopulation)
        : AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM>(pCellPopulation)
{

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ObstructionBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::ImposeBoundaryCondition(const std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> >& rOldLocations)
{
    if (dynamic_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(this->mpCellPopulation)==nullptr)
    {
        EXCEPTION("ObstructionBoundaryCondition requires a subclass of AbstractOffLatticeCellPopulation.");
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
            unsigned num_nodes = this->mpCellPopulation->GetNumNodes();
            for (unsigned node_index=0; node_index<num_nodes; node_index++)
            {
                Node<SPACE_DIM>* p_node = this->mpCellPopulation->GetNode(node_index);
                c_vector<double, SPACE_DIM> node_location = p_node->rGetLocation();

                double obsruction_radius = 2.0;

                c_vector<double, 2> obstruction_centre;
                obstruction_centre[0] = 2.0;
                obstruction_centre[1] = 4.8;

                double radius = norm_2(node_location-obstruction_centre);
                if (radius < obsruction_radius)
                {
                    p_node->rGetModifiableLocation() = obstruction_centre + obsruction_radius/radius*(node_location-obstruction_centre);
                }
            }

        }
    }
    else
    {
        // SPACE_DIM == 1
        NEVER_REACHED;
        //ObstructionBoundaryCondition::ImposeBoundaryCondition is not implemented in 1D
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool ObstructionBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::VerifyBoundaryCondition()
{
    return true;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ObstructionBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{

    // Call method on direct parent class
    AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

// Explicit instantiation
template class ObstructionBoundaryCondition<1,1>;
template class ObstructionBoundaryCondition<1,2>;
template class ObstructionBoundaryCondition<2,2>;
template class ObstructionBoundaryCondition<1,3>;
template class ObstructionBoundaryCondition<2,3>;
template class ObstructionBoundaryCondition<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(ObstructionBoundaryCondition)
