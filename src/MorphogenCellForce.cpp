/**********************************************
 * THIS CODE WORKS WITH RELEASE 3.3 OF CHASTE *
 **********************************************/

#include "MorphogenCellForce.hpp"
#include "CellTip.hpp"
#include "VertexBasedCellPopulation.hpp"
#include <numeric>
#include <iostream>
#include <cmath>

template<unsigned DIM>
MorphogenCellForce<DIM>::MorphogenCellForce(double a)
    : AbstractForce<DIM>(),
      mStrength(a)
{
assert(mStrength > 0.0);
}

template<unsigned DIM>
MorphogenCellForce<DIM>::~MorphogenCellForce()
{
}

template<unsigned DIM>
void MorphogenCellForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Throw an exception message if not using a VertexBasedCellPopulation
    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == NULL)
    {
        EXCEPTION("MorphogenCellForce is to be used with a VertexBasedCellPopulation only");
    }
    double simulation_time = 48.0 ;
    // Helper variable that is a static cast of the cell population
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)

    	{
        if (cell_iter->template HasCellProperty<CellTip>())
        {
            double xmoy = cell_iter->GetCellData()->GetItem("mass_center_x");
            double ymoy = cell_iter->GetCellData()->GetItem("mass_center_y");
            c_vector<double, 2> centre_of_cell2 = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
            double x = centre_of_cell2[0];
            double y = centre_of_cell2[1];
            double norme = sqrt((y - ymoy) * (y - ymoy) + (x - xmoy) * (x - xmoy));
            c_vector<double,DIM> force;
            force(1) =   - mStrength * (1 + (SimulationTime::Instance()->GetTime()/simulation_time)) * (y - ymoy) / norme; // force that increse over time (simulate the increse in morphogen concentration)
            force(0) =   - mStrength * (1 + (SimulationTime::Instance()->GetTime()/simulation_time)) * (x - xmoy) / norme;

	    VertexElement<DIM, DIM>* p_element = p_cell_population->GetElementCorrespondingToCell(*cell_iter);

		// Iterate over nodes owned by this VertexElement
		unsigned num_nodes_in_vertex_element = p_element->GetNumNodes();
		for (unsigned local_index=0; local_index<num_nodes_in_vertex_element; local_index++)
       		 {
        		unsigned node_index = p_element->GetNodeGlobalIndex(local_index);
			p_cell_population->GetNode(node_index)->AddAppliedForceContribution(force);

       		 }
        }
    }
}

template<unsigned DIM>
void MorphogenCellForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class MorphogenCellForce<1>;
template class MorphogenCellForce<2>;
template class MorphogenCellForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MorphogenCellForce)
