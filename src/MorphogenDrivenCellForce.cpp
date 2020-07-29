#include "MorphogenDrivenCellForce.hpp"
#include "CellTip.hpp"
#include "CellEpi.hpp"
#include "VertexBasedCellPopulation.hpp"
#include <numeric>
#include <iostream>
#include <cmath>

template<unsigned DIM>
MorphogenDrivenCellForce<DIM>::MorphogenDrivenCellForce(double a, double k)
    : AbstractForce<DIM>(),
      mStrength(a),
      mTreshold(k)
{
assert(mStrength > 0.0);
}

template<unsigned DIM>
MorphogenDrivenCellForce<DIM>::~MorphogenDrivenCellForce()
{
}

template<unsigned DIM>
void MorphogenDrivenCellForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Throw an exception message if not using a VertexBasedCellPopulation
    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == NULL)
    {
        EXCEPTION("MorphogenDrivenCellForce is to be used with a VertexBasedCellPopulation only");
    }

    // Helper variable that is a static cast of the cell population
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);

    double xmoy = 0.0;
    double ymoy = 0.0;

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)

    	{
        if (cell_iter->template HasCellProperty<CellTip>())
        {

            double x =  cell_iter->GetCellData()->GetItem("morphogen_grad_x_moy");
            double y =  cell_iter->GetCellData()->GetItem("morphogen_grad_y_moy");

            

            //c_vector<double, 2> centre_of_cell2 = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
            //double x = cell_iter->GetCellData()->GetItem("morphogen_grad_x");
            //double y = cell_iter->GetCellData()->GetItem("morphogen_grad_y");

	          c_vector<double,DIM> force;

            force(0) =    mStrength * (x - xmoy) / (mTreshold + std::abs(x - xmoy)) ;
            force(1) =    mStrength * (y - ymoy) / (mTreshold + std::abs(y - ymoy)) ;

            /*
            double norme = sqrt((y - ymoy) * (y - ymoy) + (x - xmoy) * (x - xmoy));
            c_vector<double,DIM> force;
            force(1) =    mStrength  * (y - ymoy) / norme;
            force(0) =    mStrength  * (x - xmoy) / norme;
            */

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
void MorphogenDrivenCellForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class MorphogenDrivenCellForce<1>;
template class MorphogenDrivenCellForce<2>;
template class MorphogenDrivenCellForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MorphogenDrivenCellForce)
