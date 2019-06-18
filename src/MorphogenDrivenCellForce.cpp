#include "MorphogenDrivenCellForce.hpp"
#include "CellTip.hpp"
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

    // Iterate over cells in the cell population
    
    std::vector<double> moy_x_morphogen_grad;
    std::vector<double> moy_y_morphogen_grad;

    double morphogen_grad_x = 0.0 ;
    double morphogen_grad_y = 0.0 ;

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
     {

    morphogen_grad_x = cell_iter->GetCellData()->GetItem("morphogen_grad_x");
    morphogen_grad_y = cell_iter->GetCellData()->GetItem("morphogen_grad_y");
    moy_x_morphogen_grad.push_back(morphogen_grad_x) ;
    moy_y_morphogen_grad.push_back(morphogen_grad_y) ;
    
     }

    // double xmoy = std::accumulate(moy_x_morphogen_grad.begin(), moy_x_morphogen_grad.end(), 0)/moy_x_morphogen_grad.size();
    // double ymoy = std::accumulate(moy_y_morphogen_grad.begin(), moy_y_morphogen_grad.end(), 0)/moy_y_morphogen_grad.size();
    double xmoy = 0.0;
    double ymoy = 0.0;

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)

    	{
        if (cell_iter->template HasCellProperty<CellTip>())
        {
            c_vector<double, 2> centre_of_cell2 = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
            double x = cell_iter->GetCellData()->GetItem("morphogen_grad_x");
            double y = cell_iter->GetCellData()->GetItem("morphogen_grad_y");
	    c_vector<double,DIM> force;
            force(1) =    mStrength * (y - ymoy) / (mTreshold + std::abs(y - ymoy)) ;
            force(0) =    mStrength * (x - xmoy) / (mTreshold + std::abs(x - xmoy)) ;

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
