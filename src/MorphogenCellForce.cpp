/*

Copyright (c) 2005-2015, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

/**********************************************
 * THIS CODE WORKS WITH RELEASE 3.3 OF CHASTE *
 **********************************************/

#include "MorphogenCellForce.hpp"
#include "CellLabel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include <numeric> 
#include <iostream> 
#include <cmath>

template<unsigned DIM>
MorphogenCellForce<DIM>::MorphogenCellForce(double a, double k)
    : AbstractForce<DIM>(),
      mStrength(a),
      mTreshold(k)
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

    // Helper variable that is a static cast of the cell population
    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);

    // Iterate over cells in the cell population
    
    c_vector<double,DIM> mass_x_centre;
    c_vector<double,DIM> mass_y_centre;

    double i = 0.0 ;

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)
     {
	c_vector<double, 2> centre_of_cell = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
	mass_x_centre(i) = centre_of_cell[0];
	mass_y_centre(i) = centre_of_cell[1];
	i = i + 1.0 ;

     }

    double xmoy = std::accumulate(mass_x_centre.begin(), mass_x_centre.end(), 0)/ mass_x_centre.size();
    double ymoy = std::accumulate(mass_y_centre.begin(), mass_y_centre.end(), 0)/ mass_x_centre.size();


    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
             cell_iter != rCellPopulation.End();
             ++cell_iter)

    	{
        if (cell_iter->template HasCellProperty<CellLabel>())
        {
            c_vector<double, 2> centre_of_cell2 = rCellPopulation.GetLocationOfCellCentre(*cell_iter);
            double x = centre_of_cell2[0];
            double y = centre_of_cell2[1];
	    c_vector<double,DIM> force;
            force(1) =   - mStrength * (y - ymoy) / (mTreshold + std::abs(y - ymoy)) ;
            force(0) =   - mStrength * (x - xmoy) / (mTreshold + std::abs(x - xmoy));

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
