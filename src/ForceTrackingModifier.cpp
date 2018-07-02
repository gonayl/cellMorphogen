/*

Copyright (c) 2005-2017, University of Oxford.
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

#include "ForceTrackingModifier.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "NodeBasedCellPopulation.hpp"
#include <cmath>
#include "Debug.hpp"
#include <iostream>
#include "CellLabel.hpp"
#include "CellEndo.hpp"


template<unsigned DIM>
ForceTrackingModifier<DIM>::ForceTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
ForceTrackingModifier<DIM>::~ForceTrackingModifier()
{
}

template<unsigned DIM>
void ForceTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void ForceTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
    MARK;
}

template<unsigned DIM>
void ForceTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
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

  /* /  if (bool(dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation)))
    {
        static_cast<MeshBasedCellPopulation<DIM>*>(&(rCellPopulation))->CreateVoronoiTessellation();
        MARK;
    } */



    NodeBasedCellPopulation<DIM>* p_cell_population = static_cast<NodeBasedCellPopulation<DIM>*>(&(rCellPopulation));

    // Iterate over cell population
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        unsigned node_index = p_cell_population->GetLocationIndexUsingCell(*cell_iter);
        c_vector<double,DIM> force;
       //  c_vector<double,2> force_on_spring;
        // MAKE_PTR(GeneralisedLinearSpringForce<2>, linear_force);
        // DifferentialAdhesionGeneralisedLinearSpringForce<2> linear_force;

        force(0) = p_cell_population->rGetMesh().GetNode(node_index)->rGetAppliedForce()[0];
        force(1) = p_cell_population->rGetMesh().GetNode(node_index)->rGetAppliedForce()[1];

        double forcex = force(0);
        double forcey = force(1);
        double forcetot = sqrt((forcex * forcex) + (forcey * forcey));
        // Store the cell's volume in CellData
        cell_iter->GetCellData()->SetItem("forcex", forcex);
        cell_iter->GetCellData()->SetItem("forcey", forcey);
        cell_iter->GetCellData()->SetItem("forcetot", forcetot);


        /* std::set<unsigned> node_neighbours = p_cell_population->GetNeighbouringNodeIndices(node_index);

        for (std::set<unsigned>::iterator neighbour_iter = node_neighbours.begin();
        neighbour_iter != node_neighbours.end();
        ++neighbour_iter)
         {
          unsigned neighbour_index = *neighbour_iter;
          CellPtr p_neighbour_cell = p_cell_population->GetCellUsingLocationIndex(neighbour_index);
          if (p_neighbour_cell->HasCellProperty<CellEndo>())
          {
            Node<DIM>* p_node_a = rCellPopulation.GetNode(node_index);
            Node<DIM>* p_node_b = rCellPopulation.GetNode(neighbour_index);

            const c_vector<double, DIM>& node_a_location = p_node_a->rGetLocation();
            const c_vector<double, DIM>& node_b_location = p_node_b->rGetLocation();

            // Get the node radii for a NodeBasedCellPopulation
            double node_a_radius = 0.0;
            double node_b_radius = 0.0;

            if (bool(dynamic_cast<NodeBasedCellPopulation<DIM>*>(&rCellPopulation)))
            {
                node_a_radius = p_node_a->GetRadius();
                node_b_radius = p_node_b->GetRadius();
            }

            // Get the unit vector parallel to the line joining the two nodes
            c_vector<double, DIM> unit_difference;
            unit_difference = rCellPopulation.rGetMesh().GetVectorFromAtoB(node_a_location, node_b_location);

            // Calculate the distance between the two nodes
            double distance_between_nodes = norm_2(unit_difference);
            assert(distance_between_nodes > 0);
            assert(!std::isnan(distance_between_nodes));

            unit_difference /= distance_between_nodes;


            // Calculate the rest length of the spring connecting the two nodes with a default value of 1.0.
            double rest_length_final = 1.0;
            assert(node_a_radius > 0 && node_b_radius > 0);
            rest_length_final = node_a_radius+node_b_radius;
            double rest_length = rest_length_final;

            double overlap = distance_between_nodes - rest_length;
            bool is_closer_than_rest_length = (overlap <= 0);
            double multiplication_factor = 1.0 ;
            double spring_stiffness = 15.0;
            c_vector<double, DIM> force_on_spring;
            double llog = 1.0 ;
            double eexp = 1.0 ;

            // Get the force between these two cells, and its magnitude

            if (is_closer_than_rest_length) //overlap is negative
            {
                //log(x+1) is undefined for x<=-1
                assert(overlap > -rest_length_final);
                c_vector<double,DIM> force_on_spring = multiplication_factor*spring_stiffness * unit_difference * rest_length_final* log(1.0 + overlap/rest_length_final);
                // llog = log(1.0 + overlap/rest_length_final) ;
            }
            else
            {
                double alpha = 5.0;
                c_vector<double,DIM> force_on_spring = multiplication_factor*spring_stiffness * unit_difference * overlap * exp(-alpha * overlap/rest_length_final);
                // eexp = exp(-alpha * overlap/rest_length_final) ;
            }

            // double strength = norm_2(force_on_spring);

            // The next two lines only make sense if using a NodeBasedCellPopulation
            CellPtr p_cell_A =  rCellPopulation.GetCellUsingLocationIndex(node_index);
            CellPtr p_cell_B =  rCellPopulation.GetCellUsingLocationIndex(neighbour_index);

            // p_cell_A->GetCellData()->SetItem("force_on_spring", strength);
            //std::cout << strength ;
            //std::cout << "force x :" << forcex << "force y:" << forcey ;
            std::cout << force_on_spring.size() << "&" << force_on_spring[0] << "&" << force_on_spring[1];
          */ 
    }
}

template<unsigned DIM>
void ForceTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class ForceTrackingModifier<1>;
template class ForceTrackingModifier<2>;
template class ForceTrackingModifier<3>;
// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ForceTrackingModifier)
