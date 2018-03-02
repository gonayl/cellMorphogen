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

#ifndef TESTVERTEXACTIVEMIGRATION_HPP_
#define TESTVERTEXACTIVEMIGRATION_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "AbstractCellBasedWithTimingsTestSuite.hpp"

#include "SmartPointers.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "UniformCellCycleModel.hpp"
// #include "EndoCellCycle.hpp"
// #include "MorphogenDependentCellCycleModel.hpp"
#include "PerimeterDependentCellCycleModel.hpp"

#include "VertexBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "NagaiHondaForce.hpp"
#include "NagaiHondaDifferentialAdhesionForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "CellLabel.hpp"
#include "CellAgesWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "CellLabelWriter.hpp"
#include "CellIdWriter.hpp"
#include "NodeVelocityWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellVolumesWriter.hpp"
// #include "CellPerimetersWriter.hpp"

#include "CellPopulationElementWriter.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "StemCellProliferativeType.hpp"

#include "VoronoiDataWriter.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "HoneycombMeshGenerator.hpp"

#include "VertexMeshReader.hpp"
#include "VertexMeshWriter.hpp"

class TestCellCycle: public AbstractCellBasedTestSuite
{
private:

    double mLastStartTime;
    void setUp()
    {
        mLastStartTime = std::clock();
        AbstractCellBasedTestSuite::setUp();
    }
    void tearDown()
    {
        double time = std::clock();
        double elapsed_time = (time - mLastStartTime)/(CLOCKS_PER_SEC);
        std::cout << "Elapsed time: " << elapsed_time << std::endl;
        AbstractCellBasedTestSuite::tearDown();
    }


public:

    void TestCaseI() throw (Exception)
    {
        // Create a simple 2D MutableVertexMesh
        HoneycombVertexMeshGenerator generator(5, 5);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        p_mesh->SetCellRearrangementThreshold(0.1);

        // Set up cells, one for each VertexElement; specify cells as differentiated to prevent cell division
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        MAKE_PTR(WildTypeCellMutationState, p_state);


   for (unsigned i=0; i<p_mesh->GetNumElements(); i++)
   {
       PerimeterDependentCellCycleModel* p_model = new PerimeterDependentCellCycleModel();
       CellPtr p_cell(new Cell(p_state, p_model));
       p_cell->SetCellProliferativeType(p_transit_type);
       p_cell->SetBirthTime(-20);
       cells.push_back(p_cell);
   }

        VertexBasedCellPopulation<2> population(*p_mesh, cells);

        population.AddCellWriter<CellLabelWriter>();
        population.AddCellWriter<CellVolumesWriter>();

        OffLatticeSimulation<2> simulator(population);

        simulator.SetOutputDirectory("TestCellCycle/PerimeterDependent");
        simulator.SetEndTime(50.0);
     	  simulator.SetDt(0.01);
        simulator.SetSamplingTimestepMultiple(1);

        MAKE_PTR(NagaiHondaForce<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(55.0);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(0.0);
        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(5.0);
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(10.0);
        simulator.AddForce(p_force);

        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        MAKE_PTR(PerimeterTrackingModifier<2>, p_stretch_modifier);
        simulator.AddSimulationModifier(p_stretch_modifier);

	    std::cout << "Growing monolayer " << endl ;
        simulator.Solve();

        // Record mesh
        VertexMeshWriter<2,2> vertex_mesh_writer("TestCellCycle/PerimeterDependent", "vertex_mesh");
        MutableVertexMesh<2,2> vertex_mesh = static_cast<VertexBasedCellPopulation<2>(simulator.rGetCellPopulation()).rGetMesh();
        vertex_mesh_writer.WriteFilesUsingMesh(vertex_mesh);
    }

        void TestReadingInMesh() throw (Exception)
        {
            // Create mesh
            VertexMeshReader<2,2> mesh_reader("testoutput/TestCellCycle/PerimeterDependent/vertex_mesh");
            MutableVertexMesh<2,2> mesh;
            mesh.ConstructFromMeshReader(mesh_reader);
            mesh.SetCellRearrangementThreshold(0.1);

            // Set up cells, one for each VertexElement; specify cells as differentiated to prevent cell division
            std::vector<CellPtr> cells;
            MAKE_PTR(TransitCellProliferativeType, p_transit_type);
            MAKE_PTR(StemCellProliferativeType, p_stem_type);
            MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
            MAKE_PTR(WildTypeCellMutationState, p_state);


       for (unsigned i=0; i<mesh.GetNumElements(); i++)
       {
           PerimeterDependentCellCycleModel* p_model = new PerimeterDependentCellCycleModel();
           CellPtr p_cell(new Cell(p_state, p_model));
           p_cell->SetCellProliferativeType(p_transit_type);
           p_cell->SetBirthTime(-20);
           cells.push_back(p_cell);
       }

            VertexBasedCellPopulation<2> population(mesh, cells);

            population.AddCellWriter<CellLabelWriter>();
            population.AddCellWriter<CellVolumesWriter>();

            OffLatticeSimulation<2> simulator(population);

            simulator.SetOutputDirectory("TestCellCycle/ReadingInMesh");
            simulator.SetEndTime(50.0);
         	  simulator.SetDt(0.01);
            simulator.SetSamplingTimestepMultiple(1);

            MAKE_PTR(NagaiHondaForce<2>, p_force);
            p_force->SetNagaiHondaDeformationEnergyParameter(55.0);
            p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(0.0);
            p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(5.0);
            p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(10.0);
            simulator.AddForce(p_force);

            MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
            simulator.AddSimulationModifier(p_growth_modifier);

            MAKE_PTR(PerimeterTrackingModifier<2>, p_stretch_modifier);
            simulator.AddSimulationModifier(p_stretch_modifier);

    	    std::cout << "Growing monolayer " << endl ;
            simulator.Solve();
        }

    void TestMeshBasedMonolayer()
    {
      std::cout << "TestingCenterBasedModel" << endl ;
        HoneycombMeshGenerator generator(2, 2);
        MutableMesh<2,2>* p_mesh = generator.GetMesh();  //**Changed**//

        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        CellsGenerator<UniformG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumNodes(), p_transit_type);

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells); //**Changed**//

        cell_population.AddPopulationWriter<VoronoiDataWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestCenterBasedSimulation"); //**Changed**//
        simulator.SetSamplingTimestepMultiple(12);
        simulator.SetEndTime(20.0);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force); //**Changed**//
        simulator.AddForce(p_force);

        simulator.Solve();

        TS_ASSERT_EQUALS(cell_population.GetNumRealCells(), 16u);
        TS_ASSERT_DELTA(SimulationTime::Instance()->GetTime(), 20.0, 1e-10);
    }

};

#endif /*TESTVERTEXACTIVEMIGRATION_HPP_*/
