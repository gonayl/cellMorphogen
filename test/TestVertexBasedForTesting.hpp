#ifndef TESTMORPHOGENMONOLAYER_HPP_
#define TESTMORPHOGENMONOLAYER_HPP_

#include <cxxtest/TestSuite.h>

#include "CellBasedSimulationArchiver.hpp"

#include "SmartPointers.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"

#include "DifferentiatedCellProliferativeType.hpp"

#include "CellIdWriter.hpp"
#include "CellAgesWriter.hpp"
#include "VoronoiDataWriter.hpp"
#include "CellMutationStatesWriter.hpp"

#include "ParabolicGrowingDomainPdeModifier.hpp"
#include "MorphogenCellwiseSourceParabolicPde.hpp"
#include "VolumeTrackingModifier.hpp"

#include "PerimeterDependentCellCycleModel.hpp"
#include "MorphogenDependentCellCycleModel.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "UniformCellCycleModel.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "CellDataItemWriter.hpp"
#include "OffLatticeSimulation.hpp"
#include "OnLatticeSimulation.hpp"
#include "CellsGenerator.hpp"
#include "RandomCellKiller.hpp"

#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"

#include "NodeBasedCellPopulation.hpp"
#include "RepulsionForce.hpp"

#include "VertexBasedCellPopulation.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "DifferentialAdhesionForce.hpp"

#include "PottsBasedCellPopulation.hpp"
#include "PottsMeshGenerator.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "SurfaceAreaConstraintPottsUpdateRule.hpp"

#include "CaBasedCellPopulation.hpp"
#include "DiffusionCaUpdateRule.hpp"
#include "ShovingCaBasedDivisionRule.hpp"
#include "AdhesionCaSwitchingUpdateRule.hpp"

#include "PetscSetupAndFinalize.hpp"

#include "MorphogenDrivenCellForce.hpp"
#include "CellLabel.hpp"
#include "CellEndo.hpp"
#include "CellLumen.hpp"
#include "CellEpi.hpp"
#include "CellLabelWriter.hpp"
#include "CellTypeWriter.hpp"

#include "VertexMeshWriter.hpp"
#include "MorphogenTrackingModifier.hpp"
#include "PerimeterTrackingModifier.hpp"
#include "PerimeterDependentCellCycleModel.hpp"
#include "StochasticLumenCellCycleModel.hpp"

#include <stdlib.h>

class TestMorphogenMonolayer : public AbstractCellBasedWithTimingsTestSuite
{
private:

    void GenerateCells(unsigned num_cells, std::vector<CellPtr>& rCells)
    {
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);

      //  RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

        for (unsigned i=0; i<num_cells; i++)
        {

            // StochasticLumenCellCycleModel* p_cycle_model = new StochasticLumenCellCycleModel();
            UniformG1GenerationalCellCycleModel* p_cycle_model = new UniformG1GenerationalCellCycleModel() ;

            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_transit_type);

            //p_cell->SetBirthTime(0.0);
            p_cell->SetBirthTime(-20);
            p_cell->InitialiseCellCycleModel();

            // Initial Condition for Morphogen PDE
            p_cell->GetCellData()->SetItem("morphogen",0.0);
            p_cell->GetCellData()->SetItem("morphogen_grad_x",0.0);
            p_cell->GetCellData()->SetItem("morphogen_grad_y",0.0);

            // Set Target Area so dont need to use a growth model in vertex simulations
            p_cell->GetCellData()->SetItem("target area", 1.0);
            rCells.push_back(p_cell);
        }
     }

public:

    void TestVertexBasedMorphogenMonolayerDirichletMotile()
    {
        // Create Mesh

        /* HoneycombVertexMeshGenerator generator(10, 10);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        p_mesh->SetCellRearrangementThreshold(0.1); */

	std::cout << "Creating mesh" << endl ;

	 VertexMeshReader<2,2> mesh_reader("testoutput/VertexModel/morphogen_mesh");
        MutableVertexMesh<2,2> p_mesh;
        p_mesh.ConstructFromMeshReader(mesh_reader);
        p_mesh.SetCellRearrangementThreshold(0.1);

        // Create Cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasic(cells, p_mesh.GetNumElements(), std::vector<unsigned>(), p_transit_type);

        VertexBasedCellPopulation<2> cell_population(p_mesh, cells);

        //Make cell data writer so can pass in variable name
        cell_population.AddCellWriter<CellTypeWriter>();

        MAKE_PTR(CellEndo, p_endo);
        MAKE_PTR(CellLumen, p_lumen);
	      MAKE_PTR(CellEpi, p_epi);

        // Create Simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("VertexModel/FromHalo/TestWriter");
        /* simulator.SetDt(1.0/5.0);
        simulator.SetSamplingTimestepMultiple(5);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION); */

        simulator.SetOutputDivisionLocations(true);

        std::cout << "Adding passive force" << endl ;
        // Create Forces and pass to simulation
        MAKE_PTR(DifferentialAdhesionForce<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(55.0);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(1.0);

        p_force->SetEndoEndoAdhesionEnergyParameter(5.0);
        p_force->SetLumenLumenAdhesionEnergyParameter(10.0);
        p_force->SetEndoEpiAdhesionEnergyParameter(5.0);
        p_force->SetLumenEpiAdhesionEnergyParameter(1.0);
        p_force->SetLumenEndoAdhesionEnergyParameter(1.0);

        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(10.0);
        p_force->SetEndoBoundaryAdhesionEnergyParameter(10.0);
        p_force->SetLumenBoundaryAdhesionEnergyParameter(10.0);

        simulator.AddForce(p_force);

        std::cout << "Labelling Random Cell " << endl ;

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)

             {
              if (RandomNumberGenerator::Instance()->ranf() < 0.33)
              {
                  cell_iter->AddCellProperty(p_endo);
              }
              else if (RandomNumberGenerator::Instance()->ranf() >= 0.33 && RandomNumberGenerator::Instance()->ranf() < 0.66)
              {
                  cell_iter->AddCellProperty(p_lumen);
              }
	            else
	            {
		              cell_iter->AddCellProperty(p_epi);
	            }
             }

        MAKE_PTR(SimpleTargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        std::cout << "Growing Monolayer" << endl ;

        simulator.SetEndTime(30.0);
        simulator.SetDt(1.0/100.0);
        simulator.SetSamplingTimestepMultiple(1);

        simulator.Solve();




    }
};

#endif /* TESTMORPHOGENMONOLAYER_HPP_ */
