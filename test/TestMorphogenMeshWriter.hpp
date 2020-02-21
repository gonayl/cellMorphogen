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

#include "VolumeTrackingModifier.hpp"

#include "PerimeterDependentCellCycleModel.hpp"
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
#include "NagaiHondaDifferentialAdhesionForce.hpp"

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

#include "CellLabel.hpp"
#include "CellEndo.hpp"
#include "CellLabelWriter.hpp"
#include "CellLocationIndexWriter.hpp"

#include "VertexMeshWriter.hpp"

#include "PerimeterTrackingModifier.hpp"
#include "PerimeterDependentCellCycleModel.hpp"

static const double M_TIME_FOR_SIMULATION = 10;
static const double M_NUM_CELLS_ACROSS = 10;
static const double M_UPTAKE_RATE = 10.0 ;
static const double M_DIFFUSION_CONSTANT = 5e-1;
static const double M_DUDT_COEFFICIENT = 1.0;
static const double M_DECAY_COEFFICIENT = 9.0;
static const double M_RADIUS = 9.0;

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
            FixedG1GenerationalCellCycleModel* p_cycle_model = new FixedG1GenerationalCellCycleModel();

            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_transit_type);

            //p_cell->SetBirthTime(0.0);
            p_cell->SetBirthTime(-20);
            p_cell->InitialiseCellCycleModel();

            // Set Target Area so dont need to use a growth model in vertex simulations
            p_cell->GetCellData()->SetItem("target area", 1.0);
            rCells.push_back(p_cell);
        }
     }

public:

    void TestMorphogenMeshWriter()
    {
        // Create Mesh
        HoneycombVertexMeshGenerator generator(6.0, 6.0);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        p_mesh->SetCellRearrangementThreshold(0.1);

        // Create Cells
        std::vector<CellPtr> cells;
        GenerateCells(p_mesh->GetNumElements(),cells);

        // Create Population
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        MAKE_PTR(CellLabel, p_label);
        MAKE_PTR(CellEndo, p_endo);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        cell_population.AddCellWriter<CellLabelWriter>();
        cell_population.AddCellWriter<CellLocationIndexWriter>();

        // Create Simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("MorphogenMonolayer/MeshWriter");
        simulator.SetDt(1.0/5.0);
        simulator.SetSamplingTimestepMultiple(5);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION);

        simulator.SetOutputDivisionLocations(true);

        std::cout << "Adding passive force" << endl ;
        // Create Forces and pass to simulation
        MAKE_PTR(NagaiHondaDifferentialAdhesionForce<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(55.0);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(1.0);

        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(5.0);
        p_force->SetNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter(4.0);
        p_force->SetNagaiHondaLabelledCellCellAdhesionEnergyParameter(7.0);
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(15.0);
        p_force->SetNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter(10.0);
        simulator.AddForce(p_force);

        std::cout << "Growing Monolayer" << endl ;
        simulator.Solve();

        // Record mesh
        VertexMeshWriter<2,2> vertex_mesh_writer("TestMorphogenMeshWriter", "morphogen_mesh");
        vertex_mesh_writer.WriteFilesUsingMesh(*p_mesh);


    }
};

#endif /* TESTMORPHOGENMONOLAYER_HPP_ */
