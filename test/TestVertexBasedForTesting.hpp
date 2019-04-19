#ifndef TESTVERTEXBASEDTESTING_HPP_
#define TESTVERTEXBASEDTESTING_HPP_

#include <cxxtest/TestSuite.h>

#include "CellBasedSimulationArchiver.hpp"

#include "SmartPointers.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"

#include "DifferentiatedCellProliferativeType.hpp"

#include "CellIdWriter.hpp"
#include "CellAgesWriter.hpp"
#include "VoronoiDataWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "NeighbourWriter.hpp"

#include "ParabolicGrowingDomainPdeModifier.hpp"
#include "MorphogenCellwiseSourceParabolicPde.hpp"
#include "VolumeTrackingModifier.hpp"

#include "PerimeterDependentCellCycleModel.hpp"
#include "MorphogenDependentCellCycleModel.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "UniformCellCycleModel.hpp"
#include "UniformCellCycleModel2.hpp"
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
#include "TargetAreaModifier.hpp"
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
#include "CellStalk.hpp"
#include "CellTip.hpp"
#include "CellPolar.hpp"
#include "CellBoundary.hpp"
#include "CellLabelWriter.hpp"
#include "CellTypeWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellAncestorWriter.hpp"


#include "VertexMeshWriter.hpp"
#include "MorphogenTrackingModifier.hpp"
#include "PerimeterTrackingModifier.hpp"
#include "LabelTrackingModifier.hpp"
#include "NeighbourTrackingModifier.hpp"
#include "PolarTrackingModifier.hpp"
#include "MassCenterTrackingModifier.hpp"
#include "PositionWeightTrackingModifier.hpp"
#include "PositionWeightConstantTrackingModifier.hpp"
#include "BorderTrackingModifier.hpp"
#include "CellFixingModifier.hpp"
#include "ForceTrackingModifier.hpp"
// #include "MergeNodeModifier.hpp"

#include "PlaneBoundaryCondition.hpp"
#include "FixedBoundaryCondition.hpp"

#include "PerimeterDependentCellCycleModel.hpp"
#include "StochasticLumenCellCycleModel.hpp"
#include "SimplePositionBasedCellCycleModel.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
using namespace std ;

class TestVertexBasedTesting : public AbstractCellBasedWithTimingsTestSuite
{
private:

    void GenerateCells(unsigned num_cells, std::vector<CellPtr>& rCells)
    {
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(WildTypeCellMutationState, p_state);

        for (unsigned i=0; i<num_cells; i++)
        {
            //StochasticLumenCellCycleModel* p_cycle_model = new StochasticLumenCellCycleModel();
            //UniformCellCycleModel2* p_cycle_model = new UniformCellCycleModel2();
            UniformG1GenerationalCellCycleModel* p_cycle_model = new UniformG1GenerationalCellCycleModel();
            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            // double birth_time = - RandomNumberGenerator::Instance()->ranf() * p_cycle_model->GetAverageTransitCellCycleTime();
            double birth_time = rand() % 10 + 5 ; // random double in the range 5 - 9
            p_cell->SetBirthTime(-birth_time);
            p_cell->InitialiseCellCycleModel();
            rCells.push_back(p_cell);
        }
     }




public:

    void TestVertexBasedThreeCellTypes()
    {
        // Create Mesh

        HoneycombVertexMeshGenerator generator(4, 4);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        p_mesh->SetCellRearrangementThreshold(0.1);

	      std::cout << "Creating mesh" << endl ;

        /*VertexMeshReader<2,2> mesh_reader("testoutput/mesh/vertex_based_mesh");
        MutableVertexMesh<2,2> p_mesh;
        p_mesh.ConstructFromMeshReader(mesh_reader);
        p_mesh.SetCellRearrangementThreshold(0.1);*/

        // Create Cells
        std::vector<CellPtr> cells;
        //MAKE_PTR(WildTypeCellMutationState, p_state);
        //MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        //MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        //CellsGenerator<StochasticLumenCellCycleModel, 2> cells_generator;
        //cells_generator.GenerateBasic(cells, p_mesh.GetNumElements(), std::vector<unsigned>(), p_transit_type);
        //cells_generator.GenerateBasic(cells, p_mesh->GetNumElements(), std::vector<unsigned>(), p_transit_type);

        //GenerateCells(p_mesh.GetNumElements(),cells);
        GenerateCells(p_mesh->GetNumElements(),cells);

        //VertexBasedCellPopulation<2> cell_population(p_mesh, cells);
        VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        // MAKE_PTR(MergeNodeModifier<2>, p_merge_modifier);
        // simulator.AddSimulationModifier(p_merge_modifier);

        //Make cell data writer so can pass in variable name
        cell_population.AddCellWriter<CellTypeWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellAncestorWriter>();

        MAKE_PTR(CellEndo, p_endo);
        MAKE_PTR(CellLumen, p_lumen);
	      MAKE_PTR(CellEpi, p_epi);
	      MAKE_PTR(CellLabel, p_label);
	      MAKE_PTR(CellStalk, p_stalk);

        // Create Simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("VertexModel/TestMeshReader/TestBoundaryConditions");


        simulator.SetOutputDivisionLocations(true);

        std::cout << "Adding passive force" << endl ;
        // Create Forces and pass to simulation
        MAKE_PTR(DifferentialAdhesionForce<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(55.0);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(1.0);

        p_force->SetEndoEndoAdhesionEnergyParameter(5.0);
        p_force->SetLumenLumenAdhesionEnergyParameter(5.0);
        p_force->SetEpiEpiAdhesionEnergyParameter(5.0);
        p_force->SetEndoEpiAdhesionEnergyParameter(5.0);
        p_force->SetEpiLumenAdhesionEnergyParameter(5.0);
        p_force->SetEndoLumenAdhesionEnergyParameter(5.0);

        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(10.0);
        p_force->SetEndoBoundaryAdhesionEnergyParameter(10.0);
        p_force->SetLumenBoundaryAdhesionEnergyParameter(10.0);
        p_force->SetEpiBoundaryAdhesionEnergyParameter(10.0);

        simulator.AddForce(p_force);

        MAKE_PTR(TargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        /*

        MAKE_PTR(MassCenterTrackingModifier<2>, p_center_modifier);
        simulator.AddSimulationModifier(p_center_modifier);

        MAKE_PTR(TestTrackingModifier<2>, p_test_modifier);
        simulator.AddSimulationModifier(p_test_modifier);

        MAKE_PTR(PositionWeightTrackingModifier<2>, p_distance_modifier);
        simulator.AddSimulationModifier(p_distance_modifier);

        MAKE_PTR(PositionWeightConstantTrackingModifier<2>, p_weight_modifier);
        simulator.AddSimulationModifier(p_weight_modifier);

        MAKE_PTR(NeighbourTrackingModifier<2>, p_neighbour_modifier);
        simulator.AddSimulationModifier(p_neighbour_modifier);

        MAKE_PTR(PolarTrackingModifier<2>, p_polar_modifier);
        simulator.AddSimulationModifier(p_polar_modifier);

        MAKE_PTR(LabelTrackingModifier<2>, p_label_modifier);
        simulator.AddSimulationModifier(p_label_modifier);
*/
        MAKE_PTR(ForceTrackingModifier<2>, p_force_modifier);
        simulator.AddSimulationModifier(p_force_modifier);

        MAKE_PTR(BorderTrackingModifier<2>, p_boundary_modifier);
        simulator.AddSimulationModifier(p_boundary_modifier);

        MAKE_PTR(PositionWeightTrackingModifier<2>, p_distance_modifier);
        simulator.AddSimulationModifier(p_distance_modifier);

        MAKE_PTR(VolumeTrackingModifier<2>, p_volume_modifier);
        simulator.AddSimulationModifier(p_volume_modifier);


        std::cout << "Importing label data from csv" << std::endl ;

        ifstream inFile ;
        int x ;


        inFile.open("testoutput/test_label.txt") ;
        std::vector<double> label_input;

        if(!inFile)
        {
          cout << "Unable to open file" << endl ;
        }
        while(inFile >> x)
        {
          label_input.push_back(x) ;

        }
        inFile.close() ;

        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        std::cout << "Labelling Epi cells" << endl ;

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {

          if (label_input[cell_population.GetLocationIndexUsingCell(*cell_iter)] == 1)
          {
              cell_iter->AddCellProperty(p_epi);
          }
          else if (label_input[cell_population.GetLocationIndexUsingCell(*cell_iter)] == 2)
          {
              cell_iter->AddCellProperty(p_endo);
              cell_iter->AddCellProperty(p_label);
              cell_iter->SetCellProliferativeType(p_diff_type);
           }
          else if (label_input[cell_population.GetLocationIndexUsingCell(*cell_iter)] == 3)
          {
              cell_iter->AddCellProperty(p_lumen);
          }

        }

        MAKE_PTR_ARGS(FixedBoundaryCondition<2>, p_bc, (&cell_population));
        simulator.AddCellPopulationBoundaryCondition(p_bc);


        std::cout << "Growing Monolayer" << endl ;

        simulator.SetEndTime(500.0);
        simulator.SetDt(1.0/100.0);
        simulator.SetSamplingTimestepMultiple(1.0);


        MAKE_PTR(CellFixingModifier<2>, p_fixing_modifier);
        simulator.AddSimulationModifier(p_fixing_modifier);


        simulator.Solve();




    }
};

#endif /* TESTVERTEXBASEDFORTESTING_HPP_ */
