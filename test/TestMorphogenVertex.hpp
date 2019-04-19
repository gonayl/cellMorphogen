#ifndef TESTMORPHOGENVERTEX_HPP_
#define TESTMORPHOGENVERTEX_HPP_

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

#include "MorphogenDrivenCellForce.hpp"
#include "CellLabel.hpp"
#include "CellEndo.hpp"
#include "CellEpi.hpp"
#include "CellStalk.hpp"
#include "CellLabelWriter.hpp"

#include "VertexMeshWriter.hpp"
#include "MorphogenTrackingModifier.hpp"
#include "TargetAreaModifier.hpp"
#include "CellFixingModifier.hpp"

#include "PerimeterTrackingModifier.hpp"
#include "PerimeterDependentCellCycleModel.hpp"

#include <stdlib.h>

static const double M_TIME_FOR_SIMULATION = 40;
static const double M_NUM_CELLS_ACROSS = 10;
static const double M_UPTAKE_RATE = 10.0 ;
static const double M_DIFFUSION_CONSTANT = 5e-1;
static const double M_DUDT_COEFFICIENT = 1.0;
static const double M_DECAY_COEFFICIENT = 9.0;
static const double M_RADIUS = 100.0;

class TestMorphogenVertex : public AbstractCellBasedWithTimingsTestSuite
{

public:

    void TestVertexBasedMorphogenMonolayerDirichletMotile()
    {
        // Create Mesh
        /*
        HoneycombVertexMeshGenerator generator(3.0, 3.0);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        p_mesh->SetCellRearrangementThreshold(0.1);
        */

        VertexMeshReader<2,2> mesh_reader("testoutput/mesh/vertex_based_mesh");
        MutableVertexMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.SetCellRearrangementThreshold(0.1);

        // Create Cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
           PerimeterDependentCellCycleModel* p_elong_model = new PerimeterDependentCellCycleModel();
           //StochasticLumenCellCycleModel* p_uni_model = new StochasticLumenCellCycleModel();
           UniformG1GenerationalCellCycleModel* p_uni_model = new UniformG1GenerationalCellCycleModel();
           unsigned elem_index = mesh.GetElement(i)->GetIndex();
           // std::cout << elem_index << endl ;

           if (elem_index == 3 or elem_index == 17 or elem_index == 64 or elem_index == 54 or elem_index == 30 or elem_index == 49 or elem_index == 0 or elem_index == 2 ){
             CellPtr p_cell(new Cell(p_state, p_elong_model));
             p_cell->SetCellProliferativeType(p_transit_type);
             p_cell->SetBirthTime(-20);
             // Initial Condition for Morphogen PDE
             p_cell->GetCellData()->SetItem("morphogen",0.0);
             p_cell->GetCellData()->SetItem("morphogen_grad_x",0.0);
             p_cell->GetCellData()->SetItem("morphogen_grad_y",0.0);

             // Set Target Area so dont need to use a growth model in vertex simulations
             p_cell->GetCellData()->SetItem("target area", 1.0);
             cells.push_back(p_cell);
           } else if (elem_index == 42 or elem_index == 18 or elem_index == 61 or elem_index == 1){
             CellPtr p_cell(new Cell(p_state, p_elong_model));
             p_cell->SetCellProliferativeType(p_diff_type);
             p_cell->SetBirthTime(-20);
             // Initial Condition for Morphogen PDE
             p_cell->GetCellData()->SetItem("morphogen",0.0);
             p_cell->GetCellData()->SetItem("morphogen_grad_x",0.0);
             p_cell->GetCellData()->SetItem("morphogen_grad_y",0.0);

             // Set Target Area so dont need to use a growth model in vertex simulations
             p_cell->GetCellData()->SetItem("target area", 1.0);
             cells.push_back(p_cell);
           } else {
             CellPtr p_cell(new Cell(p_state, p_uni_model));
             p_cell->SetCellProliferativeType(p_transit_type);
             double birth_time = rand() % 15 + 0 ;
             p_cell->SetBirthTime(-birth_time);
             // Initial Condition for Morphogen PDE
             p_cell->GetCellData()->SetItem("morphogen",0.0);
             p_cell->GetCellData()->SetItem("morphogen_grad_x",0.0);
             p_cell->GetCellData()->SetItem("morphogen_grad_y",0.0);

             // Set Target Area so dont need to use a growth model in vertex simulations
             p_cell->GetCellData()->SetItem("target area", 1.0);
             cells.push_back(p_cell);
           }

        }

        VertexBasedCellPopulation<2> cell_population(mesh, cells);

        //Make cell data writer so can pass in variable name
        boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer(new CellDataItemWriter<2,2>("morphogen"));
        cell_population.AddCellWriter(p_cell_data_item_writer);
        cell_population.AddCellWriter<CellLabelWriter>();

        MAKE_PTR(CellLabel, p_label);
        MAKE_PTR(CellEndo, p_endo);
        MAKE_PTR(CellEndo, p_epi);
        MAKE_PTR(CellStalk, p_stalk);

        // Create Simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("CellMorphogen/TestMotileForce/Coupe7/Motile/TestFix/");
        /* simulator.SetDt(1.0/5.0);
        simulator.SetSamplingTimestepMultiple(5);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION); */

        simulator.SetOutputDivisionLocations(true);


        MAKE_PTR(PerimeterTrackingModifier<2>, p_stretch_modifier);
        simulator.AddSimulationModifier(p_stretch_modifier);

        MAKE_PTR(CellFixingModifier<2>, p_fixing_modifier);
        simulator.AddSimulationModifier(p_fixing_modifier);

        // MAKE_PTR(TargetAreaModifier<2>, p_growth_modifier);
        // simulator.AddSimulationModifier(p_growth_modifier);

        std::cout << "Adding passive force" << endl ;
        // Create Forces and pass to simulation
        MAKE_PTR(NagaiHondaDifferentialAdhesionForce<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(55.0);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(1.0);

        p_force->SetNagaiHondaCellCellAdhesionEnergyParameter(5.0);
        p_force->SetNagaiHondaLabelledCellLabelledCellAdhesionEnergyParameter(2.0);
        p_force->SetNagaiHondaLabelledCellCellAdhesionEnergyParameter(8.0);
        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(20.0);
        p_force->SetNagaiHondaLabelledCellBoundaryAdhesionEnergyParameter(8.0);
        simulator.AddForce(p_force);

        // std::cout << "Growing Monolayer" << endl ;
        // simulator.Solve();


        std::cout << "VeGF diffusion" << endl ;

        // Create a parabolic PDE object - see the header file for what the constructor arguments mean
        MAKE_PTR_ARGS(MorphogenCellwiseSourceParabolicPde<2>, p_pde, (cell_population, M_DUDT_COEFFICIENT,M_DIFFUSION_CONSTANT,M_RADIUS));

        // Create a constant boundary conditions object, taking the value zero
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (0.0));

        // Create a ParabolicGrowingDomainPdeModifier, which is a simulation modifier, using the PDE and BC objects, and use 'false' to specify that you want a Dirichlet (rather than Neumann) BC
        MAKE_PTR_ARGS(ParabolicGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false)); // Change this last argument to 'false' for fixed  BCs
        // Optionally, name the PDE state variable (for visualization purposes)
        p_pde_modifier->SetDependentVariableName("morphogen");

        // Add the modifier to the simulation and run it

        std::cout << "Adding Modifier" << endl ;

        p_pde_modifier->SetOutputGradient(true);
        p_pde_modifier->GetOutputGradient();

        simulator.AddSimulationModifier(p_pde_modifier);

        // boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer2(new CellDataItemWriter<2,2>("morphogen_grad_x"));
        // cell_population.AddCellWriter(p_cell_data_item_writer2);


        std::cout << "Labelling Motile Cell " << endl ;

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
            if (cell_population.GetLocationIndexUsingCell(*cell_iter) == 42 or cell_population.GetLocationIndexUsingCell(*cell_iter) == 1 or cell_population.GetLocationIndexUsingCell(*cell_iter) == 18 or cell_population.GetLocationIndexUsingCell(*cell_iter) == 61 )
            {
                cell_iter->AddCellProperty(p_label);
                cell_iter->AddCellProperty(p_endo);
            }
            else if (cell_population.GetLocationIndexUsingCell(*cell_iter) == 3 or cell_population.GetLocationIndexUsingCell(*cell_iter) == 0 or cell_population.GetLocationIndexUsingCell(*cell_iter) == 2 or cell_population.GetLocationIndexUsingCell(*cell_iter) == 17 or cell_population.GetLocationIndexUsingCell(*cell_iter) == 64
            or cell_population.GetLocationIndexUsingCell(*cell_iter) == 54 or cell_population.GetLocationIndexUsingCell(*cell_iter) == 30 or cell_population.GetLocationIndexUsingCell(*cell_iter) == 49)
            {
                cell_iter->AddCellProperty(p_label);
                cell_iter->AddCellProperty(p_stalk);
            }

        }

        MAKE_PTR(MorphogenTrackingModifier<2>, morphogen_modifier);
        simulator.AddSimulationModifier(morphogen_modifier);

        cell_population.AddCellWriter<CellAgesWriter>();

        std::cout << "Adding active force" << endl ;
        MAKE_PTR_ARGS(MorphogenDrivenCellForce<2>, p_motile_force, (15,0.55));
        simulator.AddForce(p_motile_force);

        simulator.SetEndTime(50.0);
        simulator.SetDt(1.0/50.0);
        simulator.SetSamplingTimestepMultiple(1);

        std::cout << "Growing monolayer" << endl ;

        simulator.Solve();


    }
};

#endif /* TESTMORPHOGENVERTEX_HPP_ */
