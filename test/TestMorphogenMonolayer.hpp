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
#include "MotileCellForce.hpp"
#include "CellLabel.hpp"
#include "CellEndo.hpp"
#include "CellLabelWriter.hpp"

#include "VertexMeshWriter.hpp"
#include "MorphogenTrackingModifier.hpp"
#include "PerimeterTrackingModifier.hpp"
#include "PerimeterDependentCellCycleModel.hpp"

#include <stdlib.h>

static const double M_TIME_FOR_SIMULATION = 40;
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
            //UniformCellCycleModel* p_cycle_model = new UniformCellCycleModel();
            FixedG1GenerationalCellCycleModel* p_cycle_model = new FixedG1GenerationalCellCycleModel();
            // UniformG1GenerationalCellCycleModel* p_cycle_model = new UniformG1GenerationalCellCycleModel();
            //MorphogenDependentCellCycleModel* p_cycle_model = new MorphogenDependentCellCycleModel();
            //p_cycle_model->SetDimension(2);
            //p_cycle_model->SetCurrentMass(0.5*(p_gen->ranf()+1.0));
            //p_cycle_model->SetMorphogenInfluence(1.0);

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

    void TestVertexBasedMorphogenMonolayerDirichletMotile() throw (Exception)
    {
        // Create Mesh
        /*
        HoneycombVertexMeshGenerator generator(3.0, 3.0);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        p_mesh->SetCellRearrangementThreshold(0.1);
        */

        VertexMeshReader<2,2> mesh_reader("testoutput/TestMorphogenMeshWriter/morphogen_mesh");
        MutableVertexMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);
        mesh.SetCellRearrangementThreshold(0.1);

       // p_mesh->Translate(-M_NUM_CELLS_ACROSS,-sqrt(3.0)*M_NUM_CELLS_ACROSS+ sqrt(3.0)/6.0);

        /* Remove all elements outside the specified initial radius
        for (VertexMesh<2,2>::VertexElementIterator elem_iter = p_mesh->GetElementIteratorBegin();
                 elem_iter != p_mesh->GetElementIteratorEnd();
                 ++elem_iter)
        {
            unsigned elem_index = elem_iter->GetIndex();
            c_vector<double,2> element_centre = p_mesh->GetCentroidOfElement(elem_index);

            if (norm_2(element_centre)>0.5*M_NUM_CELLS_ACROSS + 1e-5)
            {
                p_mesh->DeleteElementPriorToReMesh(elem_index);
            }
        }
        p_mesh->ReMesh(); */


        // Create Cells
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
           PerimeterDependentCellCycleModel* p_elong_model = new PerimeterDependentCellCycleModel();
           // FixedG1GenerationalCellCycleModel* p_model = new FixedG1GenerationalCellCycleModel();
           UniformG1GenerationalCellCycleModel* p_uni_model = new UniformG1GenerationalCellCycleModel();
           unsigned elem_index = mesh.GetElement(i)->GetIndex();
           // std::cout << elem_index << endl ;

           if (elem_index == 38 or elem_index == 2 or elem_index == 105 or elem_index == 141 or elem_index == 90 or elem_index == 72 or elem_index == 26 or elem_index == 131){
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
           } else if (elem_index == 74 or elem_index == 69 or elem_index == 18 or elem_index == 68){
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
             // p_cell->SetBirthTime(-20);
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
        MAKE_PTR(PerimeterTrackingModifier<2>, p_stretch_modifier);

        // Create Simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("MorphogenMonolayer/Gradient");
        /* simulator.SetDt(1.0/5.0);
        simulator.SetSamplingTimestepMultiple(5);
        simulator.SetEndTime(M_TIME_FOR_SIMULATION); */

        simulator.SetOutputDivisionLocations(true);

        simulator.AddSimulationModifier(p_stretch_modifier);

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

        // std::cout << "Growing Monolayer" << endl ;
        // simulator.Solve();


        std::cout << "VeGF diffusion" << endl ;

        // Create a parabolic PDE object - see the header file for what the constructor arguments mean
        MAKE_PTR_ARGS(CellwiseSourceParabolicPde<2>, p_pde, (cell_population, M_DUDT_COEFFICIENT,M_DIFFUSION_CONSTANT,M_UPTAKE_RATE,M_DECAY_COEFFICIENT,M_RADIUS));

        // Create a constant boundary conditions object, taking the value zero
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (0.0));

        // Create a ParabolicGrowingDomainPdeModifier, which is a simulation modifier, using the PDE and BC objects, and use 'true' to specify that you want a Dirichlet (rather than Neumann) BC
        MAKE_PTR_ARGS(ParabolicGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, true)); // Change this last argument to 'false' for no-flux BCs (Dirichlet??)

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
            if (cell_population.GetLocationIndexUsingCell(*cell_iter) == 74 or cell_population.GetLocationIndexUsingCell(*cell_iter) == 69 or cell_population.GetLocationIndexUsingCell(*cell_iter) == 68 or cell_population.GetLocationIndexUsingCell(*cell_iter) == 18)
            {
                cell_iter->AddCellProperty(p_label);
                cell_iter->AddCellProperty(p_endo);
            }
            else if (cell_population.GetLocationIndexUsingCell(*cell_iter) == 38 or cell_population.GetLocationIndexUsingCell(*cell_iter) == 2 or cell_population.GetLocationIndexUsingCell(*cell_iter) == 105 or cell_population.GetLocationIndexUsingCell(*cell_iter) == 141 or cell_population.GetLocationIndexUsingCell(*cell_iter) == 90 or cell_population.GetLocationIndexUsingCell(*cell_iter) == 72 \
            or cell_population.GetLocationIndexUsingCell(*cell_iter) == 26 or cell_population.GetLocationIndexUsingCell(*cell_iter) == 131)
            {
                cell_iter->AddCellProperty(p_label);
            }

        }

        MAKE_PTR(MorphogenTrackingModifier<2>, morphogen_modifier);
        simulator.AddSimulationModifier(morphogen_modifier);
        std::cout << "Adding active force" << endl ;
        MAKE_PTR_ARGS(MorphogenDrivenCellForce<2>, p_motile_force, (32,0.55));
        simulator.AddForce(p_motile_force);

        simulator.SetEndTime(50.0);
        simulator.SetDt(1.0/100.0);
        simulator.SetSamplingTimestepMultiple(1);

        simulator.Solve();


    }
};

#endif /* TESTMORPHOGENMONOLAYER_HPP_ */
