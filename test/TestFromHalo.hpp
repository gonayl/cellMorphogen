#ifndef TESTFROMHALO_HPP_
#define TESTFROMHALO_HPP_

#include <cxxtest/TestSuite.h>

#include "CellBasedSimulationArchiver.hpp"

#include "SmartPointers.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "DifferentialAdhesionGeneralisedLinearSpringForce.hpp"
#include "TrianglesMeshWriter.hpp"
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
#include "UniformG1GenerationalBoundaryCellCycleModel.hpp"
#include "CellDataItemWriter.hpp"
#include "OffLatticeSimulation.hpp"
#include "OnLatticeSimulation.hpp"
#include "CellsGenerator.hpp"
#include "RandomCellKiller.hpp"
#include "DifferentialAdhesionGeneralisedLinearSpringForce.hpp"

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
#include "CellStalk.hpp"
#include "CellTip.hpp"
#include "CellEpi.hpp"
#include "CellPolar.hpp"
#include "CellVessel.hpp"
#include "CellBoundary.hpp"
#include "CellCore.hpp"
#include "CellPeriph.hpp"
#include "CellLabelWriter.hpp"
#include "CellTypeWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "CellPosWriter.hpp"
#include "CellLabelWriter.hpp"
#include "CellPosWriter.hpp"
#include "CellAdjacencyMatrixWriter.hpp"
#include "CalibrationErrorWriter.hpp"


#include "VertexMeshWriter.hpp"
#include "MorphogenTrackingModifier.hpp"
#include "PerimeterTrackingModifier.hpp"
#include "LabelTrackingModifier.hpp"
#include "NeighbourTrackingModifier.hpp"
#include "PolarTrackingModifier.hpp"
#include "MassCenterTrackingModifier.hpp"
//#include "PositionWeightTrackingModifier.hpp"
//#include "PositionWeightConstantTrackingModifier.hpp"
#include "BorderTrackingModifier.hpp"
#include "CellFixingModifier.hpp"
#include "AdhesionCoefModifier.hpp"
// #include "MergeNodeModifier.hpp"

#include "FixedBoundaryCondition.hpp"
#include "PerimeterDependentCellCycleModel.hpp"
#include "StochasticLumenCellCycleModel.hpp"
//#include "SimplePositionBasedCellCycleModel.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
using namespace std ;

static const double M_TIME_FOR_SIMULATION = 48;
static const double M_TIME_STEPS = 10;
static const double M_NUM_CELLS_ACROSS = 10;
static const double M_UPTAKE_RATE = 5.0 ;
static const double M_DIFFUSION_CONSTANT = 5e-1;
static const double M_DUDT_COEFFICIENT = 1.0;
static const double M_DECAY_COEFFICIENT = 9.0;
static const double M_RADIUS = 100.0;
static const double M_EPI = 5.0 ;
static const double M_EPIBND = 6.0 ;
static const double M_ENDOBND = 10.0 ;
static const double M_ENDOEPI = 10.0 ;
static const double M_MOTILITY = 15.0 ;


class TestFromHalo : public AbstractCellBasedWithTimingsTestSuite
{
private:

    void GenerateCells(unsigned num_cells, std::vector<CellPtr>& rCells, std::vector<double> label_input, std::vector<double> boundary_input)
    {
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        MAKE_PTR(WildTypeCellMutationState, p_state);

        // ICI : modifier durée du cycle cellulaire MAIS dans code source (cfr UniformG1UniformG1GenerationalBoundaryCellCycleModel.cpp)!

        for (unsigned i=0; i<num_cells; i++)
        {
            //StochasticLumenCellCycleModel* p_cycle_model = new StochasticLumenCellCycleModel();
            //UniformG1GenerationalCellCycleModel* p_cycle_model = new UniformG1GenerationalCellCycleModel();
            PerimeterDependentCellCycleModel* p_elong_model = new PerimeterDependentCellCycleModel();
            UniformG1GenerationalBoundaryCellCycleModel* p_cycle_model = new UniformG1GenerationalBoundaryCellCycleModel();

          if (label_input[i] == 0)
          {
            if (boundary_input[i] == 0)
            {
              p_cycle_model->SetCycleDuration(32) ;
            }
            else if (boundary_input[i] == 1 )
            {
              p_cycle_model->SetCycleDuration(11);
            }
            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            double birth_time = rand() % 5 + 1 ;
            p_cell->SetBirthTime(-birth_time);
            p_cell->InitialiseCellCycleModel();
            p_cell->GetCellData()->SetItem("target area", 0.8);
            // Initial Condition for Morphogen PDE
            p_cell->GetCellData()->SetItem("morphogen",0.0);
            p_cell->GetCellData()->SetItem("morphogen_grad_x",0.0);
            p_cell->GetCellData()->SetItem("morphogen_grad_y",0.0);
            rCells.push_back(p_cell);
          }
          else if (label_input[i] == 1 )
          {
            CellPtr p_cell(new Cell(p_state, p_elong_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            double birth_time = rand() % 10 + 1 ;
            p_cell->SetBirthTime(-birth_time);
            p_cell->InitialiseCellCycleModel();
            p_cell->GetCellData()->SetItem("target area", 0.8);
            // Initial Condition for Morphogen PDE
            p_cell->GetCellData()->SetItem("morphogen",0.0);
            p_cell->GetCellData()->SetItem("morphogen_grad_x",0.0);
            p_cell->GetCellData()->SetItem("morphogen_grad_y",0.0);
            rCells.push_back(p_cell);
          }
          else if (label_input[i] == 2 or label_input[i] == 3)
          {
            CellPtr p_cell(new Cell(p_state, p_elong_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            double birth_time = rand() % 10 + 1 ;
            p_cell->SetBirthTime(-birth_time);
            p_cell->InitialiseCellCycleModel();
            p_cell->GetCellData()->SetItem("target area", 0.8);
            // Initial Condition for Morphogen PDE
            p_cell->GetCellData()->SetItem("morphogen",0.0);
            p_cell->GetCellData()->SetItem("morphogen_grad_x",0.0);
            p_cell->GetCellData()->SetItem("morphogen_grad_y",0.0);
            rCells.push_back(p_cell);
          }
        }
     }

public:

    void TestVertexBasedMeshReader()
    {

        // Permet d'importer un fichier test et s'en servir pour taguer les cellules, voir ci-après
        std::cout << "Importing label data from txt" << std::endl ;
        ifstream inFile ;
        int x ;
        inFile.open("testoutput/test_label_vessel.txt") ;
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

        ifstream inFileBnd ;
        int x_bnd ;
        inFileBnd.open("testoutput/boundary_input_vessel.txt") ;
        std::vector<double> boundary_input;
        if(!inFileBnd)
        {
          cout << "Unable to open file" << endl ;
        }
        while(inFileBnd >> x_bnd)
        {
          boundary_input.push_back(x_bnd) ;
        }
        inFileBnd.close() ;
        cout << boundary_input.size() << endl ;

        // Create Mesh

        std::cout << "Creating mesh" << endl ;

        VertexMeshReader<2,2> mesh_reader("testoutput/mesh/vertex_based_mesh_vessel");
        MutableVertexMesh<2,2> p_mesh;
        p_mesh.ConstructFromMeshReader(mesh_reader);
        p_mesh.SetCellRearrangementThreshold(0.1);

        /*std::cout << "Creating mesh" << endl ;

        HoneycombVertexMeshGenerator generator(4, 4);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        p_mesh->SetCellRearrangementThreshold(0.1);*/

        // Create Cells
        std::vector<CellPtr> cells;
        GenerateCells(p_mesh.GetNumElements(),cells,label_input,boundary_input);

        VertexBasedCellPopulation<2> cell_population(p_mesh, cells);

        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellPosWriter>();
        cell_population.AddCellWriter<CellTypeWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();
        //cell_population.AddPopulationWriter<CellAdjacencyMatrixWriter>();
        cell_population.AddPopulationWriter<CalibrationErrorWriter>();

        MAKE_PTR(CellEpi, p_epi);
        MAKE_PTR(CellTip, p_tip);
        MAKE_PTR(CellEndo, p_endo);
        MAKE_PTR(CellLabel, p_label);
        MAKE_PTR(CellStalk, p_stalk);
        MAKE_PTR(CellVessel, p_vessel);

        // Create Simulation
        OffLatticeSimulation<2> simulator(cell_population);

        simulator.SetOutputDivisionLocations(true);

        // ICI : MODIFIER PARAMETRE D'ADHESION

        // double adhesion_coef = 200.0 ;

        std::cout << "Adding passive force" << endl ;
        // Create Forces and pass to simulation
        MAKE_PTR(DifferentialAdhesionForce<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(55.0);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(1.0);

        p_force->SetEndoEndoAdhesionEnergyParameter(4.0);
        p_force->SetLumenLumenAdhesionEnergyParameter(5.0);
        p_force->SetCoreCoreAdhesionEnergyParameter(M_EPI);
        p_force->SetCorePeriphAdhesionEnergyParameter(M_EPI);
        p_force->SetPeriphPeriphAdhesionEnergyParameter(M_EPI);
        p_force->SetEndoEpiAdhesionEnergyParameter(M_ENDOEPI);
        p_force->SetEpiLumenAdhesionEnergyParameter(5.0);
        p_force->SetEndoLumenAdhesionEnergyParameter(5.0);

        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(10.0);
        p_force->SetEndoBoundaryAdhesionEnergyParameter(M_ENDOBND);
        p_force->SetLumenBoundaryAdhesionEnergyParameter(10.0);
        p_force->SetEpiBoundaryAdhesionEnergyParameter(M_EPIBND);

        simulator.AddForce(p_force);


        MAKE_PTR(VolumeTrackingModifier<2>, p_volume_modifier);
        simulator.AddSimulationModifier(p_volume_modifier);
        MAKE_PTR(PerimeterTrackingModifier<2>, p_stretch_modifier);
        simulator.AddSimulationModifier(p_stretch_modifier);
        MAKE_PTR(BorderTrackingModifier<2>, p_border_modifier);
        simulator.AddSimulationModifier(p_border_modifier);
        MAKE_PTR(MassCenterTrackingModifier<2>, p_center_modifier);
        simulator.AddSimulationModifier(p_center_modifier) ;
        //MAKE_PTR(AdhesionCoefModifier<2>, p_coef_modifier);
        //simulator.AddSimulationModifier(p_coef_modifier);

        //MAKE_PTR(NeighbourTrackingModifier<2>, p_neighbour_modifier) ;
        //simulator.AddSimulationModifier(p_neighbour_modifier) ;
        //MAKE_PTR(PolarTrackingModifier<2>, p_polar_modifier) ;
        //simulator.AddSimulationModifier(p_polar_modifier) ;
        //MAKE_PTR(LabelTrackingModifier<2>, p_lumen_modifier) ;
        //simulator.AddSimulationModifier(p_lumen_modifier) ;


        // Diffusion de gradient, pas encore utile à ce stade (besoin pour simuler la motilité des cellules endo)

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


        // Labelling cells from file (Epi/Endo/Lumen)

        // On pourrait utiliser les informations de HALO pour taguer les cellules endo
        // (1 si positif au marqueur, 0 si non)

        std::cout << "Labelling Epi cells" << endl ;
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {

          if (label_input[cell_population.GetLocationIndexUsingCell(*cell_iter)] == 0)
          {
              cell_iter->AddCellProperty(p_epi);
          }
          else if (label_input[cell_population.GetLocationIndexUsingCell(*cell_iter)] == 1)
          {
              cell_iter->AddCellProperty(p_endo);
              cell_iter->AddCellProperty(p_stalk);
          }
          else if (label_input[cell_population.GetLocationIndexUsingCell(*cell_iter)] == 2)
          {
              cell_iter->AddCellProperty(p_endo);
              cell_iter->AddCellProperty(p_tip);
          }
          else if (label_input[cell_population.GetLocationIndexUsingCell(*cell_iter)] == 3)
          {
              cell_iter->AddCellProperty(p_endo);
              cell_iter->AddCellProperty(p_vessel);
          }
        }


        MAKE_PTR(MorphogenTrackingModifier<2>, morphogen_modifier);
        simulator.AddSimulationModifier(morphogen_modifier);


        // NE PAS DECOMMENTER LA SECTION SUIVANTE (bugs à régler)


        std::cout << "Adding active force" << endl ;
        MAKE_PTR_ARGS(MorphogenDrivenCellForce<2>, p_motile_force, (16,0.55));
        simulator.AddForce(p_motile_force);


        //MAKE_PTR(TargetAreaModifier<2>, p_growth_modifier);
        //simulator.AddSimulationModifier(p_growth_modifier);

        MAKE_PTR_ARGS(FixedBoundaryCondition<2>, p_fixed_bc, (&cell_population));
        simulator.AddCellPopulationBoundaryCondition(p_fixed_bc);


        std::cout << "Growing Monolayer" << endl ;

        simulator.SetEndTime(48.0);
        simulator.SetDt(1.0/10.0);
        simulator.SetSamplingTimestepMultiple(1.0);
        simulator.SetOutputDirectory("CellMorphogen/VertexModel/TestMeeting/6");

        simulator.Solve();

    }
};

#endif /* TESTFROMHALO_HPP_ */
