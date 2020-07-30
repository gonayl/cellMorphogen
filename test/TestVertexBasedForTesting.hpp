#ifndef TESTVERTEXBASEDFORTESTING_HPP_
#define TESTVERTEXBASEDFORTESTING_HPP_

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

#include "CounterSingletonEndo.hpp"

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

#include "MorphogenCellForce.hpp"
#include "MorphogenDrivenCellForce.hpp"
#include "CellLabel.hpp"
#include "CellEndo.hpp"
#include "CellLumen.hpp"
#include "CellStalk.hpp"
#include "CellTip.hpp"
#include "CellEpi.hpp"
#include "CellPolar.hpp"
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
#include "NewEndoGeneratorModifier.hpp"
#include "ObstructionWriterModifier.hpp"

//FOR LUMEN
#include "SimulationParameters.hpp"

#include "LumenGenerationModifier.hpp"
#include "LumenModifier.hpp"
#include "PolarisationModifier.hpp"
#include "SimuInfoModifier.hpp"

#include "CellPolarityXWriter.hpp"
#include "CellPolarityYWriter.hpp"

// #include "MergeNodeModifier.hpp"

#include "FixedBoundaryCondition.hpp"
#include "ObstructionBoundaryCondition.hpp"
#include "PerimeterDependentCellCycleModel.hpp"
#include "StochasticLumenCellCycleModel.hpp"
#include "SimplePositionBasedCellCycleModel.hpp"
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
static const double M_PEERIPHPERIPH = 5.0 ;
static const double M_EPIBND = 5.0 ;
static const double M_EPILUMEN = 10.0 ;
static const double M_ENDOBND = 4.0 ;//chang
static const double M_ENDOEPI = 6.0 ;
static const double M_ENDOENDO = 1.7 ;
static const double M_LUMENBND = 8.0 ;
static const double M_MOTILITY = 15.0 ;
static const double M_EPIEPI_INI = -0.008 ;
static const double M_ENDOEPI_INI = 0.024 ;
static const double M_LUMENEPI_INI = -0.015 ;
static const double M_POLARDEC_INI = 0.0075 ;
static const double M_LUMEN_SIZE_FACTOR_INI = 0.007 ;
static const double M_DURATION2_INI = 48.0 ;
static const double M_SIM_NB = 6 ;
static const double M_STRETCH_INI = 0.2 ;
static const double M_MAX_STRETCH = 2.5 ;
static const double M_STRETCH_PERIPH_INI = 0.2 ;

static const double M_PIXEL_COEF = 256 ;

class TestVertexBasedForTesting : public AbstractCellBasedWithTimingsTestSuite
{
private:

    void GenerateCells(unsigned num_cells, std::vector<CellPtr>& rCells)
    {
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        MAKE_PTR(WildTypeCellMutationState, p_state);


        // ICI : modifier durée du cycle cellulaire MAIS dans code source (cfr UniformG1UniformG1GenerationalBoundaryCellCycleModel.cpp)!

        for (unsigned i=0; i<num_cells; i++)
        {
            //StochasticLumenCellCycleModel* p_cycle_model = new StochasticLumenCellCycleModel();
            //UniformG1GenerationalCellCycleModel* p_cycle_model = new UniformG1GenerationalCellCycleModel();
            UniformG1GenerationalBoundaryCellCycleModel* p_cycle_model = new UniformG1GenerationalBoundaryCellCycleModel();
            PerimeterDependentCellCycleModel* p_elong_model = new PerimeterDependentCellCycleModel();
            if (i == 188 || i == 239 || i == 35 || i == 160)
            {
              CellPtr p_cell(new Cell(p_state, p_elong_model));
              p_cell->SetCellProliferativeType(p_diff_type);

              p_cycle_model->SetMaxTransitGenerations(10);
              double birth_time = 10 ;
              p_cell->SetBirthTime(-birth_time);
              p_cell->InitialiseCellCycleModel();
              p_cell->GetCellData()->SetItem("target area", 0.8);
              p_elong_model->SetMaxStretch(3.5);
              p_elong_model->SetMaxStretchPeriph(6.0);
              // Initial Condition for Morphogen PDE
              p_cell->GetCellData()->SetItem("morphogen",0.0);
              p_cell->GetCellData()->SetItem("morphogen_grad_x",0.0);
              p_cell->GetCellData()->SetItem("morphogen_grad_y",0.0);

            rCells.push_back(p_cell);
            }
            /*
            else if (i == 212 || i == 109 || i == 103 || i == 81)
            {
              CellPtr p_cell(new Cell(p_state, p_cycle_model));
              p_cell->SetCellProliferativeType(p_diff_type);
              p_cycle_model->SetMaxTransitGenerations(10);
              double birth_time = 10 ;
              p_cell->SetBirthTime(-birth_time);
              p_cell->InitialiseCellCycleModel();
              p_cell->GetCellData()->SetItem("target area", 0.5);
              // Initial Condition for Morphogen PDE
              p_cell->GetCellData()->SetItem("morphogen",0.0);
              p_cell->GetCellData()->SetItem("morphogen_grad_x",0.0);
              p_cell->GetCellData()->SetItem("morphogen_grad_y",0.0);

              rCells.push_back(p_cell);
            }
            */
            else
            {
              CellPtr p_cell(new Cell(p_state, p_cycle_model));
              p_cell->SetCellProliferativeType(p_transit_type);
              p_cycle_model->SetMaxTransitGenerations(10);
              double birth_time = rand() % 10 + 5 ;
              p_cell->SetBirthTime(-birth_time);
              p_cell->InitialiseCellCycleModel();
              p_cell->GetCellData()->SetItem("target area", 0.5);
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

        // Create Mesh

        std::cout << "Creating mesh" << endl ;
        /*
        HoneycombVertexMeshGenerator generator(4, 4);
        MutableVertexMesh<2,2>* p_mesh = generator.GetMesh();
        p_mesh->SetCellRearrangementThreshold(0.1);
        */

        VertexMeshReader<2,2> mesh_reader("testoutput/mesh/morphogen_mesh");
        MutableVertexMesh<2,2> p_mesh;
        p_mesh.ConstructFromMeshReader(mesh_reader);
        p_mesh.SetCellRearrangementThreshold(0.1);

	      std::cout << "Creating mesh" << endl ;

        // Create Cells
        std::vector<CellPtr> cells;

        GenerateCells(p_mesh.GetNumElements(),cells);
        //GenerateCells(p_mesh->GetNumElements(),cells);

        VertexBasedCellPopulation<2> cell_population(p_mesh, cells);
        //VertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellPosWriter>();
        cell_population.AddCellWriter<CellTypeWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();

        MAKE_PTR(CellEpi, p_epi);
        MAKE_PTR(CellEndo, p_endo);
        MAKE_PTR(CellLabel, p_label);
        MAKE_PTR(CellStalk, p_stalk);
        MAKE_PTR(CellTip, p_tip);
        MAKE_PTR(CellLumen, p_lumen);
        MAKE_PTR(CellPeriph, p_periph);
        // Create Simulation
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("CellMorphogen/VertexModel/TestRepulsionForce/MeshReader/WithNewObstru/WithNewIC/2"); // 1 with deformationparam = 55, 2 = 99


        simulator.SetOutputDivisionLocations(true);

        // ICI : MODIFIER PARAMETRE D'ADHESION
        std::cout << "Adding passive force" << endl ;
        // Create Forces and pass to simulation
        MAKE_PTR(DifferentialAdhesionForce<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(99.0);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(1.0);

        p_force->SetEndoEndoAdhesionEnergyParameter(M_ENDOENDO*1.0);
        p_force->SetStalkStalkAdhesionEnergyParameter(M_ENDOENDO*1.0);
        p_force->SetStalkTipAdhesionEnergyParameter(M_ENDOENDO*1.0);
        p_force->SetTipTipAdhesionEnergyParameter(M_ENDOENDO*1.0);
        p_force->SetLumenLumenAdhesionEnergyParameter(5.0);
        p_force->SetCoreCoreAdhesionEnergyParameter(M_EPI);
        p_force->SetCorePeriphAdhesionEnergyParameter(M_EPI);
        p_force->SetPeriphPeriphAdhesionEnergyParameter(M_PEERIPHPERIPH);
        p_force->SetEndoEpiAdhesionEnergyParameter(M_ENDOEPI);
        p_force->SetEpiLumenAdhesionEnergyParameter(6.0);
        p_force->SetEndoLumenAdhesionEnergyParameter(35.0);

        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(10.0);
        p_force->SetEndoBoundaryAdhesionEnergyParameter(5.0);
        p_force->SetLumenBoundaryAdhesionEnergyParameter(7.0);
        p_force->SetEpiBoundaryAdhesionEnergyParameter(M_EPIBND);

        simulator.AddForce(p_force);

        /*MAKE_PTR(TargetAreaModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);
        p_growth_modifier->SetReferenceTargetArea(0.9);*/
        MAKE_PTR(VolumeTrackingModifier<2>, p_volume_modifier);
        simulator.AddSimulationModifier(p_volume_modifier);
        MAKE_PTR(BorderTrackingModifier<2>, p_border_modifier);
        simulator.AddSimulationModifier(p_border_modifier);
        MAKE_PTR(PerimeterTrackingModifier<2>, p_stretch_modifier);
        simulator.AddSimulationModifier(p_stretch_modifier);
        //MAKE_PTR(MassCenterTrackingModifier<2>, p_center_modifier);
        //simulator.AddSimulationModifier(p_center_modifier) ;
        MAKE_PTR(ObstructionWriterModifier<2>, p_endopos_modifier);
        simulator.AddSimulationModifier(p_endopos_modifier);


        std::cout << "Labelling Epi cells" << endl ;
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
          //FOR LUMEN
          cell_iter->GetCellData()->SetItem("cellIndex",SimulationParameters::getNextIndex());
          cell_iter->GetCellData()->SetItem("timeFromLastLumenGeneration",0);
          cell_iter->GetCellData()->SetItem("hadALumenDivision",0);
          cell_iter->GetCellData()->SetItem("lumenNearby",1);
          cell_iter->GetCellData()->SetItem("vecPolaX",0);
          cell_iter->GetCellData()->SetItem("vecPolaY",0);


          if (cell_population.GetLocationIndexUsingCell(*cell_iter) == 188 || cell_population.GetLocationIndexUsingCell(*cell_iter) == 239 || cell_population.GetLocationIndexUsingCell(*cell_iter) == 160 || cell_population.GetLocationIndexUsingCell(*cell_iter) == 35 )
          {
              cell_iter->AddCellProperty(p_endo);
              cell_iter->AddCellProperty(p_stalk);
              //cell_iter->AddCellProperty(p_periph); // PLANTE LA SIMU

              cell_iter->GetCellData()->SetItem("tagVessel",-1);
          }
          /*else if (cell_population.GetLocationIndexUsingCell(*cell_iter) == 212 || cell_population.GetLocationIndexUsingCell(*cell_iter) == 109 || cell_population.GetLocationIndexUsingCell(*cell_iter) == 103 || cell_population.GetLocationIndexUsingCell(*cell_iter) == 81 )
          {
              cell_iter->AddCellProperty(p_endo);
              cell_iter->AddCellProperty(p_tip);

              cell_iter->GetCellData()->SetItem("tagVessel",-1);
          }*/
          else
          {
              cell_iter->AddCellProperty(p_epi);
          }
        }



        // NE PAS DECOMMENTER LA SECTION SUIVANTE (bugs à régler)
        /*
        std::cout << "Adding active force" << endl ;
        MAKE_PTR_ARGS(MorphogenDrivenCellForce<2>, p_motile_force, (15,0.55));
        simulator.AddForce(p_motile_force);

        MAKE_PTR(CellFixingModifier<2>, p_fixing_modifier);
        simulator.AddSimulationModifier(p_fixing_modifier);
        */
/*
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

        std::cout << "Adding Gradient Modifier" << endl ;

        p_pde_modifier->SetOutputGradient(true);
        p_pde_modifier->GetOutputGradient();

        simulator.AddSimulationModifier(p_pde_modifier);
        */

        //MAKE_PTR(MorphogenTrackingModifier<2>, p_morphoendo_modifier);
        //simulator.AddSimulationModifier(p_morphoendo_modifier);



        std::cout << "Adding repulsion force" << endl ;
        //MAKE_PTR_ARGS(RepulsionForce<2>, p_repulsion_force, (0.5));
        //simulator.AddForce(p_repulsion_force);


        MAKE_PTR_ARGS(ObstructionBoundaryCondition<2>, p_repuls_bc, (&cell_population));
        simulator.AddCellPopulationBoundaryCondition(p_repuls_bc);



        std::cout << "Adding lumen" << endl ;
        MAKE_PTR(SimuInfoModifier<2>, p_simuInfoModifier);
        simulator.AddSimulationModifier(p_simuInfoModifier);

        MAKE_PTR(PolarisationModifier<2>, p_polarisation_modifier);
        simulator.AddSimulationModifier(p_polarisation_modifier);
        p_polarisation_modifier->SetEpiEpiPolarisationParameter(M_EPIEPI_INI*5.0);
        p_polarisation_modifier->SetEndoEpiPolarisationParameter(M_ENDOEPI_INI*5.0);
        p_polarisation_modifier->SetLumenEpiPolarisationParameter(M_LUMENEPI_INI);
        p_polarisation_modifier->SetVecPolarisationDecrease(M_POLARDEC_INI);

        //MAKE_PTR(NewEndoGeneratorModifier<2>, p_newendo_modifier);
        //simulator.AddSimulationModifier(p_newendo_modifier);

        MAKE_PTR_ARGS(FixedBoundaryCondition<2>, p_fixed_bc, (&cell_population));
        simulator.AddCellPopulationBoundaryCondition(p_fixed_bc);

        //MAKE_PTR(LumenGenerationModifier<2>, p_lumen_generation_modifier);
        //simulator.AddSimulationModifier(p_lumen_generation_modifier);
        //MAKE_PTR(DifferentialTargetAreaModifier<2>, p_growth_modifier);
        //simulator.AddSimulationModifier(p_growth_modifier);

        MAKE_PTR(LumenModifier<2>, p_lumen_modifier);
        simulator.AddSimulationModifier(p_lumen_modifier);
        p_lumen_modifier->SetLumenSizeFactor(M_LUMEN_SIZE_FACTOR_INI) ;
        p_lumen_modifier->SetlumenDuration2TargetArea(M_DURATION2_INI) ;

/*
        std::cout << "Adding active force" << endl ;
        MAKE_PTR_ARGS(MorphogenCellForce<2>, p_motile_force, (5.0));//force initiale, croissante
        simulator.AddForce(p_motile_force);
*/

        std::cout << "Growing Monolayer" << endl ;

        simulator.SetEndTime(96.0);
        simulator.SetDt(0.01);
        simulator.SetSamplingTimestepMultiple(1.0);

        simulator.Solve();

        // Record mesh
      //  VertexMeshWriter<2,2> vertex_mesh_writer("TestRepulsionForcenMeshWriter", "repulsion_mesh");
      //  vertex_mesh_writer.WriteFilesUsingMesh(*p_mesh);

    }
};

#endif /* TESTVERTEXBASEDFORTESTING_HPP_ */
