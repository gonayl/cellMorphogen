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
#include "RepulsionForce.hpp"

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


#include "MorphogenCellForce.hpp"
#include "MassCenterTrackingModifier.hpp"

//FOR LUMEN
#include "SimulationParameters.hpp"

#include "LumenGenerationModifier.hpp"
#include "LumenModifier.hpp"
#include "PolarisationModifier.hpp"
#include "SimuInfoModifier.hpp"

#include "CellPolarityXWriter.hpp"
#include "CellPolarityYWriter.hpp"

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
static const double M_EPIBND = 5.0 ;
static const double M_ENDOBND = 4.0 ;//chang
static const double M_ENDOEPI = 5.0 ;
static const double M_MOTILITY = 15.0 ;


// Program option includes for handling command line arguments
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/lexical_cast.hpp>

/*
 * Prototype functions
 */

void SetupSingletons();
void SetupAndRunSimulation(unsigned mScaling, out_stream overall_results_file);
void DestroySingletons();

int main(int argc, char *argv[])
{
    // This sets up PETSc and prints out copyright information, etc.
    ExecutableSupport::StartupWithoutShowingCopyright(&argc, &argv);

    // Define command line options
    boost::program_options::options_description general_options("This is the lumen calibration executable.\n");
    general_options.add_options()
                    ("help", "produce help message")
                    ("X", boost::program_options::value<unsigned>()->default_value(1),"scaling factor");

    // Define parse command line into variables_map
    boost::program_options::variables_map variables_map;
    boost::program_options::store(parse_command_line(argc, argv, general_options), variables_map);

    // Print help message if wanted
    if (variables_map.count("help"))
    {
        std::cout << std::setprecision(3) << general_options << "\n";
        std::cout << general_options << "\n";
        return 1;
    }

    // Get ID and name from command line
    unsigned scaling = variables_map["X"].as<unsigned>();

    SetupSingletons();
    // Create results file handler
    OutputFileHandler results_handler("TestLumenCalibraiton/lumen_calibration_results", false);
    // Create overall results file
    std::string overall_results_filename = "overall_results" + "_x" + boost::lexical_cast<std::string>(scaling) + ".dat";
    out_stream overall_results_file = results_handler.OpenOutputFile(overall_results_filename);

    SetupAndRunSimulation(scaling, overall_results_file);

    *overall_results_file << "SIMULATIONS COMPLETE\n" << std::flush;
    overall_results_file->close();
    DestroySingletons();
}

void SetupSingletons()
{
    SimulationTime::Instance()->SetStartTime(0.0);
    CellPropertyRegistry::Instance()->Clear();
}

void DestroySingletons()
{
    SimulationTime::Destroy();
    CellPropertyRegistry::Instance()->Clear();
}



void SetupAndRunSimulation(unsigned mScaling, out_stream overall_results_file)
{

  // Permet d'importer un fichier test et s'en servir pour taguer les cellules, voir ci-après
  std::cout << "Importing label data from txt" << std::endl ;
  ifstream inFile ;
  int x ;
  //inFile.open("testoutput/SimpleConditionInit/test_label_simple.txt") ;
  inFile.open("testoutput/test_label_simple.txt") ;
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
  //inFileBnd.open("testoutput/SimpleConditionInit/boundary_input.txt") ;
  inFileBnd.open("testoutput/boundary_input.txt") ;
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
        VertexMeshReader<2,2> mesh_reader("/home/gonayl/Chaste/chaste-src/testoutput/mesh/vertex_based_mesh");
        MutableVertexMesh<2,2> p_mesh;
        p_mesh.ConstructFromMeshReader(mesh_reader);
        p_mesh.SetCellRearrangementThreshold(0.1);

        unsigned num_cells = p_mesh.GetNumElements();

        std::cout << "mesh generated" << endl ;

        // Creat cells

        std::vector<CellPtr> cells;

        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        MAKE_PTR(WildTypeCellMutationState, p_state);

        for (unsigned i=0; i<num_cells; i++)
        {
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



        VertexBasedCellPopulation<2> cell_population(p_mesh, cells);

        std::cout << "cells generated" << endl ;

        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellPosWriter>();
        cell_population.AddCellWriter<CellTypeWriter>();
        //cell_population.AddCellWriter<CellVolumesWriter>();                   COMMENTE PAR MOI
        //cell_population.AddPopulationWriter<CellAdjacencyMatrixWriter>();
        //cell_population.AddPopulationWriter<CalibrationErrorWriter>();        COMMENTE PAR MOI


        //FOR LUMEN
        cell_population.AddCellWriter<CellPolarityXWriter>();
        cell_population.AddCellWriter<CellPolarityYWriter>();


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

        p_force->SetEndoEndoAdhesionEnergyParameter(1.7);
        p_force->SetLumenLumenAdhesionEnergyParameter(5.0);
        p_force->SetCoreCoreAdhesionEnergyParameter(M_EPI);
        p_force->SetCorePeriphAdhesionEnergyParameter(M_EPI);
        p_force->SetPeriphPeriphAdhesionEnergyParameter(M_EPI);
        p_force->SetEndoEpiAdhesionEnergyParameter(M_ENDOEPI);
        p_force->SetEpiLumenAdhesionEnergyParameter(4.0);
        p_force->SetEndoLumenAdhesionEnergyParameter(35.0);

        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(10.0);
        p_force->SetEndoBoundaryAdhesionEnergyParameter(M_ENDOBND);
        p_force->SetLumenBoundaryAdhesionEnergyParameter(5.0);
        p_force->SetEpiBoundaryAdhesionEnergyParameter(M_EPIBND);


        p_force->SetEpiEpiAdhesionEnergyParameter(6);

        simulator.AddForce(p_force);

        MAKE_PTR(VolumeTrackingModifier<2>, p_volume_modifier);
        simulator.AddSimulationModifier(p_volume_modifier);
        MAKE_PTR(PerimeterTrackingModifier<2>, p_stretch_modifier);
        simulator.AddSimulationModifier(p_stretch_modifier);
        MAKE_PTR(BorderTrackingModifier<2>, p_border_modifier);
        simulator.AddSimulationModifier(p_border_modifier);
        MAKE_PTR(MassCenterTrackingModifier<2>, p_center_modifier);
        simulator.AddSimulationModifier(p_center_modifier) ;

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


          if (label_input[cell_population.GetLocationIndexUsingCell(*cell_iter)] == 0)
          {
              cell_iter->AddCellProperty(p_epi);
          }
          else if (label_input[cell_population.GetLocationIndexUsingCell(*cell_iter)] == 1)
          {
              cell_iter->AddCellProperty(p_endo);
              cell_iter->AddCellProperty(p_stalk);
              cell_iter->GetCellData()->SetItem("tagVessel",-1);
          }
          else if (label_input[cell_population.GetLocationIndexUsingCell(*cell_iter)] == 2)
          {
              cell_iter->AddCellProperty(p_endo);
              cell_iter->AddCellProperty(p_tip);
              cell_iter->GetCellData()->SetItem("tagVessel",cell_iter->GetCellData()->GetItem("cellIndex"));

          }
          else if (label_input[cell_population.GetLocationIndexUsingCell(*cell_iter)] == 3)
          {
              cell_iter->AddCellProperty(p_endo);
              cell_iter->AddCellProperty(p_vessel);

          }
        }


        //FOR LUMEN
        std::cout << "Adding lumen" << endl ;
        MAKE_PTR(SimuInfoModifier<2>, p_simuInfoModifier);
        simulator.AddSimulationModifier(p_simuInfoModifier);
        MAKE_PTR(PolarisationModifier<2>, p_polarisation_modifier);
        simulator.AddSimulationModifier(p_polarisation_modifier);

        MAKE_PTR(LumenGenerationModifier<2>, p_lumen_generation_modifier);
        simulator.AddSimulationModifier(p_lumen_generation_modifier);

        MAKE_PTR(LumenModifier<2>, p_lumen_modifier);
        simulator.AddSimulationModifier(p_lumen_modifier);


        //MAKE_PTR(MorphogenTrackingModifier<2>, morphogen_modifier);
        //simulator.AddSimulationModifier(morphogen_modifier);


        // NE PAS DECOMMENTER LA SECTION SUIVANTE (bugs à régler)

        std::cout << "Adding active force" << endl ;
        MAKE_PTR_ARGS(MorphogenCellForce<2>, p_motile_force, (7.2));//force initiale, decroissement
        simulator.AddForce(p_motile_force);


        std::cout << "Growing Monolayer" << endl ;

        simulator.SetEndTime(SimulationParameters::TIME_OF_SIMULATION);
        simulator.SetDt(SimulationTime::Instance()->GetTimeStep());
        simulator.SetSamplingTimestepMultiple(1);


        /* Calibration part */


        *overall_results_file << "Simulation with scaling x " << mScaling << std::flush;
        cout << "Testing scaling x" << mScaling << endl;

        // Specify output directory (unique to each simulation)
        std::string output_directory = std::string("TestLumenCalibraiton/")
        + std::string("_X") + boost::lexical_cast<std::string>(mScaling);

        /* END of Calibration part */
        simulator.SetOutputDirectory(output_directory);

        simulator.Solve();


}
