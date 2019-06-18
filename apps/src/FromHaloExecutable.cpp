#include <iomanip>

// Includes from trunk
#include "ExecutableSupport.hpp"
#include "CellBasedSimulationArchiver.hpp"

#include "SmartPointers.hpp"
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
#include "CellAdjacencyMatrixWriter.hpp"

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

#include "FixedBoundaryCondition.hpp"
#include "PerimeterDependentCellCycleModel.hpp"
#include "StochasticLumenCellCycleModel.hpp"
#include "SimplePositionBasedCellCycleModel.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
using namespace std ;

// Program option includes for handling command line arguments
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/lexical_cast.hpp>

/*
 * Prototype functions
 */

void SetupSingletons();
void SetupAndRunSimulation(unsigned mEpiEpi, unsigned mEpiBnd, unsigned mEpiEndo, unsigned mMemb);
void DestroySingletons();

int main(int argc, char *argv[])
{
    // This sets up PETSc and prints out copyright information, etc.
    ExecutableSupport::StartupWithoutShowingCopyright(&argc, &argv);

    // Define command line options
    boost::program_options::options_description general_options("This is the adhesion calibration executable.\n");
    general_options.add_options()
                    ("help", "produce help message")
                    ("E", boost::program_options::value<unsigned>()->default_value(5),"The energy parameter for the Epi-Epi adhesion")
                    ("B", boost::program_options::value<unsigned>()->default_value(5),"The energy parameter for the Epi-Boundary adhesion")
                    ("D", boost::program_options::value<unsigned>()->default_value(10),"The energy parameter for the Epi-Endo adhesion")
                    ("M", boost::program_options::value<unsigned>()->default_value(1),"The energy parameter for the Membrane Surface");

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
    unsigned epiepi = variables_map["E"].as<unsigned>();
    unsigned epibnd = variables_map["B"].as<unsigned>();
    unsigned epiendo = variables_map["D"].as<unsigned>();
    unsigned memb = variables_map["M"].as<unsigned>();

    SetupSingletons();
    SetupAndRunSimulation(epiepi, epibnd, epiendo, memb);
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

void SetupAndRunSimulation(unsigned mEpiEpi, unsigned mEpiBnd, unsigned mEpiEndo, unsigned mMemb)
{
        // Specify output directory (unique to each simulation)
        std::string output_directory = std::string("FromHaloCalibration")
        + std::string("_E") + boost::lexical_cast<std::string>(mEpiEpi)
        + std::string("_B") + boost::lexical_cast<std::string>(mEpiBnd)
        + std::string("_D") + boost::lexical_cast<std::string>(mEpiEndo)
        + std::string("_M") + boost::lexical_cast<std::string>(mMemb);

        // Permet d'importer un fichier test et s'en servir pour taguer les cellules, voir ci-après
        std::cout << "Importing label data from txt" << std::endl ;
        ifstream inFile ;
        int x ;
        inFile.open("/home/gonayl/Chaste/chaste-src/testoutput/test_label.txt") ;
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


        // Create Mesh

        std::cout << "Creating mesh" << endl ;

        VertexMeshReader<2,2> mesh_reader("/home/gonayl/Chaste/chaste-src/testoutput/mesh/vertex_based_mesh");
        MutableVertexMesh<2,2> p_mesh;
        p_mesh.ConstructFromMeshReader(mesh_reader);
        p_mesh.SetCellRearrangementThreshold(0.1);
        unsigned num_cells = p_mesh.GetNumElements();

        // Create Cells
        std::vector<CellPtr> cells;
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
        MAKE_PTR(WildTypeCellMutationState, p_state);

        // ICI : modifier durée du cycle cellulaire MAIS dans code source (cfr UniformG1UniformG1GenerationalBoundaryCellCycleModel.cpp)!

        for (unsigned i=0; i<num_cells; i++)
        {
            UniformG1GenerationalBoundaryCellCycleModel* p_cycle_model = new UniformG1GenerationalBoundaryCellCycleModel();
          if (label_input[i] == 0)
          {
            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            double birth_time = rand() % 10 + 1 ;
            p_cell->SetBirthTime(-birth_time);
            p_cell->InitialiseCellCycleModel();
            p_cell->GetCellData()->SetItem("target area", 0.5);
            cells.push_back(p_cell);
          }
          else if (label_input[i] == 1)
          {
            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            double birth_time = rand() % 10 + 1 ;
            p_cell->SetBirthTime(-birth_time);
            p_cell->InitialiseCellCycleModel();
            p_cell->GetCellData()->SetItem("target area", 0.5);
            cells.push_back(p_cell);
          }
        }



        VertexBasedCellPopulation<2> cell_population(p_mesh, cells);

        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellPosWriter>();
        cell_population.AddCellWriter<CellTypeWriter>();
        cell_population.AddCellWriter<CellVolumesWriter>();
        cell_population.AddPopulationWriter<CellAdjacencyMatrixWriter>();

        MAKE_PTR(CellEpi, p_epi);
        MAKE_PTR(CellTip, p_tip);
        MAKE_PTR(CellEndo, p_endo);
        MAKE_PTR(CellLabel, p_label);
        MAKE_PTR(CellStalk, p_stalk);

        // Create Simulation
        OffLatticeSimulation<2> simulator(cell_population);

        simulator.SetOutputDivisionLocations(true);

        // ICI : MODIFIER PARAMETRE D'ADHESION
        std::cout << "Adding passive force" << endl ;
        // Create Forces and pass to simulation
        MAKE_PTR(DifferentialAdhesionForce<2>, p_force);
        p_force->SetNagaiHondaDeformationEnergyParameter(55.0);
        p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(mMemb);

        p_force->SetEndoEndoAdhesionEnergyParameter(3.0);
        p_force->SetLumenLumenAdhesionEnergyParameter(3.0);
        p_force->SetCoreCoreAdhesionEnergyParameter(mEpiEpi);
        p_force->SetCorePeriphAdhesionEnergyParameter(mEpiEpi);
        p_force->SetPeriphPeriphAdhesionEnergyParameter(mEpiEpi);
        p_force->SetEndoEpiAdhesionEnergyParameter(mEpiEndo);
        p_force->SetEpiLumenAdhesionEnergyParameter(6.0);
        p_force->SetEndoLumenAdhesionEnergyParameter(6.0);

        p_force->SetNagaiHondaCellBoundaryAdhesionEnergyParameter(10.0);
        p_force->SetEndoBoundaryAdhesionEnergyParameter(10.0);
        p_force->SetLumenBoundaryAdhesionEnergyParameter(10.0);
        p_force->SetEpiBoundaryAdhesionEnergyParameter(mEpiBnd);

        simulator.AddForce(p_force);

        MAKE_PTR(VolumeTrackingModifier<2>, p_volume_modifier);
        simulator.AddSimulationModifier(p_volume_modifier);
        MAKE_PTR(PerimeterTrackingModifier<2>, p_stretch_modifier);
        simulator.AddSimulationModifier(p_stretch_modifier);
        MAKE_PTR(BorderTrackingModifier<2>, p_border_modifier);
        simulator.AddSimulationModifier(p_border_modifier);
        //MAKE_PTR(NeighbourTrackingModifier<2>, p_neighbour_modifier) ;
        //simulator.AddSimulationModifier(p_neighbour_modifier) ;
        //MAKE_PTR(PolarTrackingModifier<2>, p_polar_modifier) ;
        //simulator.AddSimulationModifier(p_polar_modifier) ;
        //MAKE_PTR(LabelTrackingModifier<2>, p_lumen_modifier) ;
        //simulator.AddSimulationModifier(p_lumen_modifier) ;

        std::cout << "Labelling Epi cells" << endl ;
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
        {
          if (label_input[cell_population.GetLocationIndexUsingCell(*cell_iter)] == 0)
          {
              cell_iter->AddCellProperty(p_epi);
          }
          else
          {
              cell_iter->AddCellProperty(p_endo);
              cell_iter->AddCellProperty(p_stalk);
          }
        }

        MAKE_PTR_ARGS(FixedBoundaryCondition<2>, p_fixed_bc, (&cell_population));
        simulator.AddCellPopulationBoundaryCondition(p_fixed_bc);

        std::cout << "Growing Monolayer" << endl ;

        simulator.SetEndTime(48.0);
        simulator.SetDt(1.0/20.0);
        simulator.SetSamplingTimestepMultiple(1.0);
        simulator.SetOutputDirectory(output_directory);

        cout << "Testing parameters set : (" << mEpiEpi << " ; " << mEpiBnd << " ; " << mEpiEndo << " ; " << mMemb << " )" << endl;
        simulator.Solve();

}
