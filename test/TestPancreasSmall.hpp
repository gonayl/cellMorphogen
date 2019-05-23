#ifndef TESTVERTEXBASEDTESTING_HPP_
#define TESTVERTEXBASEDTESTING_HPP_

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
// #include "MergeNodeModifier.hpp"


#include "PerimeterDependentCellCycleModel.hpp"
#include "StochasticLumenCellCycleModel.hpp"
#include "SimplePositionBasedCellCycleModel.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
using namespace std ;

static const double M_TIME_FOR_SIMULATION = 40;
static const double M_NUM_CELLS_ACROSS = 10;
static const double M_UPTAKE_RATE = 10.0 ;
static const double M_DIFFUSION_CONSTANT = 5e-1;
static const double M_DUDT_COEFFICIENT = 1.0;
static const double M_DECAY_COEFFICIENT = 9.0;
static const double M_RADIUS = 100.0;

class TestCellMorphogen : public AbstractCellBasedWithTimingsTestSuite
{
private:

    void GenerateCells(unsigned num_cells, std::vector<CellPtr>& rCells)
    {
        MAKE_PTR(TransitCellProliferativeType, p_transit_type);
        MAKE_PTR(WildTypeCellMutationState, p_state);

        for (unsigned i=0; i<num_cells; i++)
        {
            //StochasticLumenCellCycleModel* p_cycle_model = new StochasticLumenCellCycleModel();
            //UniformCellCycleModel* p_cycle_model = new UniformCellCycleModel();
            UniformG1GenerationalBoundaryCellCycleModel* p_cycle_model = new UniformG1GenerationalBoundaryCellCycleModel();

            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            // double birth_time = - RandomNumberGenerator::Instance()->ranf() * p_cycle_model->GetAverageTransitCellCycleTime();
            double birth_time = 0 ;
            p_cell->SetBirthTime(-birth_time);
            p_cell->InitialiseCellCycleModel();
            p_cell->GetCellData()->SetItem("target area", 0.5);
            rCells.push_back(p_cell);
        }
     }




public:

  void TestNodeBasedMeshReader()
  {
    // Create a simple mesh
    /*
    HoneycombMeshGenerator generator(1, 4);
    MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
    NodesOnlyMesh<2> p_mesh;
    p_mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);
    std::vector<CellPtr> cells;
    GenerateCells(p_mesh.GetNumNodes(),cells);
    // Create cell population
    NodeBasedCellPopulation<2> cell_population(p_mesh, cells);
    */
    TrianglesMeshReader<2,2> mesh_reader("testoutput/mesh/node_based_mesh");
    TetrahedralMesh<2,2> p_generating_mesh ;
    p_generating_mesh.ConstructFromMeshReader(mesh_reader);
    //Extended to allow sorting for longer distances
    double cut_off_length = 2.5;
    // Convert this to a NodesOnlyMesh
    NodesOnlyMesh<2> p_mesh;
    p_mesh.ConstructNodesWithoutMesh(p_generating_mesh, cut_off_length);
    std::vector<CellPtr> cells;
    GenerateCells(p_mesh.GetNumNodes(),cells);

    // Create cell population
    NodeBasedCellPopulation<2> cell_population(p_mesh, cells);

    //Make cell data writer so can pass in variable name
    cell_population.AddCellWriter<CellAgesWriter>();
    // cell_population.AddCellWriter<CellAppliedForceWriter>();
    // MAKE_PTR(ForceTrackingModifier<2>, force_modifier);


    // Set up cell-based simulation and output directory
    OffLatticeSimulation<2> simulator(cell_population);
    simulator.SetOutputDirectory("CellMorphogen/OSModel/TestMeshReader/Thy155");

    // Create a force law and pass it to the simulation
    MAKE_PTR(DifferentialAdhesionGeneralisedLinearSpringForce<2>, p_differential_adhesion_force);
    p_differential_adhesion_force->SetMeinekeSpringStiffness(15.0);
    p_differential_adhesion_force->SetHomotypicLabelledSpringConstantMultiplier(1.0);
    p_differential_adhesion_force->SetHeterotypicSpringConstantMultiplier(0.5);
    p_differential_adhesion_force->SetCutOffLength(2.5);
    simulator.AddForce(p_differential_adhesion_force);


    cell_population.AddCellWriter<CellAgesWriter>();
    cell_population.AddCellWriter<CellPosWriter>();
    cell_population.AddCellWriter<CellTypeWriter>();
    cell_population.AddCellWriter<CellVolumesWriter>();

    MAKE_PTR(CellEpi, p_epi);
    MAKE_PTR(CellEndo, p_endo);

    std::cout << "Importing label data from txt" << std::endl ;
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
       }

    }

    MAKE_PTR(VolumeTrackingModifier<2>, p_volume_modifier);
    simulator.AddSimulationModifier(p_volume_modifier);
    MAKE_PTR(TargetAreaModifier<2>, p_growth_modifier);
    simulator.AddSimulationModifier(p_growth_modifier);
    p_growth_modifier->SetReferenceTargetArea(0.8);


    /* simulator.AddSimulationModifier(force_modifier);

    MAKE_PTR(PolarityTrackingModifier<2>, polarity_modifier);
    simulator.AddSimulationModifier(polarity_modifier);

    MAKE_PTR_ARGS(MyForce, p_force, (3.0));
    simulator.AddForce(p_force);
    std::cout << "Active force" << endl ; */


    // FileFinder file_finder("blah", RelativeTo::ChasteTestOutput);
    // std::cout << file_finder.GetAbsolutePath();

    simulator.SetEndTime(1.0);
    simulator.SetDt(1.0/100.0);
    simulator.SetSamplingTimestepMultiple(1);
    std::cout << "Growing Monolayer" << endl ;
    simulator.Solve();

    /*
    // Record mesh
    TrianglesMeshWriter<2,2> triangles_mesh_writer("TestNodeBasedMeshWriter", "node_based_mesh");
    triangles_mesh_writer.WriteFilesUsingMesh(p_mesh);
    */

  }

  void TestVertexBasedMeshReader()
  {
      // Create Mesh

      std::cout << "Creating mesh" << endl ;

      VertexMeshReader<2,2> mesh_reader("testoutput/mesh/vertex_based_mesh");
      MutableVertexMesh<2,2> p_mesh;
      p_mesh.ConstructFromMeshReader(mesh_reader);
      p_mesh.SetCellRearrangementThreshold(0.1);

      // Create Cells
      std::vector<CellPtr> cells;
      GenerateCells(p_mesh.GetNumElements(),cells);

      VertexBasedCellPopulation<2> cell_population(p_mesh, cells);

      cell_population.AddCellWriter<CellAgesWriter>();
      cell_population.AddCellWriter<CellPosWriter>();
      cell_population.AddCellWriter<CellTypeWriter>();

      MAKE_PTR(CellEpi, p_epi);
      MAKE_PTR(CellEndo, p_endo);

      // Create Simulation
      OffLatticeSimulation<2> simulator(cell_population);
      simulator.SetOutputDirectory("CellMorphogen/VertexModel/TestMeshReader/E135A1Coupe2_2");


      simulator.SetOutputDivisionLocations(true);

      std::cout << "Adding passive force" << endl ;
      // Create Forces and pass to simulation
      MAKE_PTR(DifferentialAdhesionForce<2>, p_force);
      p_force->SetNagaiHondaDeformationEnergyParameter(55.0);
      p_force->SetNagaiHondaMembraneSurfaceEnergyParameter(1.0);

      p_force->SetEndoEndoAdhesionEnergyParameter(5.0);
      p_force->SetLumenLumenAdhesionEnergyParameter(5.0);
      p_force->SetCoreCoreAdhesionEnergyParameter(5.0);
      p_force->SetCorePeriphAdhesionEnergyParameter(5.0);
      p_force->SetPeriphPeriphAdhesionEnergyParameter(5.0);
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
      p_growth_modifier->SetReferenceTargetArea(0.9);
      MAKE_PTR(VolumeTrackingModifier<2>, p_volume_modifier);
      simulator.AddSimulationModifier(p_volume_modifier);

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


      // Labeling all cells are epithelial


      /*std::cout << "Labelling Epi cells" << endl ;

      for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
           cell_iter != cell_population.End();
           ++cell_iter)
      {

        cell_iter->AddCellProperty(p_epi);

      }*/

      // Labelling cells from file (Epi/Endo/Lumen)

      std::cout << "Importing label data from txt" << std::endl ;
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
         }

      }

      MAKE_PTR(MorphogenTrackingModifier<2>, morphogen_modifier);
      simulator.AddSimulationModifier(morphogen_modifier);

      std::cout << "Adding active force" << endl ;
      MAKE_PTR_ARGS(MorphogenDrivenCellForce<2>, p_motile_force, (15,0.55));
      simulator.AddForce(p_motile_force);

      std::cout << "Growing Monolayer" << endl ;

      simulator.SetEndTime(360.0);
      simulator.SetDt(1.0/100.0);
      simulator.SetSamplingTimestepMultiple(1.0);

      simulator.Solve();

  }

};

#endif /* TESTVERTEXBASEDFORTESTING_HPP_ */
