#include <cxxtest/TestSuite.h>
#include "CellBasedSimulationArchiver.hpp"
#include "SmartPointers.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "AdhesionPottsUpdateRule.hpp"
#include "CellsGenerator.hpp"
#include "CylindricalHoneycombMeshGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "MeshBasedCellPopulationWithGhostNodes.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "OnLatticeSimulation.hpp"
#include "PlaneBoundaryCondition.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "PottsMeshGenerator.hpp"
#include "RandomCellKiller.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "SurfaceAreaConstraintPottsUpdateRule.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "VolumeConstraintPottsUpdateRule.hpp"
#include "VoronoiDataWriter.hpp"
#include "DifferentialAdhesionGeneralisedLinearSpringForce.hpp"
#include "TrianglesMeshWriter.hpp"
#include "CellAgesWriter.hpp"
#include "CellAppliedForceWriter.hpp"
#include "FakePetscSetup.hpp"
#include "CounterSingleton.hpp"

#include "CellDataItemWriter.hpp"

#include "ParabolicGrowingDomainPdeModifier.hpp"
#include "MorphogenCellwiseSourceParabolicPde.hpp"
#include "PolarityTrackingModifier.hpp"
#include "CellAddingModifier.hpp"

#include "AbstractForce.hpp"
#include <numeric>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <complex>      // std::complex, std::polar

#include "CellLabel.hpp"
#include "CellEndo.hpp"
#include "CellLabelWriter.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "MorphogenTrackingModifier.hpp"
#include "PolarityTrackingModifier.hpp"
#include "PerimeterTrackingModifier.hpp"
#include "ForceTrackingModifier.hpp"
#include "PositionTrackingModifier.hpp"
#include "RadiusTrackingModifier.hpp"

#include "PerimeterDependentCellCycleModel.hpp"
#include "ForceDependentCellCycleModel.hpp"
#include "PositionDependentCellCycleModel.hpp"
#include "Debug.hpp"
// #include "EndothelialAdhesionGeneralisedLinearSpringForce.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CounterSingleton.hpp"
#include "DifferentialAdhesionForce.hpp"

class MyForce : public AbstractForce<2>
{
private:

    double mStrength;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<2> >(*this);
        archive & mStrength;
    }

public:
    MyForce(double strength=1.0)
        : AbstractForce<2>(),
          mStrength(strength)
    {
        assert(mStrength > 0.0);
    }

    void AddForceContribution(AbstractCellPopulation<2>& rCellPopulation)
    {


        c_vector<double, 2> force = zero_vector<double>(2);
        force(1) = 0.0 ;
        force(0) = 0.0 ;

        for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
                cell_iter != rCellPopulation.End();
                ++cell_iter)
           {
             if (cell_iter->HasCellProperty<CellLabel>())
             {
             // double polarity_x = cell_iter->GetCellData()->GetItem("polarityx");
             // double polarity_y = cell_iter->GetCellData()->GetItem("polarityy");
             // force(1) =    mStrength * polarity_x ;
             // force(0) =    mStrength * polarity_y ;

             force(1) =    0.0 ;
             force(0) =    mStrength  ;

             unsigned node_index = rCellPopulation.GetLocationIndexUsingCell(*cell_iter);
             rCellPopulation.GetNode(node_index)->AddAppliedForceContribution(force);
             }
           }

    }

    double GetStrength()
    {
        return mStrength;
    }

    void OutputForceParameters(out_stream& rParamsFile)
    {
        *rParamsFile << "\t\t\t<Strength>" << mStrength << "</Strength>\n";
        AbstractForce<2>::OutputForceParameters(rParamsFile);
    }
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(MyForce)
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(MyForce)



class TestCenterBasedModel : public AbstractCellBasedTestSuite
{

private:

  void GenerateCells(unsigned num_cells, std::vector<CellPtr>& cells)
  {

      // RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

      for (unsigned i=0; i<num_cells; i++)
      {
          //UniformlyDistributedCellCycleModel* p_cycle_model = new UniformlyDistributedCellCycleModel();
          //MorphogenDependentCellCycleModel* p_cycle_model = new MorphogenDependentCellCycleModel();
          UniformG1GenerationalCellCycleModel* p_cycle_model = new UniformG1GenerationalCellCycleModel();
          MAKE_PTR(WildTypeCellMutationState, p_state);
          MAKE_PTR(TransitCellProliferativeType, p_transit_type);
          MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
          // double phase = rand() % 20 + 1;
          // unsigned node_index = i;
          // double magnitude = 1.0 ;

          /* std::complex<double> coord = std::polar (magnitude,phase);
          //std::cout << node_index << endl ;
          double polarityx = std::abs(coord)*cos(std::arg(coord));
          double polarityy = std::abs(coord)*sin(std::arg(coord));

          std::cout << coord << "corresponds to (" << polarityx << "," << polarityy << ")" ; */

          CellPtr p_cell(new Cell(p_state, p_cycle_model));
          p_cell->SetCellProliferativeType(p_transit_type);
          /* p_cell->GetCellData()->SetItem("magnitude", magnitude);
          p_cell->GetCellData()->SetItem("phase", phase);
          p_cell->GetCellData()->SetItem("polarityx", polarityx);
          p_cell->GetCellData()->SetItem("polarityy", polarityy); */
          p_cell->GetCellData()->SetItem("target area", 1.0);
          cells.push_back(p_cell);

      }
   }

public:
  void TestNodeBasedMonolayer()
  {
        // Create a simple mesh

        HoneycombMeshGenerator generator(1, 4);
        MutableMesh<2,2>* p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2> p_mesh;
        p_mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1.5);
        std::vector<CellPtr> cells;
        GenerateCells(p_mesh.GetNumNodes(),cells);
        // Create cell population
        NodeBasedCellPopulation<2> cell_population(p_mesh, cells);
        /*
        TrianglesMeshReader<2,2> mesh_reader("testoutput/TestNodeBasedMeshWriter/node_based_mesh");
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
        */

        //Make cell data writer so can pass in variable name
        cell_population.AddCellWriter<CellAgesWriter>();
        // cell_population.AddCellWriter<CellAppliedForceWriter>();
        // MAKE_PTR(ForceTrackingModifier<2>, force_modifier);


        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestCenterBasedModel/ForTesting/Remesh");

        // Create a force law and pass it to the simulation
        MAKE_PTR(DifferentialAdhesionGeneralisedLinearSpringForce<2>, p_differential_adhesion_force);
        p_differential_adhesion_force->SetMeinekeSpringStiffness(15.0);
        p_differential_adhesion_force->SetHomotypicLabelledSpringConstantMultiplier(5.0);
        p_differential_adhesion_force->SetHeterotypicSpringConstantMultiplier(0.5);
        p_differential_adhesion_force->SetCutOffLength(2.5);
        simulator.AddForce(p_differential_adhesion_force);

        MAKE_PTR(CellLabel, p_label);
        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)

             {
              if (cell_population.GetLocationIndexUsingCell(*cell_iter) == 1)
              {
                  cell_iter->AddCellProperty(p_label);
              }
             }

        /* simulator.AddSimulationModifier(force_modifier);

        MAKE_PTR(PolarityTrackingModifier<2>, polarity_modifier);
        simulator.AddSimulationModifier(polarity_modifier);

        MAKE_PTR_ARGS(MyForce, p_force, (3.0));
        simulator.AddForce(p_force);
        std::cout << "Active force" << endl ; */

        MAKE_PTR(CellAddingModifier<2>, cell_adding_modifier);
        simulator.AddSimulationModifier(cell_adding_modifier);

        // FileFinder file_finder("blah", RelativeTo::ChasteTestOutput);
        // std::cout << file_finder.GetAbsolutePath();

        simulator.SetEndTime(1.0);
        simulator.SetDt(1.0/100.0);
        simulator.SetSamplingTimestepMultiple(1);
        std::cout << "Growing Monolayer" << endl ;
        simulator.Solve();


        // Record mesh
        TrianglesMeshWriter<2,2> triangles_mesh_writer("TestNodeBasedMeshWriter", "node_based_mesh");
        triangles_mesh_writer.WriteFilesUsingMesh(p_mesh);
        



    }
};
