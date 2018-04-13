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
#include "FakePetscSetup.hpp"

#include "CellDataItemWriter.hpp"

#include "ParabolicGrowingDomainPdeModifier.hpp"
#include "MorphogenCellwiseSourceParabolicPde.hpp"

#include "AbstractForce.hpp"
#include <numeric>
#include <iostream>
#include <cmath>
#include <stdlib.h>

#include "CellLabel.hpp"
#include "CellEndo.hpp"
#include "CellLabelWriter.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "PerimeterTrackingModifier.hpp"
#include "PerimeterDependentCellCycleModel.hpp"

#include "DifferentiatedCellProliferativeType.hpp"

static const double M_TIME_FOR_SIMULATION = 40;
static const double M_NUM_CELLS_ACROSS = 10;
static const double M_UPTAKE_RATE = 10.0 ;
static const double M_DIFFUSION_CONSTANT = 5e-1;
static const double M_DUDT_COEFFICIENT = 1.0;
static const double M_DECAY_COEFFICIENT = 9.0;
static const double M_RADIUS = 2.0;


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

        std::vector<double> moy_x_morphogen_grad;
        std::vector<double> moy_y_morphogen_grad;

        double morphogen_grad_x = 0.0 ;
        double morphogen_grad_y = 0.0 ;

        for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
                cell_iter != rCellPopulation.End();
                ++cell_iter)
         {

        morphogen_grad_x = cell_iter->GetCellData()->GetItem("morphogen_grad_x");
        morphogen_grad_y = cell_iter->GetCellData()->GetItem("morphogen_grad_y");
        moy_x_morphogen_grad.push_back(morphogen_grad_x) ;
        moy_y_morphogen_grad.push_back(morphogen_grad_y) ;

         }

        // double xmoy = std::accumulate(moy_x_morphogen_grad.begin(), moy_x_morphogen_grad.end(), 0)/moy_x_morphogen_grad.size();
        // double ymoy = std::accumulate(moy_y_morphogen_grad.begin(), moy_y_morphogen_grad.end(), 0)/moy_y_morphogen_grad.size();
        double xmoy = 0.0;
        double ymoy = 0.0;


        for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
                cell_iter != rCellPopulation.End();
                ++cell_iter)
           {
               if (cell_iter->HasCellProperty<CellEndo>())
               {

                   double x = cell_iter->GetCellData()->GetItem("morphogen_grad_x");
                   double y = cell_iter->GetCellData()->GetItem("morphogen_grad_y");
                   force(1) =    mStrength * (y - ymoy) / (0.55 + std::abs(y - ymoy)) ;
                   force(0) =    mStrength * (x - xmoy) / (0.55 + std::abs(x - xmoy)) ;

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
          int age = rand() % 20 + 1;
          unsigned node_index = i;
          //std::cout << node_index << endl ;

          if (node_index == 25 or node_index == 88){
            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            // p_cell->SetBirthTime(-20);
            // Initial Condition for Morphogen PDE
            p_cell->GetCellData()->SetItem("morphogen",0.0);
            p_cell->GetCellData()->SetItem("morphogen_grad_x",0.0);
            p_cell->GetCellData()->SetItem("morphogen_grad_y",0.0);

            // Set Target Area so dont need to use a growth model in vertex simulations
            p_cell->GetCellData()->SetItem("target area", 1.0);
            cells.push_back(p_cell);
          } else if (node_index == 113){
            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_diff_type);
            // Initial Condition for Morphogen PDE
            p_cell->GetCellData()->SetItem("morphogen",0.0);
            p_cell->GetCellData()->SetItem("morphogen_grad_x",0.0);
            p_cell->GetCellData()->SetItem("morphogen_grad_y",0.0);

            // Set Target Area so dont need to use a growth model in vertex simulations
            p_cell->GetCellData()->SetItem("target area", 1.0);
            cells.push_back(p_cell);
          } else {
            CellPtr p_cell(new Cell(p_state, p_cycle_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            p_cell->SetBirthTime(-age);
            // Initial Condition for Morphogen PDE
            p_cell->GetCellData()->SetItem("morphogen",0.0);
            p_cell->GetCellData()->SetItem("morphogen_grad_x",0.0);
            p_cell->GetCellData()->SetItem("morphogen_grad_y",0.0);

            // Set Target Area so dont need to use a growth model in vertex simulations
            p_cell->GetCellData()->SetItem("target area", 1.0);
            cells.push_back(p_cell);
          }
      }
   }

public:
  void TestNodeBasedMonolayer()
  {
        // Create a simple mesh

        //HoneycombMeshGenerator generator(3, 3, 0);
        TrianglesMeshReader<2,2> mesh_reader("testoutput/TestNodeMeshWriter/triangles_mesh");
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
        boost::shared_ptr<CellDataItemWriter<2,2> > p_cell_data_item_writer(new CellDataItemWriter<2,2>("morphogen"));
        cell_population.AddCellWriter(p_cell_data_item_writer);

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestCenterBasedModel/Adhesion");

        // Create a force law and pass it to the simulation
        MAKE_PTR(DifferentialAdhesionGeneralisedLinearSpringForce<2>, p_differential_adhesion_force);
        p_differential_adhesion_force->SetMeinekeSpringStiffness(15.0);
        p_differential_adhesion_force->SetHomotypicLabelledSpringConstantMultiplier(0.5);
        p_differential_adhesion_force->SetHeterotypicSpringConstantMultiplier(1.0);
        p_differential_adhesion_force->SetCutOffLength(cut_off_length);
        p_differential_adhesion_force->SetHomotypicLabelledSpringConstantMultiplier(5.0);
        simulator.AddForce(p_differential_adhesion_force);

        std::cout << "VeGF diffusion" << endl ;

        // Create a parabolic PDE object - see the header file for what the constructor arguments mean
        MAKE_PTR_ARGS(CellwiseSourceParabolicPde<2>, p_pde, (cell_population, M_DUDT_COEFFICIENT,M_DIFFUSION_CONSTANT,M_UPTAKE_RATE,M_DECAY_COEFFICIENT,M_RADIUS));

        // Create a constant boundary conditions object, taking the value zero
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (0.0));

        // Create a ParabolicGrowingDomainPdeModifier, which is a simulation modifier, using the PDE and BC objects, and use 'true' to specify that you want a Dirichlet (rather than Neumann) BC
        MAKE_PTR_ARGS(ParabolicGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, false)); // Change this last argument to 'false' for no-flux BCs (Dirichlet??)

        // Optionally, name the PDE state variable (for visualization purposes)
        p_pde_modifier->SetDependentVariableName("morphogen");

        // Add the modifier to the simulation and run it

        std::cout << "Adding Modifier" << endl ;

        p_pde_modifier->SetOutputGradient(true);
        p_pde_modifier->GetOutputGradient();

        simulator.AddSimulationModifier(p_pde_modifier);

        std::cout << "Cell Labelling" << endl ;

        // Record mesh
        //TrianglesMeshWriter<2,2> triangles_mesh_writer("TestNodeMeshwriter", "triangles_mesh");
        //TetrahedralMesh<2,2> triangles_mesh = static_cast<NodeBasedCellPopulation<2>(simulator.rGetCellPopulation()).rGetMesh();
        //triangles_mesh_writer.WriteFilesUsingMesh(p_mesh);


        cell_population.AddCellWriter<CellLabelWriter>();
        MAKE_PTR(CellLabel, p_label);
        MAKE_PTR(CellEndo, p_endo);

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
           {
             unsigned node_index = cell_population.GetNodeCorrespondingToCell(*cell_iter)->GetIndex();
               if (node_index == 113 )
               {
                 cell_iter->AddCellProperty(p_label);
                 cell_iter->AddCellProperty(p_endo);
               }
               if (node_index == 25 or node_index == 88)
               {
                 cell_iter->AddCellProperty(p_label);
               }
           }

        MAKE_PTR_ARGS(MyForce, p_force, (10.0));
        simulator.AddForce(p_force);
        std::cout << "Active force" << endl ;

        simulator.SetEndTime(20.0);
        simulator.SetDt(1.0/100.0);
        simulator.SetSamplingTimestepMultiple(1);
        std::cout << "Growing Monolayer" << endl ;
        simulator.Solve();



    }
};
