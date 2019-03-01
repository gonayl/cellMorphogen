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
#include "FakePetscSetup.hpp"
#include "CounterSingleton.hpp"

#include "CellDataItemWriter.hpp"

#include "ParabolicGrowingDomainPdeModifier.hpp"
#include "MorphogenCellwiseSourceParabolicPde.hpp"
#include "PolarityTrackingModifier.hpp"

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
#include <complex>      // std::complex, std::polar
#include "CounterSingleton.hpp"


static const double M_TIME_FOR_SIMULATION = 40;
static const double M_NUM_CELLS_ACROSS = 10;
static const double M_UPTAKE_RATE = 10.0 ;
static const double M_DIFFUSION_CONSTANT = 5e-1;
static const double M_DUDT_COEFFICIENT = 1.0;
static const double M_DECAY_COEFFICIENT = 9.0;
static const double M_RADIUS = 100.0;


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
        double count = CounterSingleton::Instance()->GetCount() ;
        double mult = 3.0 + sqrt(count) ;

        for (AbstractCellPopulation<2>::Iterator cell_iter = rCellPopulation.Begin();
                cell_iter != rCellPopulation.End();
                ++cell_iter)
           {
               if (cell_iter->HasCellProperty<CellEndo>())
               {
                   double x = cell_iter->GetCellData()->GetItem("morphogen_grad_x");
                   double y = cell_iter->GetCellData()->GetItem("morphogen_grad_y");
                   force(1) =    mStrength  * mult * (y - ymoy) / (0.55 + std::abs(y - ymoy)) ;
                   force(0) =    mStrength  * mult * (x - xmoy) / (0.55 + std::abs(x - xmoy)) ;

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



class TestMorphogenNode : public AbstractCellBasedTestSuite
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
          PositionDependentCellCycleModel* p_position_model = new PositionDependentCellCycleModel();
          MAKE_PTR(WildTypeCellMutationState, p_state);
          MAKE_PTR(TransitCellProliferativeType, p_transit_type);
          MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
          int age = rand() % 20 + 1;
          unsigned node_index = i;

          double magnitude = 1.0 ;
          double phase = 0.5 ;

          double xdiff = 0.0 ;
          double ydiff = 0.0 ;

          std::complex<double> polar (magnitude,phase);
          int polabs = std::abs(polar) ;
          //std::cout << node_index << endl ;

          if (node_index == 25){
            CellPtr p_cell(new Cell(p_state, p_position_model));
            p_cell->SetCellProliferativeType(p_transit_type);
            // p_cell->SetBirthTime(-20);
            // Initial Condition for Morphogen PDE
            p_cell->GetCellData()->SetItem("morphogen",0.0);
            p_cell->GetCellData()->SetItem("morphogen_grad_x",0.0);
            p_cell->GetCellData()->SetItem("morphogen_grad_y",0.0);

            p_cell->GetCellData()->SetItem("magnitude", magnitude);
            p_cell->GetCellData()->SetItem("phase", phase);
            p_cell->GetCellData()->SetItem("polarity", polabs);
            p_cell->GetCellData()->SetItem("xdiff", xdiff);
            p_cell->GetCellData()->SetItem("ydiff", ydiff);


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

            p_cell->GetCellData()->SetItem("magnitude", magnitude);
            p_cell->GetCellData()->SetItem("phase", phase);
            p_cell->GetCellData()->SetItem("polarity", polabs);
            p_cell->GetCellData()->SetItem("xdiff", xdiff);
            p_cell->GetCellData()->SetItem("ydiff", ydiff);

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

            p_cell->GetCellData()->SetItem("magnitude", 0.0);
            p_cell->GetCellData()->SetItem("phase", 0.0);
            p_cell->GetCellData()->SetItem("polarity", 0.0);
            p_cell->GetCellData()->SetItem("xdiff", xdiff);
            p_cell->GetCellData()->SetItem("ydiff", ydiff);
            p_cell->GetCellData()->SetItem("dist", 0.0);

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
        cell_population.AddCellWriter<CellAgesWriter>();
        MAKE_PTR(CellLabel, p_label);
        MAKE_PTR(CellEndo, p_endo);
        MAKE_PTR(ForceTrackingModifier<2>, force_modifier);

        // Set up cell-based simulation and output directory
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("CEllMorphogen/OSModel/TestMotileForce");

        // Create a force law and pass it to the simulation
        MAKE_PTR(DifferentialAdhesionGeneralisedLinearSpringForce<2>, p_differential_adhesion_force);
        p_differential_adhesion_force->SetMeinekeSpringStiffness(15.0);
        p_differential_adhesion_force->SetHomotypicLabelledSpringConstantMultiplier(4.0);
        p_differential_adhesion_force->SetHeterotypicSpringConstantMultiplier(0.4);
        p_differential_adhesion_force->SetCutOffLength(cut_off_length);
        simulator.AddForce(p_differential_adhesion_force);

        std::cout << "VeGF diffusion" << endl ;

        // Create a parabolic PDE object - see the header file for what the constructor arguments mean
        MAKE_PTR_ARGS(MorphogenCellwiseSourceParabolicPde<2>, p_pde, (cell_population, M_DUDT_COEFFICIENT,M_DIFFUSION_CONSTANT,M_RADIUS));

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
        simulator.AddSimulationModifier(force_modifier);

        std::cout << "Cell Labelling" << endl ;

        // Record mesh
        //TrianglesMeshWriter<2,2> triangles_mesh_writer("TestNodeMeshwriter", "triangles_mesh");
        //TetrahedralMesh<2,2> triangles_mesh = static_cast<NodeBasedCellPopulation<2>(simulator.rGetCellPopulation()).rGetMesh();
        //triangles_mesh_writer.WriteFilesUsingMesh(p_mesh);


        cell_population.AddCellWriter<CellLabelWriter>();

        for (AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
             cell_iter != cell_population.End();
             ++cell_iter)
           {

             unsigned node_index = cell_population.GetNodeCorrespondingToCell(*cell_iter)->GetIndex();
             c_vector<double, 2> location = cell_population.GetLocationOfCellCentre(*cell_iter);

             //cell_iter->GetCellData()->SetItem("xpos", location[0]);
             //cell_iter->GetCellData()->SetItem("ypos", location[1]);
             cell_iter->GetCellData()->SetItem("xposini", location[0]);
             cell_iter->GetCellData()->SetItem("yposini", location[1]);

               if (node_index == 113 )
               {
                 cell_iter->AddCellProperty(p_label);
                 cell_iter->AddCellProperty(p_endo);
               }
               if (node_index == 25 )
               {
                 cell_iter->AddCellProperty(p_label);
               }
           }


        //MAKE_PTR(MorphogenTrackingModifier<2>, morphogen_modifier);
        //simulator.AddSimulationModifier(morphogen_modifier);
        MAKE_PTR(RadiusTrackingModifier<2>, radius_modifier);
        simulator.AddSimulationModifier(radius_modifier);
        MAKE_PTR(PositionTrackingModifier<2>, position_modifier);
        simulator.AddSimulationModifier(position_modifier);
        // MAKE_PTR(PolarityTrackingModifier<2>, polarity_modifier);
        // simulator.AddSimulationModifier(polarity_modifier);
        MAKE_PTR_ARGS(MyForce, p_force, (32.0));
        simulator.AddForce(p_force);
        std::cout << "Active force" << endl ;

        // MAKE_PTR(EndothelialAdhesionGeneralisedLinearSpringForce<2>, p_endothelial_adhesion_force);
        // simulator.AddForce(p_endothelial_adhesion_force);

        simulator.SetEndTime(20.0);
        simulator.SetDt(1.0/100.0);
        simulator.SetSamplingTimestepMultiple(1);
        std::cout << "Growing Monolayer" << endl ;
        simulator.Solve();



    }
};
