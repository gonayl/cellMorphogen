#include "NewEndoGeneratorModifier.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CellEndo.hpp"
#include "CellStalk.hpp"
#include "CellTip.hpp"
#include "CellEpi.hpp"
#include "CellLumen.hpp"
#include "CellPolar.hpp"
#include "CellPeriph.hpp"
#include "CellCore.hpp"
#include "CellMotile.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "PerimeterDependentCellCycleModel.hpp"

#include "CounterSingletonEndo.hpp"

#include <stdlib.h>
#include <numeric>
#include <iostream>
#include <cmath>
using namespace std ;


template<unsigned DIM>
NewEndoGeneratorModifier<DIM>::NewEndoGeneratorModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
NewEndoGeneratorModifier<DIM>::~NewEndoGeneratorModifier()
{
}

template<unsigned DIM>
void NewEndoGeneratorModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void NewEndoGeneratorModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void NewEndoGeneratorModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();
    /**
     * This hack is needed because in the case of a MeshBasedCellPopulation in which
     * multiple cell divisions have occurred over one time step, the Voronoi tessellation
     * (while existing) is out-of-date. Thus, if we did not regenerate the Voronoi
     * tessellation here, an assertion may trip as we try to access a Voronoi element
     * whose index exceeds the number of elements in the out-of-date tessellation.
     *
     * \todo work out how to properly fix this (#1986)
     */
    if (bool(dynamic_cast<MeshBasedCellPopulation<DIM>*>(&rCellPopulation)))
    {
        static_cast<MeshBasedCellPopulation<DIM>*>(&(rCellPopulation))->CreateVoronoiTessellation();
    }

    if (dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) == NULL)
    {
        EXCEPTION("Border Tracking Modifier is to be used with a VertexBasedCellPopulation only");
    }

    //cout << endo_count_global << "endo stalk fixed cell(s) "  << endl ;

    vector<double> endo_stalk_vector ;
    endo_stalk_vector.push_back(162);
    endo_stalk_vector.push_back(148);
    endo_stalk_vector.push_back(471);
    endo_stalk_vector.push_back(59);
    endo_stalk_vector.push_back(13);
    endo_stalk_vector.push_back(1200);

/*
    for (typename AbstractCellPopulation<DIM,DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)

    {
      if ( cell_iter->template HasCellProperty<CellPeriph>() && cell_iter->template HasCellProperty<CellEpi>())
      {

        int cell_endo_location = rCellPopulation.GetLocationIndexUsingCell(*cell_iter) ;
        endo_stalk_vector.push_back(cell_endo_location) ;


      }

    }
*/
    // MAKE_PTR(CellBoundary, p_border);

    //double time_max = 40.0 ;
    //double chance_2_endo ;

    //double count = CounterSingletonEndo::Instance()->GetCount() ;
    double time_elapsed = SimulationTime::Instance()->GetTime() ;
    int need_tip = 0 ;
    int need_stalk = 0 ;
    double count = CounterSingletonEndo::Instance()->GetCount() ;
    double tip_index ;
    double stalk_index ;

    double modulo = fmod(time_elapsed / 12.0 , 1.0) ;
    //std::cout << modulo << " / " << count << std::endl ;

    if (modulo == 0 && time_elapsed > 0 && count == 0)
    {
      CounterSingletonEndo::Instance()->IncrementCounter();
      std::cout << "now!" << std::endl ;
    }

    //std::cout << time_elapsed << std::endl ;
    //std::cout << modulo << std::endl ;

    int randoEndo = rand() % endo_stalk_vector.size() ;
    //chance_2_endo = time_elapsed/(time_max*count) ;

    //if (chance_2_endo >= 1)
    //{cout << randoEndo << endl ;}


    // Iterate over cell population
    //VertexBasedCellPopulation<DIM>* p_cell_population = dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) ;

    // labelling new stalk cell
    for (typename AbstractCellPopulation<DIM, DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {

      //VertexElement<DIM,DIM>* p_element = p_cell_population->GetElementCorrespondingToCell(*cell_iter);

      //RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

      //cout << chance_2_endo << endl ;
      int n_neighbour_endo = 0;
      // bool need_tip = 0 ;

      int cell_location = rCellPopulation.GetLocationIndexUsingCell(*cell_iter) ;

      set<unsigned> neighbours_of_cell = rCellPopulation.GetNeighbouringLocationIndices(*cell_iter) ;
      for (set<unsigned>::iterator neighbour_iter = neighbours_of_cell.begin();
           neighbour_iter != neighbours_of_cell.end();
           ++neighbour_iter)
      {
          // Determine whether this neighbour is also endo
          CellPtr p_neighbour_cell = rCellPopulation.GetCellUsingLocationIndex(*neighbour_iter);
          bool neighbour_is_endo = p_neighbour_cell->template HasCellProperty<CellEndo>();

          if ( neighbour_is_endo)
          {
              n_neighbour_endo++;
          }
      }




      if (cell_iter->template HasCellProperty<CellPeriph>() && cell_iter->template HasCellProperty<CellEpi>() && count == 1 && time_elapsed > 0 && n_neighbour_endo == 0 && cell_location == endo_stalk_vector[randoEndo])
      {

        cout << "cell " << rCellPopulation.GetLocationIndexUsingCell(*cell_iter) << " needs to be changed into a stalk cell" << endl ;

        stalk_index =  rCellPopulation.GetLocationIndexUsingCell(*cell_iter) ;
        need_stalk = 1 ;

        // labelling tip cell
          int is_tip = 0 ;

          set<unsigned> neighbours_of_cell = rCellPopulation.GetNeighbouringLocationIndices(*cell_iter) ;
          cout << "itering over all cells around cell " << rCellPopulation.GetLocationIndexUsingCell(*cell_iter) << endl ;


          for (set<unsigned>::iterator neighbour_iter = neighbours_of_cell.begin();
               neighbour_iter != neighbours_of_cell.end();
               ++neighbour_iter)
          {

              // Determine whether any neighbour cell is already tip
              CellPtr p_neighbour_cell = rCellPopulation.GetCellUsingLocationIndex(*neighbour_iter);
              bool neighbour_is_tip = p_neighbour_cell->template HasCellProperty<CellTip>();

              if ( neighbour_is_tip)
              {
                  is_tip++;

              }


          }

          for (set<unsigned>::iterator neighbour_iter = neighbours_of_cell.begin();
               neighbour_iter != neighbours_of_cell.end();
               ++neighbour_iter)
          {
              // Determine whether this neighbour is epi and non periph
              CellPtr p_neighbour_cell = rCellPopulation.GetCellUsingLocationIndex(*neighbour_iter);
              bool neighbour_is_epi = p_neighbour_cell->template HasCellProperty<CellEpi>();
              bool neighbour_is_core = p_neighbour_cell->template HasCellProperty<CellCore>();
              bool neighbour_is_endo = p_neighbour_cell->template HasCellProperty<CellEndo>();

              if (neighbour_is_epi && neighbour_is_core && is_tip == 0 && neighbour_is_endo == 0 && need_stalk == 1)
              {

                tip_index = rCellPopulation.GetLocationIndexUsingCell(p_neighbour_cell) ;
                is_tip++ ;
                need_tip = 1 ;
              }

            }


      // comment here

      }


    }

    if (need_stalk && need_tip == 1)
    {


      CellPtr cell_stalk = rCellPopulation.GetCellUsingLocationIndex(stalk_index);
      PerimeterDependentCellCycleModel* p_elong_model = new PerimeterDependentCellCycleModel();
      cell_stalk->template RemoveCellProperty<CellEpi>();
      cell_stalk->AddCellProperty(CellPropertyRegistry::Instance()->Get<CellEndo>());
      cell_stalk->AddCellProperty(CellPropertyRegistry::Instance()->Get<CellStalk>());
      cell_stalk->SetCellCycleModel(p_elong_model) ;
      cell_stalk->InitialiseCellCycleModel() ;

      cell_stalk->GetCellData()->SetItem("tagVessel",-1);


      //cell_iter->GetCellData()->SetItem("tagVessel",cell_iter->GetCellData()->GetItem("cellIndex"));
      //cell_iter->SetCellProliferativeType(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());


      //CounterSingletonEndo::Instance()->IncrementCounter();
/*
      // Adding the stalk cell to the obstruction file

      vector<double> boundary_nodes_pos ;
      VertexBasedCellPopulation<DIM>* p_cell_population = dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) ;
      VertexElement<DIM,DIM>* p_element ;

      p_element = p_cell_population->GetElementCorrespondingToCell(*cell_iter);

      unsigned num_nodes_in_element = p_element->GetNumNodes();
      for (unsigned local_index=0; local_index<num_nodes_in_element; local_index++)
      {
          unsigned node_index = p_element->GetNodeGlobalIndex(local_index);
          unsigned n_elements = rCellPopulation.GetNode(node_index)->GetNumContainingElements() ;

          if (rCellPopulation.GetNode(node_index)->IsBoundaryNode() && n_elements > 1  )
          {
            boundary_nodes_pos.push_back(rCellPopulation.GetNode(node_index)->rGetLocation()[0]);
            boundary_nodes_pos.push_back(rCellPopulation.GetNode(node_index)->rGetLocation()[1]);
          }
      }

      assert(boundary_nodes_pos.size() == 4) ;
      std::cout << "(" << boundary_nodes_pos[0] << " ; " << boundary_nodes_pos[1] << ")(" << boundary_nodes_pos[2] << " ; " << boundary_nodes_pos[3] << ")" << endl ;

      double a; double b; double l; double r;
      c_vector<double, 2> direct_vector ;

      c_vector<double, 2> centre_of_cell = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

      direct_vector[0] = cell_iter->GetCellData()->GetItem("morphogen_grad_x") ;
      direct_vector[1] = cell_iter->GetCellData()->GetItem("morphogen_grad_y") ;

      a =  direct_vector[1] ; // a
      b =  - direct_vector[0] ; // b
      l = - direct_vector[1]*boundary_nodes_pos[0] + direct_vector[0]*boundary_nodes_pos[1] ; // l
      r = - direct_vector[1]*boundary_nodes_pos[2] + direct_vector[0]*boundary_nodes_pos[3] ; // r
      // WARNING : the order is important !!

      cout << "pushing data for cell " << rCellPopulation.GetLocationIndexUsingCell(*cell_iter) << endl ;
      cell_iter->GetCellData()->SetItem("xpos_l", boundary_nodes_pos[0]) ;
      cell_iter->GetCellData()->SetItem("ypos_l", boundary_nodes_pos[1]) ;
      cell_iter->GetCellData()->SetItem("xpos_r", boundary_nodes_pos[2]) ;
      cell_iter->GetCellData()->SetItem("ypos_r", boundary_nodes_pos[3]) ;

      ofstream myfile ("example.txt", std::ios::app);

      if (myfile.is_open())
      {
        myfile << a << " " << b << " " << l << " " << r << " " << centre_of_cell[0] << " " << centre_of_cell[1] << " " ;
      }

      myfile.close();
      */

      cout << "cell " << tip_index << " needs to be changed into a tip cell" << endl ;
      CellPtr cell_tip = rCellPopulation.GetCellUsingLocationIndex(tip_index);
      cell_tip->template RemoveCellProperty<CellEpi>();
      cell_tip->AddCellProperty(CellPropertyRegistry::Instance()->Get<CellEndo>());
      cell_tip->AddCellProperty(CellPropertyRegistry::Instance()->Get<CellTip>());
      cell_tip->AddCellProperty(CellPropertyRegistry::Instance()->Get<CellMotile>());
      cell_tip->GetCellData()->SetItem("tagVessel",cell_tip->GetCellData()->GetItem("cellIndex"));

      need_tip = 0 ;
      CounterSingletonEndo::Instance()->ResetCountToZero();
      need_stalk = 0 ;

    }
    else if (need_stalk == 1 && need_tip == 0)
    {
      std::cout << "but there is no available tip cell, keep searching" << std::endl ;
    }

}

template<unsigned DIM>
void NewEndoGeneratorModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class NewEndoGeneratorModifier<1>;
template class NewEndoGeneratorModifier<2>;
template class NewEndoGeneratorModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NewEndoGeneratorModifier)
