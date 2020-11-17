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
#include "CellBase.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "SimulationParameters.hpp"

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
    vector<double> endo_stalk_vector ;
    vector<double> endo_base_vector ;
    double base_min_dist = 9999 ;
    int randoEndo ;

    double modulo = fmod(time_elapsed / 8.0 , 1.0) ;
    //std::cout << modulo << " / " << count << std::endl ;

    if (modulo == 0 && time_elapsed > 0 && count == 0)
    {
      CounterSingletonEndo::Instance()->IncrementCounter();
      std::cout << "now!" << std::endl ;

    }


      // iterate over cell population to create vector with potential to be transofrmed into stalk cell as well as the base cells vector

      for (typename AbstractCellPopulation<DIM,DIM>::Iterator cell_iter = rCellPopulation.Begin();
           cell_iter != rCellPopulation.End();
           ++cell_iter)
      {
        if (cell_iter->template HasCellProperty<CellPeriph>() && cell_iter->template HasCellProperty<CellEpi>())
        {
          int cell_endo_location = rCellPopulation.GetLocationIndexUsingCell(*cell_iter) ;
          endo_stalk_vector.push_back(cell_endo_location) ;
        }
        if (cell_iter->template HasCellProperty<CellBase>())
        {
          int cell_base_location = rCellPopulation.GetLocationIndexUsingCell(*cell_iter) ;
          endo_base_vector.push_back(cell_base_location) ;
        }

      }

      // iterate over base cells to compute min distance

      for (vector<double>::iterator base_iter = endo_base_vector.begin();
           base_iter != endo_base_vector.end();
           ++base_iter)
      {
          CellPtr p_current_base_cell = rCellPopulation.GetCellUsingLocationIndex(*base_iter);
          c_vector<double, DIM> current_base_cell_location = rCellPopulation.GetLocationOfCellCentre(p_current_base_cell) ;

          for (vector<double>::iterator base_iter2 = endo_base_vector.begin();
               base_iter2 != endo_base_vector.end();
               ++base_iter2)
          {
              CellPtr p_next_base_cell = rCellPopulation.GetCellUsingLocationIndex(*base_iter2);
              c_vector<double, DIM> next_base_cell_location = rCellPopulation.GetLocationOfCellCentre(p_next_base_cell) ;
              double current_dist = sqrt((current_base_cell_location[0] - next_base_cell_location[0])*(current_base_cell_location[0] - next_base_cell_location[0]) + (current_base_cell_location[1] - next_base_cell_location[1])*(current_base_cell_location[1] - next_base_cell_location[1])) ;

              if ( current_dist < base_min_dist && current_dist > 0)
              {
                  base_min_dist = current_dist ;
              }
          }
      }
      base_min_dist = base_min_dist ;
      randoEndo = rand() % endo_stalk_vector.size() ;


    //std::cout << time_elapsed << std::endl ;
    //std::cout << modulo << std::endl ;

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
      double min_edg = cell_iter->GetCellData()->GetItem("min") ;
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

      c_vector<double, DIM> current_epi_cell_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter) ;
      double min_dist = 0;
      double previous_dist = 9999 ;
      for (vector<double>::iterator base_iter3 = endo_base_vector.begin();
           base_iter3 != endo_base_vector.end();
           ++base_iter3)
      {
          CellPtr p_base_cell = rCellPopulation.GetCellUsingLocationIndex(*base_iter3);
          c_vector<double, DIM> base_cell_location = rCellPopulation.GetLocationOfCellCentre(p_base_cell) ;
          double current_dist = sqrt((current_epi_cell_location[0] - base_cell_location[0])*(current_epi_cell_location[0] - base_cell_location[0]) + (current_epi_cell_location[1] - base_cell_location[1])*(current_epi_cell_location[1] - base_cell_location[1])) ;

          if (current_dist < previous_dist)
          {
              previous_dist = current_dist ;
          }
      }

      min_dist = previous_dist ;

      if (cell_iter->template HasCellProperty<CellPeriph>() && cell_iter->template HasCellProperty<CellEpi>() && count == 1 && time_elapsed > 0 && n_neighbour_endo == 0 && cell_location == endo_stalk_vector[randoEndo] && min_dist > base_min_dist && min_edg > 0.25)
      {

        cout << "cell " << rCellPopulation.GetLocationIndexUsingCell(*cell_iter) << " needs to be changed into a stalk cell" << endl ;
        stalk_index =  rCellPopulation.GetLocationIndexUsingCell(*cell_iter) ;
        need_stalk = 1 ;
        cout << min_dist << endl ;

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

      //FOR LUMEN
      cell_stalk->GetCellData()->SetItem("cellIndex",SimulationParameters::getNextIndex());
      cell_stalk->GetCellData()->SetItem("timeFromLastLumenGeneration",0);
      cell_stalk->GetCellData()->SetItem("hadALumenDivision",0);
      cell_stalk->GetCellData()->SetItem("lumenNearby",1);
      cell_stalk->GetCellData()->SetItem("vecPolaX",0);
      cell_stalk->GetCellData()->SetItem("vecPolaY",0);

      cell_stalk->GetCellData()->SetItem("tagVessel",-1);


      //cell_iter->GetCellData()->SetItem("tagVessel",cell_iter->GetCellData()->GetItem("cellIndex"));
      //cell_iter->SetCellProliferativeType(CellPropertyRegistry::Instance()->Get<DifferentiatedCellProliferativeType>());


      //CounterSingletonEndo::Instance()->IncrementCounter();

      cout << "cell " << tip_index << " needs to be changed into a tip cell" << endl ;
      CellPtr cell_tip = rCellPopulation.GetCellUsingLocationIndex(tip_index);
      cell_tip->template RemoveCellProperty<CellEpi>();
      cell_tip->AddCellProperty(CellPropertyRegistry::Instance()->Get<CellEndo>());
      cell_tip->AddCellProperty(CellPropertyRegistry::Instance()->Get<CellTip>());
      cell_tip->AddCellProperty(CellPropertyRegistry::Instance()->Get<CellMotile>());
      cell_tip->GetCellData()->SetItem("tagVessel",cell_tip->GetCellData()->GetItem("cellIndex"));
      cell_tip->GetCellData()->SetItem("cellIndex",SimulationParameters::getNextIndex());
      cell_tip->GetCellData()->SetItem("timeFromLastLumenGeneration",0);
      cell_tip->GetCellData()->SetItem("hadALumenDivision",0);
      cell_tip->GetCellData()->SetItem("lumenNearby",1);
      cell_tip->GetCellData()->SetItem("vecPolaX",0);
      cell_tip->GetCellData()->SetItem("vecPolaY",0);


            // Adding the stalk cell to the obstruction file

            vector<double> boundary_nodes_pos ;
            VertexBasedCellPopulation<DIM>* p_cell_population = dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) ;
            VertexElement<DIM,DIM>* p_element ;

            p_element = p_cell_population->GetElementCorrespondingToCell(cell_stalk);

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

            c_vector<double, 2> centre_of_cell = rCellPopulation.GetLocationOfCellCentre(cell_stalk);

            direct_vector[0] = cell_stalk->GetCellData()->GetItem("morphogen_grad_x") ;
            direct_vector[1] = cell_stalk->GetCellData()->GetItem("morphogen_grad_y") ;

            a =  direct_vector[1] ; // a
            b =  - direct_vector[0] ; // b
            l = - direct_vector[1]*boundary_nodes_pos[0] + direct_vector[0]*boundary_nodes_pos[1] ; // l
            r = - direct_vector[1]*boundary_nodes_pos[2] + direct_vector[0]*boundary_nodes_pos[3] ; // r
            // WARNING : the order is important !!

            cout << "pushing data for cell " << rCellPopulation.GetLocationIndexUsingCell(cell_stalk) << endl ;
            cell_stalk->GetCellData()->SetItem("xpos_l", boundary_nodes_pos[0]) ;
            cell_stalk->GetCellData()->SetItem("ypos_l", boundary_nodes_pos[1]) ;
            cell_stalk->GetCellData()->SetItem("xpos_r", boundary_nodes_pos[2]) ;
            cell_stalk->GetCellData()->SetItem("ypos_r", boundary_nodes_pos[3]) ;

            ofstream myfile ("example.txt", std::ios::app);

            if (myfile.is_open())
            {
              myfile << a << " " << b << " " << l << " " << r << " " << centre_of_cell[0] << " " << centre_of_cell[1] << " " ;
            }

            myfile.close();


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
