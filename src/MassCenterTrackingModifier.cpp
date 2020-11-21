#include "MassCenterTrackingModifier.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CellEndo.hpp"
#include "CellEpi.hpp"
#include "CellMotile.hpp"
#include "CellLumen.hpp"
#include "CellPolar.hpp"
#include "TransitCellProliferativeType.hpp"
#include "Debug.hpp"
#include "CellVessel.hpp"
#include "CellPeriph.hpp"
#include "SimulationParameters.hpp"
#include <stdlib.h>
#include "CellTip.hpp"
#include "CellStalk.hpp"
using namespace std ;

#include "SimulationParameters.hpp"

double xmoy ;
double ymoy ;

template<unsigned DIM>
MassCenterTrackingModifier<DIM>::MassCenterTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
MassCenterTrackingModifier<DIM>::~MassCenterTrackingModifier()
{
}

template<unsigned DIM>
void MassCenterTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void MassCenterTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void MassCenterTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();
    //template <class CellStalk> class CellStalk;
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
        EXCEPTION("Mass Center Tracking Modifier is to be used with a VertexBasedCellPopulation only");
    }

    // Iterate over cell population

    std::vector<double> moy_x_pos;
    std::vector<double> moy_y_pos;


    for (typename AbstractCellPopulation<DIM, DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
      c_vector<double, DIM> cell_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter) ;
      moy_x_pos.push_back(cell_location[0]);
      moy_y_pos.push_back(cell_location[1]);

    }

    double xsum = std::accumulate(moy_x_pos.begin(), moy_x_pos.end(), 0.0);
    double ysum = std::accumulate(moy_y_pos.begin(), moy_y_pos.end(), 0.0);
    xmoy = xsum / moy_x_pos.size() ;
    ymoy = ysum / moy_y_pos.size() ;

    // std::cout << "x moy = " << xmoy << " " << "y moy = " << ymoy << std::endl;


    // Iterate over cells and calculate the distance between each cell and the mass center


    for (typename AbstractCellPopulation<DIM, DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
      cell_iter->GetCellData()->SetItem("mass_center_x", xmoy);
      cell_iter->GetCellData()->SetItem("mass_center_y", ymoy);
    }


    //Regarde le nomobre de tip voisin

    for (typename AbstractCellPopulation<DIM, DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
      CellPtr pCell = *cell_iter;

      //regarde si il y a des voisines des Tip
      if (pCell->HasCellProperty<CellStalk>())
      {

        pCell->GetCellData()->SetItem("have_tip_neighboor", 0);
        pCell->GetCellData()->SetItem("have_vessel_neighboor", 0);

        VertexBasedCellPopulation<DIM>* p_cell_population = dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) ;
        std::set<unsigned> neighbour_indices = p_cell_population->GetNeighbouringLocationIndices(*cell_iter);


        for (std::set<unsigned>::iterator neighbour_iter = neighbour_indices.begin();
             neighbour_iter != neighbour_indices.end();
             ++neighbour_iter)
        {
          CellPtr p_neighbour_cell = p_cell_population->GetCellUsingLocationIndex(*neighbour_iter);

          bool neighbour_is_Tip = p_neighbour_cell->template HasCellProperty<CellTip>();
          bool neighbour_is_Vessel = p_neighbour_cell->template HasCellProperty<CellVessel>();

          if ( neighbour_is_Tip == 1)
          {
            //cout << "voisine tip cell" << endl ;
            double nbrVoisin_tip = pCell->GetCellData()->GetItem("have_tip_neighboor");
            pCell->GetCellData()->SetItem("have_tip_neighboor", nbrVoisin_tip+1);
          }
          if ( neighbour_is_Vessel == 1)
          {
            //cout << "voisine vessel cell" << endl ;
            double nbrVoisin_vessel = pCell->GetCellData()->GetItem("have_vessel_neighboor");
            pCell->GetCellData()->SetItem("have_vessel_neighboor", nbrVoisin_vessel+1);
          }
        }
      }


    //Retire les tips en réalisant des vaisseaux
    if(SimulationParameters::GENERATE_VESSEL){

      for (typename AbstractCellPopulation<DIM, DIM>::Iterator cell_iter = rCellPopulation.Begin();
           cell_iter != rCellPopulation.End();
           ++cell_iter)
      {
        CellPtr pCell = *cell_iter;

        //attribue la cellule stalk à un vaisseau
        if (pCell->HasCellProperty<CellStalk>())
        {
          //seulement si elle ne faisait pas partie d'un vaisseau
          if(pCell->GetCellData()->GetItem("tagVessel")==-1)
          {
            VertexBasedCellPopulation<DIM>* p_cell_population = dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) ;
            std::set<unsigned> neighbour_indices = p_cell_population->GetNeighbouringLocationIndices(*cell_iter);



            for (std::set<unsigned>::iterator neighbour_iter = neighbour_indices.begin();
                 neighbour_iter != neighbour_indices.end();
                 ++neighbour_iter)
            {
              CellPtr p_neighbour_cell = p_cell_population->GetCellUsingLocationIndex(*neighbour_iter);

              if(p_neighbour_cell->HasCellProperty<CellTip>())
              {
                pCell->GetCellData()->SetItem("tagVessel",p_neighbour_cell->GetCellData()->GetItem("tagVessel"));
              }

            }
        }
      }

      //rtransforme une tip en non tip si elle rencontre un vaisseau autre que le siens ou la périphérie
      if (pCell->HasCellProperty<CellTip>())
      {
        bool removeTip = false;
        VertexBasedCellPopulation<DIM>* p_cell_population = dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) ;
        std::set<unsigned> neighbour_indices = p_cell_population->GetNeighbouringLocationIndices(*cell_iter);

        double n_boundary_nodes = 0 ;
        VertexElement<DIM,DIM>* p_element = p_cell_population->GetElementCorrespondingToCell(pCell);
        unsigned num_nodes_in_element = p_element->GetNumNodes();
        for (unsigned local_index=0; local_index<num_nodes_in_element; local_index++)
        {
            unsigned node_index = p_element->GetNodeGlobalIndex(local_index);

            if (rCellPopulation.GetNode(node_index)->IsBoundaryNode() )
            {

                n_boundary_nodes++ ;

            }
        }

        for (std::set<unsigned>::iterator neighbour_iter = neighbour_indices.begin();
             neighbour_iter != neighbour_indices.end();
             ++neighbour_iter)
        {
          CellPtr p_neighbour_cell = p_cell_population->GetCellUsingLocationIndex(*neighbour_iter);

          if(p_neighbour_cell->HasCellProperty<CellEndo>())
          {
            if(p_neighbour_cell->GetCellData()->GetItem("tagVessel")>-1 && p_neighbour_cell->GetCellData()->GetItem("tagVessel") !=pCell->GetCellData()->GetItem("tagVessel"))
            {
              removeTip = true;
              MARK ;
            }

          }
        }

        if ((SimulationParameters::VEGF_KO == false) && n_boundary_nodes > 0)
        {
            removeTip = true;
            MARK ;
        }


        //on la met stalk si en contact avec un autre vaisseau
        if(removeTip && (SimulationParameters::VEGF_KO == false))
        {
          pCell->RemoveCellProperty<CellTip>();
          pCell->RemoveCellProperty<CellMotile>();
          boost::shared_ptr<AbstractCellProperty> p_stalk(CellPropertyRegistry::Instance()->Get<CellStalk>());
          boost::shared_ptr<AbstractCellProperty> p_vessel(CellPropertyRegistry::Instance()->Get<CellVessel>());

          pCell->AddCellProperty(p_stalk);
          pCell->AddCellProperty(p_vessel);
          if (SimulationParameters::SMALL_SIMU)
          {
            pCell->GetCellData()->SetItem("target area", 0.8);
          }
          else if(SimulationParameters::SMALL_SIMU == false)
          {
            pCell->GetCellData()->SetItem("target area", 0.3);
          }

          boost::shared_ptr<AbstractCellProperty> p_transit_type =
              pCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
          pCell->SetCellProliferativeType(p_transit_type);

          for (typename AbstractCellPopulation<DIM, DIM>::Iterator endo_iter = rCellPopulation.Begin();
               endo_iter != rCellPopulation.End();
               ++endo_iter)
          {
            CellPtr pEndoCell = *endo_iter ;
            if (pEndoCell->HasCellProperty<CellStalk>() && pEndoCell->GetCellData()->GetItem("tagVessel") == pCell->GetCellData()->GetItem("tagVessel"))
            {
              pEndoCell->AddCellProperty(p_vessel);
              if (SimulationParameters::SMALL_SIMU)
              {
                pEndoCell->GetCellData()->SetItem("target area", 0.8);
              }
              else if(SimulationParameters::SMALL_SIMU == false)
              {
                pEndoCell->GetCellData()->SetItem("target area", 0.3);
              }
            }
          }

        }
      }


      }


    }
    //Supprime les tip cote a cote
    else{

        //change les Tip qui se touchent en stalk (les voisines seulement)
        if (pCell->HasCellProperty<CellTip>() && (SimulationParameters::VEGF_KO == false))
        {

          VertexBasedCellPopulation<DIM>* p_cell_population = dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) ;
          std::set<unsigned> neighbour_indices = p_cell_population->GetNeighbouringLocationIndices(*cell_iter);


          for (std::set<unsigned>::iterator neighbour_iter = neighbour_indices.begin();
               neighbour_iter != neighbour_indices.end();
               ++neighbour_iter)
          {
            CellPtr p_neighbour_cell = p_cell_population->GetCellUsingLocationIndex(*neighbour_iter);

            //Cell Epithéliale
            bool neighbour_is_Tip = p_neighbour_cell->template HasCellProperty<CellTip>();
            if ( neighbour_is_Tip == 1)
            {
              p_neighbour_cell->RemoveCellProperty<CellTip>();
              p_neighbour_cell->RemoveCellProperty<CellMotile>();
              boost::shared_ptr<AbstractCellProperty> p_stalk(CellPropertyRegistry::Instance()->Get<CellStalk>());
              boost::shared_ptr<AbstractCellProperty> p_vessel(CellPropertyRegistry::Instance()->Get<CellVessel>());
              p_neighbour_cell->AddCellProperty(p_stalk);
              p_neighbour_cell->AddCellProperty(p_vessel);
              if (SimulationParameters::SMALL_SIMU)
              {
                p_neighbour_cell->GetCellData()->SetItem("target area", 0.8);
              }
              else if(SimulationParameters::SMALL_SIMU == false)
              {
                p_neighbour_cell->GetCellData()->SetItem("target area", 0.3);
              }
              boost::shared_ptr<AbstractCellProperty> p_transit_type =
                  pCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
              pCell->SetCellProliferativeType(p_transit_type);

              for (typename AbstractCellPopulation<DIM, DIM>::Iterator endo_iter = rCellPopulation.Begin();
                   endo_iter != rCellPopulation.End();
                   ++endo_iter)
              {
                CellPtr pEndoCell = *endo_iter ;
                if (pEndoCell->HasCellProperty<CellStalk>() && pEndoCell->GetCellData()->GetItem("tagVessel") == p_neighbour_cell->GetCellData()->GetItem("tagVessel"))
                {
                  pEndoCell->AddCellProperty(p_vessel);
                  if (SimulationParameters::SMALL_SIMU)
                  {
                    pEndoCell->GetCellData()->SetItem("target area", 0.8);
                  }
                  else if(SimulationParameters::SMALL_SIMU == false)
                  {
                    pEndoCell->GetCellData()->SetItem("target area", 0.3);
                  }
                }
              }
            }
          }
        }


      }
  }

}

template<unsigned DIM>
double MassCenterTrackingModifier<DIM>::GetMassCenter()
{
    return xmoy ;
}

template<unsigned DIM>
void MassCenterTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class MassCenterTrackingModifier<1>;
template class MassCenterTrackingModifier<2>;
template class MassCenterTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(MassCenterTrackingModifier)
