#include "LumenGenerationModifier.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include <stdlib.h>
#include <math.h>
#include "CellEpi.hpp"
#include "CellEndo.hpp"
#include "CellLumen.hpp"
#include "SimulationParameters.hpp"
#include "SimulationParameters.hpp"


template<unsigned DIM>
LumenGenerationModifier<DIM>::LumenGenerationModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
LumenGenerationModifier<DIM>::~LumenGenerationModifier()
{
}

template<unsigned DIM>
void LumenGenerationModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void LumenGenerationModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void LumenGenerationModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
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
        EXCEPTION("Direction Modifier is to be used with a VertexBasedCellPopulation only");
    }

    // MAKE_PTR(CellBoundary, p_border);

    //double num_cells = 0 ;
    // Iterate over cell population
    //VertexBasedCellPopulation<DIM>* p_cell_population = dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) ;
    for (typename AbstractCellPopulation<DIM, DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {



      CellPtr pCell = *cell_iter;
      //IF Epithéliale
      if (pCell->HasCellProperty<CellEpi>())
      {
        double timeFromLastGen = pCell->GetCellData()->GetItem("timeFromLastLumenGeneration");
        pCell->GetCellData()->SetItem("lumenNearby",0);
        double xl = pCell->GetCellData()->GetItem("vecPolaX");
        double yl = pCell->GetCellData()->GetItem("vecPolaY");

        double norme = sqrt(xl*xl + yl*yl);
        //pCell->GetCellData()->SetItem("timeFromLastLumenGeneration", 0);
      if(norme >SimulationParameters::THRESHOLD_POLARISATION_EPI)
      {
        pCell->GetCellData()->SetItem("timeFromLastLumenGeneration", timeFromLastGen + SimulationParameters::TIMESTEP);
      }
      else
      {
        //Si le vecteur est en dessous du trheshold on le mets à 0
        pCell->GetCellData()->SetItem("timeFromLastLumenGeneration", 0);
      }


        //on cherche la cellule 2
        VertexBasedCellPopulation<DIM>* p_cell_population = dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) ;
        std::set<unsigned> neighbour_indices = p_cell_population->GetNeighbouringLocationIndices(*cell_iter);

        // If this cell has any neighbours (as defined by mesh/population/interaction distance)...
        if (!neighbour_indices.empty() )
        {

            for (std::set<unsigned>::iterator neighbour_iter = neighbour_indices.begin();
                 neighbour_iter != neighbour_indices.end();
                 ++neighbour_iter)
            {
              CellPtr p_neighbour_cell = p_cell_population->GetCellUsingLocationIndex(*neighbour_iter);

              //Cell Epithéliale
              bool neighbour_is_epiDir = p_neighbour_cell->template HasCellProperty<CellEpi>();

              if ( neighbour_is_epiDir == 1)
              {
                //on regarde si c'est deux cellules qui viennent de se diviser
                if(pCell->GetCellData()->GetItem("cellIndex")==p_neighbour_cell->GetCellData()->GetItem("cellIndex"))
                {

                  //on change l'index de la voisine
                  p_neighbour_cell->GetCellData()->SetItem("cellIndex",SimulationParameters::getNextIndex());

                  //on regarde si c'est une division pour créer un lumen
                  if(pCell->GetCellData()->GetItem("hadALumenDivision") == 1)
                  {


                    //on retire le tag de la division
                    pCell->GetCellData()->SetItem("hadALumenDivision",0);
                    p_neighbour_cell->GetCellData()->SetItem("hadALumenDivision",0);
                    pCell->GetCellData()->SetItem("timeFromLastLumenGeneration", 0);
                    p_neighbour_cell->GetCellData()->SetItem("timeFromLastLumenGeneration", 0);

                    //on regarde si elle est pointée par le vecteur

                    c_vector<double, DIM> neighbour_location = p_cell_population->GetLocationOfCellCentre(p_neighbour_cell);
                    c_vector<double, DIM> cell_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

                    double dx = neighbour_location[0] - cell_location[0];
                    double dy = neighbour_location[1] - cell_location[1];

                    double deltaX1 = dx - pCell->GetCellData()->GetItem("vecPolaX");
                    double deltaY1 = dy - pCell->GetCellData()->GetItem("vecPolaY");

                    double deltaX2 = -dx - pCell->GetCellData()->GetItem("vecPolaX");
                    double deltaY2 = -dy - pCell->GetCellData()->GetItem("vecPolaY");


                    //regarde si ce vecteur est plus proche de celui polarisé que l'autre (pas de racine, inutile)
                    double Vec1 = deltaX1 * deltaX1 + deltaY1 * deltaY1;
                    double Vec2 = deltaX2 * deltaX2 + deltaY2 * deltaY2;

                    if(Vec1 < Vec2){
                      p_neighbour_cell->RemoveCellProperty<CellEpi>();
                      MAKE_PTR(CellLumen, p_lumen);
                      p_neighbour_cell->AddCellProperty(p_lumen);
                      p_neighbour_cell->GetCellData()->SetItem("mustDie", 0);
                      p_neighbour_cell->GetCellData()->SetItem("vecPolaX", 0);
                      p_neighbour_cell->GetCellData()->SetItem("vecPolaY", 0);
                      p_neighbour_cell->GetCellData()->SetItem("target area", SimulationParameters::SIZE_MIN_LUMEN);
                    }

                    else{
                      pCell->RemoveCellProperty<CellEpi>();
                      MAKE_PTR(CellLumen, p_lumen);
                      pCell->AddCellProperty(p_lumen);
                      pCell->GetCellData()->SetItem("mustDie", 0);
                      pCell->GetCellData()->SetItem("vecPolaX", 0);
                      pCell->GetCellData()->SetItem("vecPolaY", 0);
                      pCell->GetCellData()->SetItem("target area", SimulationParameters::SIZE_MIN_LUMEN);
                    }




                }

              }

            }


            //CellLumen
            bool neighbour_is_lumen = p_neighbour_cell->template HasCellProperty<CellLumen>();
            if(neighbour_is_lumen){
              pCell->GetCellData()->SetItem("lumenNearby",SimulationParameters::NBR_CELL_BETWEEN_TWO_LUMEN);
            }



          }


      }
    }
  }
  //on repasse pour incrémenter les voisins
  for (int i=SimulationParameters::NBR_CELL_BETWEEN_TWO_LUMEN; i>0; --i)
  {

    for (typename AbstractCellPopulation<DIM, DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
      CellPtr pCell = *cell_iter;
      //IF Epithéliale

      if (pCell->HasCellProperty<CellEpi>() && SimulationParameters::NBR_CELL_BETWEEN_TWO_LUMEN-1 > pCell->GetCellData()->GetItem("lumenNearby"))
            {

              VertexBasedCellPopulation<DIM>* p_cell_population = dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) ;
              std::set<unsigned> neighbour_indices = p_cell_population->GetNeighbouringLocationIndices(*cell_iter);

              // If this cell has any neighbours (as defined by mesh/population/interaction distance)...
              if (!neighbour_indices.empty() )
              {
                for (std::set<unsigned>::iterator neighbour_iter = neighbour_indices.begin();
                     neighbour_iter != neighbour_indices.end();
                     ++neighbour_iter)
                {
                  CellPtr p_neighbour_cell = p_cell_population->GetCellUsingLocationIndex(*neighbour_iter);
                  if (p_neighbour_cell->HasCellProperty<CellEpi>())
                  {
                    if(p_neighbour_cell->GetCellData()->GetItem("lumenNearby")-1 >pCell->GetCellData()->GetItem("lumenNearby") && p_neighbour_cell->GetCellData()->GetItem("lumenNearby") > 0)
                    {
                      pCell->GetCellData()->SetItem("lumenNearby",p_neighbour_cell->GetCellData()->GetItem("lumenNearby")-1);
                      //std::cout << pCell->GetCellData()->GetItem("lumenNearby") << std::endl;
                    }
                  }
                }
              }
            }
          }
        }
      }



template<unsigned DIM>
void LumenGenerationModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class LumenGenerationModifier<1>;
template class LumenGenerationModifier<2>;
template class LumenGenerationModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(LumenGenerationModifier)
