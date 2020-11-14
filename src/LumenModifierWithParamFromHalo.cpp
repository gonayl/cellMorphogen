#include "LumenModifierWithParamFromHalo.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexElement.hpp"
#include <stdlib.h>
#include <math.h>
//#include "ApoptoticCellProperty.hpp"
#include "CellEpi.hpp"
#include "CellLumen.hpp"
#include "SimulationParameters.hpp"

template<unsigned DIM>
LumenModifierWithParamFromHalo<DIM>::LumenModifierWithParamFromHalo()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
LumenModifierWithParamFromHalo<DIM>::~LumenModifierWithParamFromHalo()
{
}

template<unsigned DIM>
void LumenModifierWithParamFromHalo<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void LumenModifierWithParamFromHalo<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void LumenModifierWithParamFromHalo<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();
    //VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);

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
        EXCEPTION("Lumen Modifier is to be used with a VertexBasedCellPopulation only");
    }

    //calcul de la taille normalisée
    double tailleVecteur = 0;
    int nombreVecteur = 0;
    for (typename AbstractCellPopulation<DIM, DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {

      CellPtr pCell = *cell_iter;

      if (pCell->HasCellProperty<CellLumen>())
      {


        //on cherche la cellule fille
        VertexBasedCellPopulation<DIM>* p_cell_population = dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) ;
        std::set<unsigned> neighbour_indices = p_cell_population->GetNeighbouringLocationIndices(*cell_iter);

        // If this cell has any neighbours (as defined by mesh/population/interaction distance)...
        if (!neighbour_indices.empty() )
        {

            for (std::set<unsigned>::iterator neighbour_iter = neighbour_indices.begin();
                 neighbour_iter != neighbour_indices.end();
                 ++neighbour_iter)
            {
              CellPtr pnCell = p_cell_population->GetCellUsingLocationIndex(*neighbour_iter);
              nombreVecteur = nombreVecteur +1;

              double xl = pnCell->GetCellData()->GetItem("vecPolaX");
              double yl = pnCell->GetCellData()->GetItem("vecPolaY");

              double normeVec = sqrt(xl*xl+yl*yl);
              tailleVecteur = tailleVecteur+normeVec;


            }
          }
        }
      }
      double apportParUniteVecteur = 0;
      if(nombreVecteur > 0)
      {
        double tailleMoyenne = tailleVecteur/ nombreVecteur;
        apportParUniteVecteur = SimulationParameters::LUMEN_SIZE_FACTOR / tailleMoyenne;
      }


    //pour la taille
    int nbrLumen =0;
    double sizeLumenByEpi = 0;
    for (typename AbstractCellPopulation<DIM, DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {

      CellPtr pCell = *cell_iter;

      if (pCell->HasCellProperty<CellLumen>())
      {
        nbrLumen = nbrLumen +1 ;
        int nbrCellEpi = 0;

        double targetSize = 0;
        if(pCell->GetAge() < SimulationParameters::AGE_TO_LUMEN_MATURITY)
        {
          targetSize = (1 - (pCell->GetAge() / SimulationParameters::AGE_TO_LUMEN_MATURITY)) *SimulationParameters::SIZE_MIN_LUMEN;
        }




        //on cherche la cellule fille
        VertexBasedCellPopulation<DIM>* p_cell_population = dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) ;
        std::set<unsigned> neighbour_indices = p_cell_population->GetNeighbouringLocationIndices(*cell_iter);

        // If this cell has any neighbours (as defined by mesh/population/interaction distance)...
        if (!neighbour_indices.empty() )
        {

            for (std::set<unsigned>::iterator neighbour_iter = neighbour_indices.begin();
                 neighbour_iter != neighbour_indices.end();
                 ++neighbour_iter)
            {

              CellPtr pnCell = p_cell_population->GetCellUsingLocationIndex(*neighbour_iter);

              bool neighbour_is_epi = pnCell->template HasCellProperty<CellEpi>();
              if(neighbour_is_epi)
              {

                c_vector<double, DIM> neighbour_location = p_cell_population->GetLocationOfCellCentre(pnCell);
                c_vector<double, DIM> cell_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);

                double dx = neighbour_location[0] - cell_location[0];
                double dy = neighbour_location[1] - cell_location[1];

                double deltaX1 = dx - pnCell->GetCellData()->GetItem("vecPolaX");
                double deltaY1 = dy - pnCell->GetCellData()->GetItem("vecPolaY");

                double deltaX2 = -dx - pnCell->GetCellData()->GetItem("vecPolaX");
                double deltaY2 = -dy - pnCell->GetCellData()->GetItem("vecPolaY");


                //regarde si ce vecteur est plus proche de celui polarisé que l'autre (pas de racine, inutile)
                double Vec1 = deltaX1 * deltaX1 + deltaY1 * deltaY1;
                double Vec2 = deltaX2 * deltaX2 + deltaY2 * deltaY2;

                if(Vec1 > Vec2){

                  double xl = pnCell->GetCellData()->GetItem("vecPolaX");
                  double yl = pnCell->GetCellData()->GetItem("vecPolaY");

                  double norme = sqrt(xl*xl + yl*yl);


                  double coefSize = 1;

                  if(pCell->GetAge() < SimulationParameters::AGE_TO_LUMEN_MATURITY)
                  {
                    coefSize = pCell->GetAge() / SimulationParameters::AGE_TO_LUMEN_MATURITY;
                  }
                  nbrCellEpi = nbrCellEpi + 1;
                  targetSize = targetSize + norme * apportParUniteVecteur * coefSize;
                }

              }


              bool neighbour_is_lumen = pnCell->template HasCellProperty<CellLumen>();
              if(neighbour_is_lumen){
                std::cout << "Kill a lumen because neighbours of a lumen" << '\n';
                pnCell->GetCellData()->SetItem("target area",0.01);
              }
            }
          }
          targetSize = targetSize * targetSize;
          pCell->GetCellData()->SetItem("target area",targetSize);

          sizeLumenByEpi = sizeLumenByEpi + sqrt(targetSize*80) / nbrCellEpi;

    }




  }
  /*if (nbrLumen > 0)
  {
    std::cout << sizeLumenByEpi/nbrLumen << '\n';
  }
  */
}



template<unsigned DIM>
void LumenModifierWithParamFromHalo<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class LumenModifierWithParamFromHalo<1>;
template class LumenModifierWithParamFromHalo<2>;
template class LumenModifierWithParamFromHalo<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(LumenModifierWithParamFromHalo)
