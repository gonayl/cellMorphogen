#include "PolarisationModifier.hpp"
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
PolarisationModifier<DIM>::PolarisationModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
      mVecPolarisationDecrease(0.075),
      mEpiEpiPolarisationParameter(-0.08),
      mEndoEpiPolarisationParameter(0.24),
      mLumenEpiPolarisationParameter(-0.15)
{
}

template<unsigned DIM>
PolarisationModifier<DIM>::~PolarisationModifier()
{
}

template<unsigned DIM>
void PolarisationModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void PolarisationModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void PolarisationModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
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
        double vecPolaX = pCell->GetCellData()->GetItem("vecPolaX");
        double vecPolaY = pCell->GetCellData()->GetItem("vecPolaY");


        vecPolaX = vecPolaX - vecPolaX * this->GetVecPolarisationDecrease() * SimulationTime::Instance()->GetTimeStep() * 2;
        vecPolaY = vecPolaY - vecPolaY * this->GetVecPolarisationDecrease() * SimulationTime::Instance()->GetTimeStep() * 2;




        //impact des voisins

        VertexBasedCellPopulation<DIM>* p_cell_population = dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) ;
        std::set<unsigned> neighbour_indices = p_cell_population->GetNeighbouringLocationIndices(*cell_iter);

        c_vector<double, DIM> cell_location = rCellPopulation.GetLocationOfCellCentre(*cell_iter);


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


                //Cell endothéliale
                bool neighbour_is_endo = p_neighbour_cell->template HasCellProperty<CellEndo>();

                //Cell Lumen
                bool neighbour_is_lumen = p_neighbour_cell->template HasCellProperty<CellLumen>();


                c_vector<double, DIM> neighbour_location = p_cell_population->GetLocationOfCellCentre(p_neighbour_cell);

                double dx = cell_location[0] - neighbour_location[0];
                double dy = cell_location[1] - neighbour_location[1];


                double normalisationInterCell = sqrt(dx*dx + dy*dy);


                if ( neighbour_is_epiDir == 1)
                {
                  //std::cout << "coucou" << std::endl;

                  double vecPolaXNeighbour = p_neighbour_cell->GetCellData()->GetItem("vecPolaX");
                  double vecPolaYNeighbour = p_neighbour_cell->GetCellData()->GetItem("vecPolaY");


                  double prodScalaire = (dx*vecPolaXNeighbour+vecPolaYNeighbour*dy);

                  double normeVecPolaNeighbour = sqrt(vecPolaXNeighbour * vecPolaXNeighbour + vecPolaYNeighbour * vecPolaYNeighbour);

                  //cos angle fr.wikihow.com/calculer-l'angle-entre-deux-vecteurs
                  double cosAngle = prodScalaire/(normeVecPolaNeighbour*normalisationInterCell);
                  cosAngle = sqrt(cosAngle*cosAngle);//valeur absolue
                  //std::cout << cosAngle << "  " << normalisationInterCell << " " << prodScalaire << "  " << normeVecPolaNeighbour << '\n';
                  if(vecPolaXNeighbour > 0){
                    vecPolaX = vecPolaX + this->GetEpiEpiPolarisationParameter() * SimulationTime::Instance()->GetTimeStep() * 2*cosAngle;
                  }
                  else if(vecPolaXNeighbour < 0){
                    vecPolaX = vecPolaX - this->GetEpiEpiPolarisationParameter() * SimulationTime::Instance()->GetTimeStep() * 2*cosAngle;
                  }

                  if(vecPolaYNeighbour > 0){
                    vecPolaY = vecPolaY + this->GetEpiEpiPolarisationParameter() * SimulationTime::Instance()->GetTimeStep() * 2*cosAngle;
                  }
                  else if(vecPolaYNeighbour < 0){
                    vecPolaY = vecPolaY - this->GetEpiEpiPolarisationParameter() * SimulationTime::Instance()->GetTimeStep() * 2*cosAngle;
                  }
                }

                if(neighbour_is_endo == 1){



                  vecPolaX = vecPolaX + dx / normalisationInterCell * this->GetEndoEpiPolarisationParameter() * SimulationTime::Instance()->GetTimeStep() * 2;
                  vecPolaY = vecPolaY + dy / normalisationInterCell * this->GetEndoEpiPolarisationParameter() * SimulationTime::Instance()->GetTimeStep() * 2;

                  //std::cout << vecPolaX << std::endl;

                }


                if(neighbour_is_lumen == 1 && p_neighbour_cell->GetCellData()->GetItem("mustDie")==0){

                  vecPolaX = vecPolaX + dx / normalisationInterCell * this->GetLumenEpiPolarisationParameter() * SimulationTime::Instance()->GetTimeStep() * 2;
                  vecPolaY = vecPolaY + dy / normalisationInterCell * this->GetLumenEpiPolarisationParameter() * SimulationTime::Instance()->GetTimeStep() * 2;

                  //std::cout << vecPolaX << std::endl;

                }

            }

        }
        //polarise la périphérie

        VertexElement<DIM,DIM>* cell_pos = p_cell_population->GetElementCorrespondingToCell(pCell);

        unsigned num_nodes_in_element = cell_pos->GetNumNodes();

        //calcul du vecteur
        double dx = 0;
        double dy = 0;
        for (unsigned i=0; i<num_nodes_in_element; i++)
        {

          if(cell_pos->GetNode(i)->IsBoundaryNode())
          {
            c_vector<double, DIM> node_pos = cell_pos->GetNodeLocation(i);
            dx = dx + cell_location[0] - node_pos[0];
            dy = dy + cell_location[1] - node_pos[1];
          }

        }
        if(dx != 0 || dy !=0)
        {
          double normalisation = sqrt(dx*dx + dy*dy);

          vecPolaX = vecPolaX + dx / normalisation * SimulationParameters::IMPACT_POLARISATION_PERIPH_ON_EPI * SimulationTime::Instance()->GetTimeStep() * 2;
          vecPolaY = vecPolaY + dy / normalisation * SimulationParameters::IMPACT_POLARISATION_PERIPH_ON_EPI * SimulationTime::Instance()->GetTimeStep() * 2;

        }



        //on assigne
        pCell->GetCellData()->SetItem("vecPolaX",vecPolaX);
        pCell->GetCellData()->SetItem("vecPolaY",vecPolaY);

      }


    }

}

template<unsigned DIM>
double PolarisationModifier<DIM>::GetVecPolarisationDecrease()
{
    return mVecPolarisationDecrease;
}
template<unsigned DIM>
double PolarisationModifier<DIM>::GetEpiEpiPolarisationParameter()
{
    return mEpiEpiPolarisationParameter;
}

template<unsigned DIM>
double PolarisationModifier<DIM>::GetEndoEpiPolarisationParameter()
{
    return mEndoEpiPolarisationParameter;
}

template<unsigned DIM>
double PolarisationModifier<DIM>::GetLumenEpiPolarisationParameter()
{
    return mLumenEpiPolarisationParameter;
}

//----*----//

template<unsigned DIM>
void PolarisationModifier<DIM>::SetVecPolarisationDecrease(double vecPolarisationDecrease)
{
    mVecPolarisationDecrease = vecPolarisationDecrease;
}

template<unsigned DIM>
void PolarisationModifier<DIM>::SetEpiEpiPolarisationParameter(double epiEpiPolarisationParameter)
{
    mEpiEpiPolarisationParameter = epiEpiPolarisationParameter;
}

template<unsigned DIM>
void PolarisationModifier<DIM>::SetEndoEpiPolarisationParameter(double endoEpiPolarisationParameter)
{
    mEndoEpiPolarisationParameter = endoEpiPolarisationParameter;
}

template<unsigned DIM>
void PolarisationModifier<DIM>::SetLumenEpiPolarisationParameter(double lumenEpiPolarisationParameter)
{
    mLumenEpiPolarisationParameter = lumenEpiPolarisationParameter;
}

template<unsigned DIM>
void PolarisationModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // Output member variables
    *rParamsFile << "\t\t\t<VecPolarisationDecrease" << mVecPolarisationDecrease << "/VecPolarisationDecrease> \n" ;

    *rParamsFile << "\t\t\t<EpiEpiPolarisationParameter" << mEpiEpiPolarisationParameter << "</EpiEpiPolarisationParameter> \n";

    *rParamsFile << "\t\t\t<EndoEpiPolarisationParameter" << mEndoEpiPolarisationParameter << "</EndoEpiPolarisationParameter> \n";

    *rParamsFile << "\t\t\t<LumenEpiPolarisationParameter" << mLumenEpiPolarisationParameter << "</LumenEpiPolarisationParameter> \n";

    //Call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class PolarisationModifier<1>;
template class PolarisationModifier<2>;
template class PolarisationModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PolarisationModifier)
