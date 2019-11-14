#include "CalibrationErrorModifier.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CellEndo.hpp"
#include "CellLumen.hpp"
#include "CellEpi.hpp"
#include "CellLabel.hpp"
#include "CellPolar.hpp"
#include "CellPeriph.hpp"
#include "RandomNumberGenerator.hpp"
#include "AbstractSimpleGenerationalCellCycleModel.hpp"
#include <iostream>
#include <stdlib.h>
#include <cmath>
using namespace std ;

template<unsigned DIM>
CalibrationErrorModifier<DIM>::CalibrationErrorModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
CalibrationErrorModifier<DIM>::~CalibrationErrorModifier()
{
}

template<unsigned DIM>
void CalibrationErrorModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void CalibrationErrorModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void CalibrationErrorModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

    VertexBasedCellPopulation<DIM>* pCellPopulation = dynamic_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation) ;

    // Permet d'importer les valeurs d'aire issue de HALO
    ifstream inFileArea ;
    double x_area ;
    inFileArea.open("/home/gonayl/Chaste/chaste-src/testoutput/area_error.txt") ;
    std::vector<double> area_error_input;
    if(!inFileArea)
    {
      cout << "Unable to open area_error file in CalibrationErrorWriter" << endl ;
    }
    while(inFileArea >> x_area)
    {
      area_error_input.push_back(x_area) ;
    }
    inFileArea.close() ;

    // Permet d'importer les valeurs de périmetre issue de HALO
    ifstream inFilePerif ;
    double x_perif ;
    inFilePerif.open("/home/gonayl/Chaste/chaste-src/testoutput/area_error.txt") ;
    vector<double> perim_error_input;
    if(!inFilePerif)
    {
      cout << "Unable to open area_error file in CalibrationErrorWriter" << endl ;
    }
    while(inFilePerif >> x_perif)
    {
      perim_error_input.push_back(x_perif) ;
    }
    inFilePerif.close() ;

    // Compute the centre of mass (average position of all the nodes)
    c_vector<double,2> centre_of_mass = zero_vector<double>(2);
    unsigned counter = 0;
    for (unsigned node_index=0; node_index<pCellPopulation->GetNumNodes(); node_index++)
    {
        double x = pCellPopulation->rGetMesh().GetNode(node_index)->rGetLocation()[0];
        double y = pCellPopulation->rGetMesh().GetNode(node_index)->rGetLocation()[1];

        centre_of_mass(0) += x;
        centre_of_mass(1) += y;

        counter++;
    }
    centre_of_mass(0) /= counter;
    centre_of_mass(1) /= counter;

      vector<double> volumes;
      vector<double> distances_x;
      vector<double> distances_y;
      vector<double> bins = {0 , 0.2 , 0.4 , 0.6 , 0.8 , 1.0} ;
      vector<double> pins = {0 , 0.5 , 1.0 , 1.5 , 2.0 , 2.5} ;
      double bin1 = 0; double bin2 = 0; double bin3 = 0; double bin4 = 0; double bin5 = 0;
      double pin1 = 0; double pin2 = 0; double pin3 = 0; double pin4 = 0; double pin5 = 0;


      // Loop over cells and find associated elements so in the same order as the cells in output files
      for (typename AbstractCellPopulation<DIM, DIM>::Iterator cell_iter = pCellPopulation->Begin();
           cell_iter != pCellPopulation->End();
           ++cell_iter)
      {
          double volume = pCellPopulation->GetVolumeOfCell(*cell_iter);
          VertexElement<DIM, DIM>* p_element = pCellPopulation->GetElementCorrespondingToCell(*cell_iter);
          c_vector<double, DIM> cell_location = pCellPopulation->GetLocationOfCellCentre(*cell_iter);
          unsigned p_index = p_element->GetIndex();
          // double perimeter = pCellPopulation->rGetMesh().GetSurfaceAreaOfElement(p_index);
          double perimeter = pCellPopulation->rGetMesh().GetElongationShapeFactorOfElement(p_index);

          double dist_x =  pow((cell_location[0] - centre_of_mass[0]),2) ;
          double dist_y =  pow((cell_location[1] - centre_of_mass[1]),2) ;
          double dist2center = sqrt((cell_location[0] - centre_of_mass[0])*(cell_location[0] - centre_of_mass[0]) + (cell_location[1] - centre_of_mass[1])*(cell_location[1] - centre_of_mass[1])) ;

          volumes.push_back(volume) ;
          distances_x.push_back(dist_x) ;
          distances_y.push_back(dist_y) ;

          if (volume < bins[1])
          {
            bin1++ ;
          }
          else if (volume < bins[2] && volume >= bins[1])
          {
            bin2++ ;
          }
          else if (volume < bins[3] && volume >= bins[2])
          {
            bin3++ ;
          }
          else if (volume < bins[4] && volume >= bins[3])
          {
            bin4++ ;
          }
          else if (volume < bins[5] && volume >= bins[4])
          {
            bin5++ ;
          }

          if  (perimeter < pins[1])
          {
            pin1++ ;
          }
          else if (perimeter < pins[2] && perimeter  >= pins[1])
          {
            pin2++ ;
          }
          else if (perimeter < pins[3] && perimeter  >= pins[2])
          {
            pin3++ ;
          }
          else if (perimeter < pins[4] && perimeter  >= pins[3])
          {
            pin4++ ;
          }
          else if (perimeter < pins[5] && perimeter  >= pins[4])
          {
            pin5++ ;
          }

          cell_iter->GetCellData()->SetItem("dist2center",dist2center) ;

      }

      double volume_total = accumulate(volumes.begin(), volumes.end(), 0.0) ;
      //double volume_average = volume_total/volumes.size();
      double distance_variance_x = accumulate(distances_x.begin(), distances_x.end(), 0.0)/distances_x.size();
      double distance_variance_y = accumulate(distances_y.begin(), distances_y.end(), 0.0)/distances_y.size();

      // std::cout << volumes.size() << " / " << volume_average << std::endl ;
      vector<double> occur_area = {bin1 , bin2 , bin3 , bin4 , bin5} ;
      vector<double> occur_perif = {pin1 ,pin2 , pin3 , pin4 , pin5} ;
      vector<double> occur_area2 = area_error_input ;
      vector<double> occur_perif2 = area_error_input ;

      double error_area = pow((occur_area[0] - occur_area2[0]),2) + pow((occur_area[1] - occur_area2[1]),2)
       + pow((occur_area[2] - occur_area2[2]),2) + pow((occur_area[3] - occur_area2[3]),2) + pow((occur_area[4] - occur_area2[4]),2)  ;

      double error_perif = pow((occur_perif[0] - occur_perif2[0]),2) + pow((occur_perif[1] - occur_perif2[1]),2)
       + pow((occur_perif[2] - occur_perif2[2]),2) + pow((occur_perif[3] - occur_perif2[3]),2) + pow((occur_perif[4] - occur_perif2[4]),2)  ;



      /*cout << "Error area = " << error_area << " Error perif = " << error_perif
       << " Error pos x = " << distance_variance_x << " Error pos y = "
       << distance_variance_y << " Volume total = " << volume_total << endl;*/
      // std::cout << occur_area << std::endl ;

      /* const float bucket_size = 0.2 ;
      double number_of_buckets = (double)ceil(1/bucket_size) ;
      vector<double> histogram(nomber_of_buckets) */

      cout << " " << error_area << " " << error_perif << " " << distance_variance_x << " " << distance_variance_y << " " << volume_total << endl;

}

template<unsigned DIM>
void CalibrationErrorModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class CalibrationErrorModifier<1>;
template class CalibrationErrorModifier<2>;
template class CalibrationErrorModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(CalibrationErrorModifier)
