#include "CalibrationErrorWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include <iostream>
#include <stdlib.h>
#include <cmath>
using namespace std ;

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CalibrationErrorWriter<ELEMENT_DIM, SPACE_DIM>::CalibrationErrorWriter()
    : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>("calibrationerror.dat")
{

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CalibrationErrorWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("CalibrationErrorWriter cannot be used with a MeshBasedCellPopulation");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CalibrationErrorWriter<ELEMENT_DIM, SPACE_DIM>::Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("CalibrationErrorWriter cannot be used with a CaBasedCellPopulation");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CalibrationErrorWriter<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("CalibrationErrorWriter cannot be used with a NodeBasedCellPopulation");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CalibrationErrorWriter<ELEMENT_DIM, SPACE_DIM>::Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("CalibrationErrorWriter cannot be used with a PottsBasedCellPopulation");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CalibrationErrorWriter<ELEMENT_DIM, SPACE_DIM>::Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{

  // Permet d'importer les valeurs d'aire issue de HALO
  ifstream inFileArea ;
  double x_area ;
  inFileArea.open("/home/gonayl/Chaste/chaste-src/testoutput/area_error.txt") ;
  vector<double> area_error_input;
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
    vector<double> bins = {26, 51, 76, 101, 126, 151, 176, 201, 226, 251, 276} ;
    vector<double> pins = {0 , 0.5 , 1.0 , 1.5 , 2.0 , 2.5} ;
    double bin1 = 0; double bin2 = 0; double bin3 = 0; double bin4 = 0; double bin5 = 0;
    double bin6 = 0; double bin7 = 0; double bin8 = 0; double bin9 = 0; double bin10 = 0;
    double bin11 = 0; double bin12 = 0;

    double pin1 = 0; double pin2 = 0; double pin3 = 0; double pin4 = 0; double pin5 = 0;


    // Loop over cells and find associated elements so in the same order as the cells in output files
    for (typename AbstractCellPopulation<SPACE_DIM, SPACE_DIM>::Iterator cell_iter = pCellPopulation->Begin();
         cell_iter != pCellPopulation->End();
         ++cell_iter)
    {
        double volume = pCellPopulation->GetVolumeOfCell(*cell_iter);
        VertexElement<SPACE_DIM, SPACE_DIM>* p_element = pCellPopulation->GetElementCorrespondingToCell(*cell_iter);
        c_vector<double, SPACE_DIM> cell_location = pCellPopulation->GetLocationOfCellCentre(*cell_iter);
        unsigned p_index = p_element->GetIndex();
        double cell_perimeter= pCellPopulation->rGetMesh().GetSurfaceAreaOfElement(p_index);
        // double cell_perimeter = pCellPopulation->rGetMesh().GetElongationShapeFactorOfElement(p_index);

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
        else if (volume < bins[6] && volume >= bins[5])
        {
          bin6++ ;
        }
        else if (volume < bins[7] && volume >= bins[6])
        {
          bin7++ ;
        }
        else if (volume < bins[8] && volume >= bins[7])
        {
          bin8++ ;
        }
        else if (volume < bins[9] && volume >= bins[8])
        {
          bin9++ ;
        }
        else if (volume < bins[10] && volume >= bins[9])
        {
          bin10++ ;
        }
        else if (volume < bins[11] && volume >= bins[10])
        {
          bin11++ ;
        }
        else if (volume >= bins[11])
        {
          bin12++ ;
        }

        if  (cell_perimeter< pins[1])
        {
          pin1++ ;
        }
        else if (cell_perimeter< pins[2] && cell_perimeter >= pins[1])
        {
          pin2++ ;
        }
        else if (cell_perimeter< pins[3] && cell_perimeter >= pins[2])
        {
          pin3++ ;
        }
        else if (cell_perimeter< pins[4] && cell_perimeter >= pins[3])
        {
          pin4++ ;
        }
        else if (cell_perimeter< pins[5] && cell_perimeter >= pins[4])
        {
          pin5++ ;
        }

        cell_iter->GetCellData()->SetItem("dist2center",dist2center) ;

    }

    double volume_total = accumulate(volumes.begin(), volumes.end(), 0.0) ;
    //double volume_total_HALO = 114900;
    //double error_volume = pow((volume_total - volume_total_HALO),2) ;
    //double volume_average = volume_total/volumes.size();
    double distance_variance_x = accumulate(distances_x.begin(), distances_x.end(), 0.0)/distances_x.size();
    double distance_variance_y = accumulate(distances_y.begin(), distances_y.end(), 0.0)/distances_y.size();

/*
    // std::cout << volumes.size() << " / " << volume_average << std::endl ;
    vector<double> occur_area = {bin1 , bin2 , bin3 , bin4 , bin5, bin6 , bin7 , bin8 , bin9 , bin10, bin11, bin12} ;
    vector<double> occur_perif = {pin1 ,pin2 , pin3 , pin4 , pin5} ;
    vector<double> occur_area2 = {0,   173,   419,   377,   219,    95,    35,    19,     4,     0,     1,     1} ;
    vector<double> occur_perif2 = area_error_input ;
    double error_area_dist = pow((occur_area[0] - occur_area2[0]),2) + pow((occur_area[1] - occur_area2[1]),2)
     + pow((occur_area[2] - occur_area2[2]),2) + pow((occur_area[3] - occur_area2[3]),2) + pow((occur_area[4] - occur_area2[4]),2)
     + pow((occur_area[5] - occur_area2[5]),2) + pow((occur_area[6] - occur_area2[6]),2) + pow((occur_area[7] - occur_area2[7]),2)
     + pow((occur_area[8] - occur_area2[8]),2) + pow((occur_area[9] - occur_area2[9]),2) + pow((occur_area[10] - occur_area2[10]),2)
     + pow((occur_area[11] - occur_area2[11]),2) + pow((occur_area[12] - occur_area2[12]),2) ;

     */
     /*
    double error_perif = pow((occur_perif[0] - occur_perif2[0]),2) + pow((occur_perif[1] - occur_perif2[1]),2)
    + pow((occur_perif[2] - occur_perif2[2]),2) + pow((occur_perif[3] - occur_perif2[3]),2) + pow((occur_perif[4] - occur_perif2[4]),2)  ;
    */


    /*cout << "Error area = " << error_area_dist << " Error perif = " << error_perif
     << " Error pos x = " << distance_variance_x << " Error pos y = "
     << distance_variance_y << " Volume total = " << volume_total << endl;*/
    // std::cout << occur_area << std::endl ;

    /* const float bucket_size = 0.2 ;
    double number_of_buckets = (double)ceil(1/bucket_size) ;
    vector<double> histogram(nomber_of_buckets) */

    //cout << "error volume total = " << error_volume << endl;

    *this->mpOutStream << " " << volume_total << " " << distance_variance_x << " " << distance_variance_y << " " ;
}

// Explicit instantiation
template class CalibrationErrorWriter<1,1>;
template class CalibrationErrorWriter<1,2>;
template class CalibrationErrorWriter<2,2>;
template class CalibrationErrorWriter<1,3>;
template class CalibrationErrorWriter<2,3>;
template class CalibrationErrorWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CalibrationErrorWriter)
