#include "ObstructionBoundaryCondition.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CellStalk.hpp"
#include "CellTip.hpp"
#include "CellVessel.hpp"
#include "CellPeriph.hpp"
#include "CounterSingletonRepulsion.hpp"
#include <stdlib.h>
#include<cmath>
#include <iostream>
#include <fstream>

using namespace std ;

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
ObstructionBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::ObstructionBoundaryCondition(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>* pCellPopulation)
        : AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM>(pCellPopulation)
{

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ObstructionBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::ImposeBoundaryCondition(const std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> >& rOldLocations)
{
    if (dynamic_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(this->mpCellPopulation)==nullptr)
    {
        EXCEPTION("ObstructionBoundaryCondition requires a subclass of AbstractOffLatticeCellPopulation.");
    }

    const int n_columns = 4;
    double endocells[][n_columns]{} ;
    int endo_count = 0 ;

    int endo_count_local = 0 ;

    int endo_count_global = CounterSingletonRepulsion::Instance()->GetCount() ;

    //cout << endo_count_global << "endo stalk fixed cell(s) "  << endl ;


    for (typename AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)

    {
      if (cell_iter->template HasCellProperty<CellStalk>() && cell_iter->template HasCellProperty<CellPeriph>())
      {

        endo_count_local++ ;

      }

    }


    if (endo_count_global != endo_count_local)
    {

      cout << "New endo cell" << endl ;

      // PHASE 0, computing mass center (will not be used after)

      vector<double> moy_x_pos;
      vector<double> moy_y_pos;


      for (typename AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
           cell_iter != this->mpCellPopulation->End();
           ++cell_iter)
      {
        c_vector<double, SPACE_DIM> cell_location = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter) ;
        moy_x_pos.push_back(cell_location[0]);
        moy_y_pos.push_back(cell_location[1]);

      }

      double xsum = accumulate(moy_x_pos.begin(), moy_x_pos.end(), 0.0);
      double ysum = accumulate(moy_y_pos.begin(), moy_y_pos.end(), 0.0);
      double xmoy = xsum / moy_x_pos.size() ;
      double ymoy = ysum / moy_y_pos.size() ;

      // For the description of phases, see lab notes pages xx

      // PHASE 1
      // iterate over cell AddPopulationWriter
    //double endocells[][4]{} ;


    CounterSingletonRepulsion::Instance()->ResetCountToZero() ;
    //int endo_count = 0 ;



    for (typename AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)

    {

      if ( cell_iter->template HasCellProperty<CellStalk>() && cell_iter->template HasCellProperty<CellPeriph>())
      {

        c_vector<double, 2> direct_vector ;

        c_vector<double, 2> centre_of_cell = this->mpCellPopulation->GetLocationOfCellCentre(*cell_iter);

        direct_vector[0] = centre_of_cell[0] - xmoy ;
        direct_vector[1] = centre_of_cell[1] - ymoy ;

        endocells[endo_count][0] =  direct_vector[1] ; // a
        endocells[endo_count][1] =  - direct_vector[0] ; // b
        endocells[endo_count][2] = - endocells[endo_count][0]*(centre_of_cell[0] - 0.4) - endocells[endo_count][1]*(centre_of_cell[1] - 0.4) ; // l
        endocells[endo_count][3] = - endocells[endo_count][0]*(centre_of_cell[0] + 0.4) - endocells[endo_count][1]*(centre_of_cell[1] + 0.4) ; // r

        endo_count++ ;

      }

    }

    double endo_count_previous = CounterSingletonRepulsion::Instance()->GetCount();
    //cout << endo_count_current << "endo cells" << endl ;

    //<< endocells[0][0] << ", b = " << endocells[0][1]
    //<< " and l and r = " << endocells[0][2] << " ; " <<  endocells[0][3] << endl ;



    // int n_line = endo_count_current - 1 ;
    ofstream myfile ("example.txt", ios::app);
    for (int i_endo = endo_count_previous; i_endo < endo_count ; i_endo++ )
    {
      if (i_endo == 0)
      {
        if (myfile.is_open())
        {
          //cout << "creating file" << endl ;
          for(int count = 0; count < n_columns; count ++)
          {
            myfile << endocells[i_endo][count] << " " ;
          }
        }
      }
      else
      {
        if (myfile.is_open())
        {
          //cout << "updating" << endl ;
          myfile << endl ;
          for(int count = 0; count < n_columns; count ++)
          {
            myfile << endocells[i_endo][count] << " " ;
          }
        }
      } // end else
    } // end for

    for (int i_counter = 0; i_counter < endo_count ; i_counter++)
    {
      CounterSingletonRepulsion::Instance()->IncrementCounter();
    }

    /*

    if (endo_count_current == 1 )
    {
      ofstream myfile ("example.txt", ios::app);
      if (myfile.is_open())
      {
        cout << "creating file" << endl ;
        for(int count = 0; count < n_columns; count ++){
          myfile << endocells[n_line][count] << " " ;
        }

        myfile.close();
      }
      else cout << "Unable to open file";
    }
    else if (endo_count_current > 1)
    {

      if (myfile.is_open())
      {
        cout << "updating file" << endl ;
        myfile << endl ;
        for(int count = 0; count < n_columns; count ++){
          myfile << endocells[n_line][count] << " " ;
        }

        myfile.close();
      }
      else cout << "Unable to open file";
    }
    else
    {
      cout << "No endo stalk cell" << endl ;
    }
    */

    }

    else
    {
      //cout << "No new endo cell" << endl ;
    }





    //cout << endo_count << "endo stalk fixed cell(s), a = " << endocells[0][0] << ", b = " << endocells[0][1]
    //<< " and l and r = " << endocells[0][2] << " ; " <<  endocells[0][3] << endl ;


    assert((dynamic_cast<AbstractCentreBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(this->mpCellPopulation))
            || (SPACE_DIM==ELEMENT_DIM && (dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>*>(this->mpCellPopulation))) );


            ifstream inFile ;
            double x ;
            inFile.open("example.txt") ;
            std::vector<double> endo_cells_fromfile;

            if(!inFile)
            {
              cout << "Unable to open file" << endl ;
            }
            while(inFile >> x)
            {
              endo_cells_fromfile.push_back(x) ;
            }
            inFile.close() ;

            //vector<double> endoLoc;
            if (SPACE_DIM != 1)
            {


            assert(SPACE_DIM == ELEMENT_DIM);
            assert(dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>*>(this->mpCellPopulation));

            // Iterate over all nodes and update their positions according to the boundary conditions
            unsigned num_nodes = this->mpCellPopulation->GetNumNodes();
            int n_lines_endo = CounterSingletonRepulsion::Instance()->GetCount() ;
            for (unsigned node_index=0; node_index<num_nodes; node_index++)
            {
                Node<SPACE_DIM>* p_node = this->mpCellPopulation->GetNode(node_index);
                c_vector<double, SPACE_DIM> node_location = p_node->rGetLocation();
                double y_coordinate = p_node->rGetLocation()[1];
                double x_coordinate = p_node->rGetLocation()[0];

                // "left" limit equation (l) : ax + by - l = 0
                // "right" limit equation : ax + by = r

                // normal to limit eq = dx + ey + f = 0


                for(int endo_index = 0; endo_index < n_lines_endo; endo_index++)
                {

                int pos_in_vec = endo_index  * 4 ;

                double a = endo_cells_fromfile[0 + pos_in_vec];
                double b = endo_cells_fromfile[1 + pos_in_vec];
                double l = endo_cells_fromfile[2 + pos_in_vec];
                double r = endo_cells_fromfile[3 + pos_in_vec];

                // cout << a << " ; " << b << " ; " << l << " ; " << r << " ; " <<endl ;

                double d = b ;
                double e = -a ;
                double f = (-b * x_coordinate) + (a * y_coordinate) ;

                // now that we have the equation for the normal line to the left limit, we need the closest point on this line to the current node (orthogonal projection)

                double x_proj_l ;
                double y_proj_l ;

                double x_proj_r ;
                double y_proj_r ;

                y_proj_l = (((d*l)/a)- f ) / (((-d*b)/a) + e) ;
                x_proj_l = (-b*y_proj_l - l) / a ;

                y_proj_r = (((d*r)/a)- f ) / (((-d*b)/a) + e) ;
                x_proj_r = (-b*y_proj_r - r) / a ;


                double dist1 = sqrt(pow((x_coordinate - x_proj_l),2) + pow((y_coordinate - y_proj_l),2)) ;
                double dist2 = sqrt(pow((x_coordinate - x_proj_r),2) + pow((y_coordinate - y_proj_r),2)) ;

                // c_vector<double, SPACE_DIM> old_node_location = rOldLocations.find(p_node)->second; // node previous location

                if (p_node->IsBoundaryNode() && (b*y_coordinate + a*x_coordinate) > -l && (b*y_coordinate + a*x_coordinate) < -r && dist1 < dist2  )
                {
                  p_node->rGetModifiableLocation()[0] = x_proj_l;
                  p_node->rGetModifiableLocation()[1] = y_proj_l;
                  //p_node->rGetModifiableLocation() = old_node_location;
                }
                else if (p_node->IsBoundaryNode() && (b*y_coordinate + a*x_coordinate) > -l && (b*y_coordinate + a*x_coordinate) < -r && dist1 > dist2  )
                {
                  p_node->rGetModifiableLocation()[0] = x_proj_r;
                  p_node->rGetModifiableLocation()[1] = y_proj_r ;
                  //p_node->rGetModifiableLocation() = old_node_location;
                }
                }

              }

    }
    else
    {
        // SPACE_DIM == 1
        NEVER_REACHED;
        //ObstructionBoundaryCondition::ImposeBoundaryCondition is not implemented in 1D
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool ObstructionBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::VerifyBoundaryCondition()
{
    return true;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void ObstructionBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::OutputCellPopulationBoundaryConditionParameters(out_stream& rParamsFile)
{

    // Call method on direct parent class
    AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM>::OutputCellPopulationBoundaryConditionParameters(rParamsFile);
}

// Explicit instantiation
template class ObstructionBoundaryCondition<1,1>;
template class ObstructionBoundaryCondition<1,2>;
template class ObstructionBoundaryCondition<2,2>;
template class ObstructionBoundaryCondition<1,3>;
template class ObstructionBoundaryCondition<2,3>;
template class ObstructionBoundaryCondition<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(ObstructionBoundaryCondition)
