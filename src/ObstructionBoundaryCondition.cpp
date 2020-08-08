#include "ObstructionBoundaryCondition.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CellStalk.hpp"
#include "CellTip.hpp"
#include "CellEndo.hpp"
#include "CellVessel.hpp"
#include "CellPeriph.hpp"
#include "CounterSingletonRepulsion.hpp"
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>

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



    c_vector<double, SPACE_DIM> population_center = this->mpCellPopulation->GetCentroidOfCellPopulation();



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
            if (SPACE_DIM != 1 && endo_cells_fromfile.size() > 0)
            {


            assert(SPACE_DIM == ELEMENT_DIM);
            assert(dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>*>(this->mpCellPopulation));

            // Iterate over all nodes and update their positions according to the boundary conditions
            unsigned num_nodes = this->mpCellPopulation->GetNumNodes();
            //cout << num_nodes << endl ;
            //int n_lines_endo = CounterSingletonRepulsion::Instance()->GetCount() ;
            int n_lines_endo = endo_cells_fromfile.size()/6 ;


            for (unsigned node_index=0; node_index<num_nodes; node_index++)
            {
                Node<SPACE_DIM>* p_node = this->mpCellPopulation->GetNode(node_index);
                c_vector<double, SPACE_DIM> node_location = p_node->rGetLocation();

                std::set<unsigned> elements_containing_node = p_node->rGetContainingElementIndices();



                int is_endo = 0 ;

                for (std::set<unsigned>::iterator element_index = elements_containing_node.begin();
                     element_index != elements_containing_node.end();
                     ++element_index)
                {
                  CellPtr p_cell = this->mpCellPopulation->GetCellUsingLocationIndex(*element_index);
                  if (p_cell->HasCellProperty<CellEndo>())
                  {
                    is_endo++ ;
                  }
                }



                double y_coordinate = p_node->rGetLocation()[1];
                double x_coordinate = p_node->rGetLocation()[0];

                // "left" limit equation (l) : ax + by - l = 0
                // "right" limit equation : ax + by = r

                // normal to limit eq = dx + ey + f = 0



                for(int endo_index = 0; endo_index < n_lines_endo; endo_index++)
                {

                int pos_in_vec = endo_index  * 6 ;

                double a = endo_cells_fromfile[0 + pos_in_vec];
                double b = endo_cells_fromfile[1 + pos_in_vec];
                double l = endo_cells_fromfile[2 + pos_in_vec];
                double r = endo_cells_fromfile[3 + pos_in_vec];
                c_vector<double, SPACE_DIM> cell_center ;
                cell_center[0] = endo_cells_fromfile[4 + pos_in_vec];
                cell_center[1] = endo_cells_fromfile[5 + pos_in_vec];

                c_vector<double, SPACE_DIM> p_node_position = p_node->rGetLocation();

                c_vector<double, SPACE_DIM> vec_from_cell_to_node = this->mpCellPopulation->rGetMesh().GetVectorFromAtoB(p_node_position, cell_center);
                c_vector<double, SPACE_DIM> vec_from_center_to_node = this->mpCellPopulation->rGetMesh().GetVectorFromAtoB(p_node_position, population_center);



                /*

                double a = endocells[endo_index][0] ;
                double b = endocells[endo_index][1] ;
                double l = endocells[endo_index][2] ;
                double r = endocells[endo_index][3] ;

                */

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


                double dist1 = sqrt(pow((x_coordinate - x_proj_l),2) + pow((y_coordinate - y_proj_l),2)) ; // SEE GetDistanceBetweenNodes
                double dist2 = sqrt(pow((x_coordinate - x_proj_r),2) + pow((y_coordinate - y_proj_r),2)) ;
                int innerproduct = inner_product(vec_from_cell_to_node.begin(), vec_from_cell_to_node.end(), vec_from_center_to_node.begin(), 0) ;
                //double dist3 = sqrt(pow((xcenter - xmoy),2) + pow((ycenter - ymoy),2)) ;
                //double dist4 = sqrt(pow((x_coordinate - xmoy),2) + pow((y_coordinate - ymoy),2)) ;
                // c_vector<double, SPACE_DIM> old_node_location = rOldLocations.find(p_node)->second; // node previous location

                if ( p_node->IsBoundaryNode() && (b*y_coordinate + a*x_coordinate) > -l && (b*y_coordinate + a*x_coordinate) < -r && dist1 < dist2 && innerproduct > 0)
                {
                  p_node->rGetModifiableLocation()[0] = x_proj_l;
                  p_node->rGetModifiableLocation()[1] = y_proj_l;
                //if (node_index == 324)  {std::cout << node_index << "1" << std::endl ; }
                  //p_node->rGetModifiableLocation() = old_node_location;
                }
                else if ( p_node->IsBoundaryNode() && (b*y_coordinate + a*x_coordinate) > -l && (b*y_coordinate + a*x_coordinate) < -r && dist1 > dist2 && innerproduct > 0 )
                {
                  p_node->rGetModifiableLocation()[0] = x_proj_r;
                  p_node->rGetModifiableLocation()[1] = y_proj_r ;
                //  if (node_index == 324)  {std::cout << node_index << "2"<< std::endl ;}
                  //p_node->rGetModifiableLocation() = old_node_location;
                }
                else if (p_node->IsBoundaryNode() && (b*y_coordinate + a*x_coordinate) < -l && (b*y_coordinate + a*x_coordinate) > -r && dist1 < dist2 && innerproduct > 0 )
                {
                  p_node->rGetModifiableLocation()[0] = x_proj_l;
                  p_node->rGetModifiableLocation()[1] = y_proj_l ;
                //  if (node_index == 324)  {std::cout << node_index<<"3" << std::endl ;}
                  //p_node->rGetModifiableLocation() = old_node_location;
                }
                else if (p_node->IsBoundaryNode() && (b*y_coordinate + a*x_coordinate) < -l && (b*y_coordinate + a*x_coordinate) > -r && dist1 > dist2 && innerproduct > 0 )
                {
                  p_node->rGetModifiableLocation()[0] = x_proj_r;
                  p_node->rGetModifiableLocation()[1] = y_proj_r ;
                //  if (node_index == 324)  {std::cout << node_index<<"4" << std::endl ;}
                  //p_node->rGetModifiableLocation() = old_node_location;
                }
                }

              }

    }
    else
    {
        // SPACE_DIM == 1
        std::cout << "no endo cells in periph" << std::endl ;
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
