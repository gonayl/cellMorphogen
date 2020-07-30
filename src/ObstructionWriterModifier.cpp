#include "ObstructionWriterModifier.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CounterSingletonRepulsion.hpp"
#include "CellEndo.hpp"
#include "CellEpi.hpp"
#include "CellLumen.hpp"
#include "CellPolar.hpp"
#include "CellPeriph.hpp"
#include "CellCore.hpp"
#include "CellStalk.hpp"
#include <stdlib.h>
using namespace std ;


template<unsigned SPACE_DIM>
ObstructionWriterModifier<SPACE_DIM>::ObstructionWriterModifier()
    : AbstractCellBasedSimulationModifier<SPACE_DIM>()
{
}

template<unsigned SPACE_DIM>
ObstructionWriterModifier<SPACE_DIM>::~ObstructionWriterModifier()
{
}

template<unsigned SPACE_DIM>
void ObstructionWriterModifier<SPACE_DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<SPACE_DIM,SPACE_DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned SPACE_DIM>
void ObstructionWriterModifier<SPACE_DIM>::SetupSolve(AbstractCellPopulation<SPACE_DIM,SPACE_DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned SPACE_DIM>
void ObstructionWriterModifier<SPACE_DIM>::UpdateCellData(AbstractCellPopulation<SPACE_DIM,SPACE_DIM>& rCellPopulation)
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
    if (bool(dynamic_cast<MeshBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation)))
    {
        static_cast<MeshBasedCellPopulation<SPACE_DIM>*>(&(rCellPopulation))->CreateVoronoiTessellation();
    }

    if (dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation) == NULL)
    {
        EXCEPTION("Border Tracking Modifier is to be used with a VertexBasedCellPopulation only");
    }

    c_vector<double, SPACE_DIM> population_center = rCellPopulation.GetCentroidOfCellPopulation();


    //const int n_columns = 4;
    //double endocells[][n_columns]{} ;
    vector<double> endo_cells_tofile;
    int endo_count = 0 ;
    int endo_count_local = 0 ;

    int endo_count_global = CounterSingletonRepulsion::Instance()->GetCount() ;

    //cout << endo_count_global << "endo stalk fixed cell(s) "  << endl ;


    for (typename AbstractCellPopulation<SPACE_DIM,SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)

    {
      if (cell_iter->template HasCellProperty<CellStalk>() && cell_iter->template HasCellProperty<CellPeriph>())
      {

        endo_count_local++ ;

      }

    }

    if (endo_count_global != endo_count_local)
    {

      CounterSingletonRepulsion::Instance()->ResetCountToZero() ;

      // PHASE 0, computing mass center (will not be used after)

    //  vector<double> moy_x_pos;
    //  vector<double> moy_y_pos;

    ofstream myfile ("example.txt", std::ofstream::out | std::ofstream::trunc);

    for (typename AbstractCellPopulation<SPACE_DIM,SPACE_DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)

    {

      if ( cell_iter->template HasCellProperty<CellStalk>() && cell_iter->template HasCellProperty<CellPeriph>())
      {


        vector<double> boundary_nodes_pos ;
        VertexBasedCellPopulation<SPACE_DIM>* p_cell_population = dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation) ;
        VertexElement<SPACE_DIM,SPACE_DIM>* p_element ;

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

        direct_vector[0] = centre_of_cell[0] - population_center[0] ;
        direct_vector[1] = centre_of_cell[1] - population_center[1] ;

        a =  direct_vector[1] ; // a
        b =  - direct_vector[0] ; // b
        l = - direct_vector[1]*boundary_nodes_pos[0] + direct_vector[0]*boundary_nodes_pos[1] ; // l
        r = - direct_vector[1]*boundary_nodes_pos[2] + direct_vector[0]*boundary_nodes_pos[3] ; // r
        // WARNING : the order is important !!

        endo_cells_tofile.push_back(a) ;
        endo_cells_tofile.push_back(b) ;
        endo_cells_tofile.push_back(l) ;
        endo_cells_tofile.push_back(r) ;

        endo_cells_tofile.push_back(centre_of_cell[0]) ;
        endo_cells_tofile.push_back(centre_of_cell[1]) ;

        cout << "pushing data for cell " << rCellPopulation.GetLocationIndexUsingCell(*cell_iter) << endl ;

        endo_count++ ;
        CounterSingletonRepulsion::Instance()->IncrementCounter();

        if (myfile.is_open())
        {
          myfile << a << " " << b << " " << l << " " << r << " " << centre_of_cell[0] << " " << centre_of_cell[1] << " " ;
        }

      }

    }

    }

}

template<unsigned SPACE_DIM>
void ObstructionWriterModifier<SPACE_DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<SPACE_DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class ObstructionWriterModifier<1>;
template class ObstructionWriterModifier<2>;
template class ObstructionWriterModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ObstructionWriterModifier)
