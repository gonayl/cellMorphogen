#include "NeighbourTrackingModifier.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "CellEndo.hpp"
#include "CellEpi.hpp"
#include "CellLumen.hpp"
#include "CellPolar.hpp"


template<unsigned DIM>
NeighbourTrackingModifier<DIM>::NeighbourTrackingModifier()
    : AbstractCellBasedSimulationModifier<DIM>()
{
}

template<unsigned DIM>
NeighbourTrackingModifier<DIM>::~NeighbourTrackingModifier()
{
}

template<unsigned DIM>
void NeighbourTrackingModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void NeighbourTrackingModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    /*
     * We must update CellData in SetupSolve(), otherwise it will not have been
     * fully initialised by the time we enter the main time loop.
     */
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void NeighbourTrackingModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
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
        EXCEPTION("Neighbour Tracking Modifier is to be used with a VertexBasedCellPopulation only");
    }

    VertexBasedCellPopulation<DIM>* p_cell_population = static_cast<VertexBasedCellPopulation<DIM>*>(&rCellPopulation);

    // Iterate over cell population

    unsigned num_cells = p_cell_population->GetNumRealCells();

    // Store a map between cells numbered 1 to n and location indices
    std::map<unsigned,unsigned> local_cell_id_location_index_map;

    unsigned local_cell_id = 0;
    for (typename AbstractCellPopulation<DIM, DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        local_cell_id_location_index_map[p_cell_population->GetLocationIndexUsingCell(*cell_iter)] = local_cell_id;
        local_cell_id++;
    }
    assert(local_cell_id = num_cells+1);

    // Iterate over cells and calculate the number of endothelial neighbouring cells
    // bool cell_is_epi ;

    for (typename AbstractCellPopulation<DIM, DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {

        double n_neighbour_endo = 0 ;
        cell_iter->template RemoveCellProperty<CellPolar>();

        // Determine whether this cell is labelled
        // bool cell_is_epi = cell_iter->template HasCellProperty<CellEpi>();

        // Get the set of neighbouring location indices
        std::set<unsigned> neighbour_indices = p_cell_population->GetNeighbouringLocationIndices(*cell_iter);

        // If this cell has any neighbours (as defined by mesh/population/interaction distance)...
        if (!neighbour_indices.empty() )
        {

            for (std::set<unsigned>::iterator neighbour_iter = neighbour_indices.begin();
                 neighbour_iter != neighbour_indices.end();
                 ++neighbour_iter)
            {
                // Determine whether this neighbour is endo or "polar"
                CellPtr p_neighbour_cell = p_cell_population->GetCellUsingLocationIndex(*neighbour_iter);
                bool neighbour_is_endo = p_neighbour_cell->template HasCellProperty<CellEndo>();

                if ( neighbour_is_endo == 1)
                {
                    n_neighbour_endo++;
                }
            }

            cell_iter->GetCellData()->SetItem("nendoneighbours", n_neighbour_endo);

            if(n_neighbour_endo != 0)
            {
              cell_iter->AddCellProperty(rCellPopulation.GetCellPropertyRegistry()->template Get<CellPolar>());
            }
        }
    }
}

template<unsigned DIM>
void NeighbourTrackingModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class NeighbourTrackingModifier<1>;
template class NeighbourTrackingModifier<2>;
template class NeighbourTrackingModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(NeighbourTrackingModifier)
