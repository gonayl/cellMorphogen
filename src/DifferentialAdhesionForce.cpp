#include "DifferentialAdhesionForce.hpp"
#include "CellEndo.hpp"
#include "CellLumen.hpp"
#include "CellEpi.hpp"

template<unsigned DIM>
DifferentialAdhesionForce<DIM>::DifferentialAdhesionForce()
    : NagaiHondaForce<DIM>(),
      mEndoEndoAdhesionEnergyParameter(1.0),
      mLumenLumenAdhesionEnergyParameter(1.0),
      mEpiEpiAdhesionEnergyParameter(1.0),
      mEndoEpiAdhesionEnergyParameter(1.0),
      mEpiLumenAdhesionEnergyParameter(1.0),
      mEndoLumenAdhesionEnergyParameter(1.0),
      mEndoBoundaryAdhesionEnergyParameter(1.0),
      mLumenBoundaryAdhesionEnergyParameter(1.0),
      mEpiBoundaryAdhesionEnergyParameter(1.0)
{
}

template<unsigned DIM>
DifferentialAdhesionForce<DIM>::~DifferentialAdhesionForce()
{
}

template<unsigned DIM>
double DifferentialAdhesionForce<DIM>::GetAdhesionParameter(Node<DIM>* pNodeA,
                                                                      Node<DIM>* pNodeB,
                                                                      VertexBasedCellPopulation<DIM>& rVertexCellPopulation)
{
    // Find the indices of the elements owned by each node
    std::set<unsigned> elements_containing_nodeA = pNodeA->rGetContainingElementIndices();
    std::set<unsigned> elements_containing_nodeB = pNodeB->rGetContainingElementIndices();

    // Find common elements
    std::set<unsigned> shared_elements;
    std::set_intersection(elements_containing_nodeA.begin(),
                          elements_containing_nodeA.end(),
                          elements_containing_nodeB.begin(),
                          elements_containing_nodeB.end(),
                          std::inserter(shared_elements, shared_elements.begin()));

    // Check that the nodes have a common edge
    assert(!shared_elements.empty());

    // If the edge corresponds to a single element, then the cell is on the boundary
    if (shared_elements.size() == 1)
    {
        unsigned element_index = *(shared_elements.begin());

        // Get cell associated with this element
        CellPtr p_cell = rVertexCellPopulation.GetCellUsingLocationIndex(element_index);

        if (p_cell->template HasCellProperty<CellEndo>())
        {
            // This cell is labelled "endo"
            return this->GetEndoBoundaryAdhesionEnergyParameter();
        }
        else if (p_cell->template HasCellProperty<CellLumen>())
        {
            // This cell is labelled "lumen"
            return this->GetLumenBoundaryAdhesionEnergyParameter();
        }
        else if (p_cell->template HasCellProperty<CellEpi>())
        {
            // This cell is labelled "lumen"
            return this->GetEpiBoundaryAdhesionEnergyParameter();
        }
        else
        {
            // This cell is not labelled
            EXCEPTION("All cells must be labelled to use DifferentialAdhesionForce ! ");
        }
    }
    else
    {
        // Work out the number of labelled cells: 0,1 or 2
        unsigned num_endo_cells = 0;
        unsigned num_lumen_cells = 0;
        unsigned num_epi_cells = 0;
        for (std::set<unsigned>::iterator iter = shared_elements.begin();
             iter != shared_elements.end();
             ++iter)
        {
            unsigned element_index = *(iter);

            // Get cell associated with this element
            CellPtr p_cell = rVertexCellPopulation.GetCellUsingLocationIndex(element_index);

            if (p_cell->template HasCellProperty<CellEndo>())
            {
                num_endo_cells++;
            }
            else if (p_cell->template HasCellProperty<CellLumen>())
            {
                num_lumen_cells++;
            }
            else if (p_cell->template HasCellProperty<CellEpi>())
            {
                num_epi_cells++;
            }
        }

        if (num_endo_cells == 2)
        {
            // Both cells are labelled "endo"
            return this->GetEndoEndoAdhesionEnergyParameter();
        }
        else if (num_lumen_cells == 2)
        {
            // Both cells are labelled "lumen"
            return this->GetLumenLumenAdhesionEnergyParameter();
        }
        else if (num_epi_cells == 2)
        {
            // Both cells are labelled "lumen"
            return this->GetEpiEpiAdhesionEnergyParameter();
        }
        else if (num_endo_cells == 1 && num_lumen_cells == 1)
        {
            // One cell is labelled "endo" and the other "lumen"
            return this->GetEndoLumenAdhesionEnergyParameter();
        }
        else if (num_endo_cells == 1 && num_epi_cells == 1)
        {
            // One cell is labelled "endo" and the other "epi""
            return this->GetEndoEpiAdhesionEnergyParameter();
        }
        else if (num_lumen_cells == 1 && num_epi_cells == 1)
        {
            // One cell is labelled "endo" and the other "epi""
            return this->GetEpiLumenAdhesionEnergyParameter();
        }
        else
        {
          // At least one cell is not labelled
          EXCEPTION("All cells must be labelled to use DifferentialAdhesionForce ! ");
        }
    }
}

template<unsigned DIM>
double DifferentialAdhesionForce<DIM>::GetEndoEndoAdhesionEnergyParameter()
{
    return mEndoEndoAdhesionEnergyParameter;
}

template<unsigned DIM>
double DifferentialAdhesionForce<DIM>::GetLumenLumenAdhesionEnergyParameter()
{
    return mLumenLumenAdhesionEnergyParameter;
}

template<unsigned DIM>
double DifferentialAdhesionForce<DIM>::GetEpiEpiAdhesionEnergyParameter()
{
    return mEpiEpiAdhesionEnergyParameter;
}

template<unsigned DIM>
double DifferentialAdhesionForce<DIM>::GetEndoEpiAdhesionEnergyParameter()
{
    return mEndoEpiAdhesionEnergyParameter;
}

template<unsigned DIM>
double DifferentialAdhesionForce<DIM>::GetEpiLumenAdhesionEnergyParameter()
{
    return mEpiLumenAdhesionEnergyParameter;
}

template<unsigned DIM>
double DifferentialAdhesionForce<DIM>::GetEndoLumenAdhesionEnergyParameter()
{
    return mEndoLumenAdhesionEnergyParameter;
}

template<unsigned DIM>
double DifferentialAdhesionForce<DIM>::GetEndoBoundaryAdhesionEnergyParameter()
{
    return mEndoBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
double DifferentialAdhesionForce<DIM>::GetLumenBoundaryAdhesionEnergyParameter()
{
    return mLumenBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
double DifferentialAdhesionForce<DIM>::GetEpiBoundaryAdhesionEnergyParameter()
{
    return mEpiBoundaryAdhesionEnergyParameter;
}



/*-------*/

template<unsigned DIM>
void DifferentialAdhesionForce<DIM>::SetEndoEndoAdhesionEnergyParameter(double endoEndoAdhesionEnergyParameter)
{
    mEndoEndoAdhesionEnergyParameter = endoEndoAdhesionEnergyParameter;
}

template<unsigned DIM>
void DifferentialAdhesionForce<DIM>::SetLumenLumenAdhesionEnergyParameter(double lumenLumenAdhesionEnergyParameter)
{
    mLumenLumenAdhesionEnergyParameter = lumenLumenAdhesionEnergyParameter;
}

template<unsigned DIM>
void DifferentialAdhesionForce<DIM>::SetEpiEpiAdhesionEnergyParameter(double epiEpiAdhesionEnergyParameter)
{
    mEpiEpiAdhesionEnergyParameter = epiEpiAdhesionEnergyParameter;
}

template<unsigned DIM>
void DifferentialAdhesionForce<DIM>::SetEndoEpiAdhesionEnergyParameter(double endoEpiAdhesionEnergyParameter)
{
    mEndoEpiAdhesionEnergyParameter = endoEpiAdhesionEnergyParameter;
}

template<unsigned DIM>
void DifferentialAdhesionForce<DIM>::SetEpiLumenAdhesionEnergyParameter(double epiLumenAdhesionEnergyParameter)
{
    mEpiLumenAdhesionEnergyParameter = epiLumenAdhesionEnergyParameter;
}

template<unsigned DIM>
void DifferentialAdhesionForce<DIM>::SetEndoLumenAdhesionEnergyParameter(double endoLumenAdhesionEnergyParameter)
{
    mEndoLumenAdhesionEnergyParameter = endoLumenAdhesionEnergyParameter;
}

template<unsigned DIM>
void DifferentialAdhesionForce<DIM>::SetEndoBoundaryAdhesionEnergyParameter(double endoBoundaryAdhesionEnergyParameter)
{
    mEndoBoundaryAdhesionEnergyParameter = endoBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
void DifferentialAdhesionForce<DIM>::SetLumenBoundaryAdhesionEnergyParameter(double lumenBoundaryAdhesionEnergyParameter)
{
    mLumenBoundaryAdhesionEnergyParameter = lumenBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
void DifferentialAdhesionForce<DIM>::SetEpiBoundaryAdhesionEnergyParameter(double epiBoundaryAdhesionEnergyParameter)
{
    mEpiBoundaryAdhesionEnergyParameter = epiBoundaryAdhesionEnergyParameter;
}

template<unsigned DIM>
void DifferentialAdhesionForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // Output member variables
    *rParamsFile << "\t\t\t<EndoEndoAdhesionEnergyParameter>" << mEndoEndoAdhesionEnergyParameter << "</EndoEndoAdhesionEnergyParameter> \n";

    *rParamsFile << "\t\t\t<LumenLumenAdhesionEnergyParameter>" << mLumenLumenAdhesionEnergyParameter << "</LumenLumenAdhesionEnergyParameter> \n";

    *rParamsFile << "\t\t\t<EpiEpiAdhesionEnergyParameter>" << mEpiEpiAdhesionEnergyParameter << "</EpiEpiAdhesionEnergyParameter> \n";

    *rParamsFile << "\t\t\t<EndoEpiAdhesionEnergyParameter>" << mEndoEpiAdhesionEnergyParameter << "</EndoEpiAdhesionEnergyParameter> \n";

    *rParamsFile << "\t\t\t<EpiLumenAdhesionEnergyParameter>" << mEpiLumenAdhesionEnergyParameter << "</EpiLumenAdhesionEnergyParameter> \n";

    *rParamsFile << "\t\t\t<EndoLumenAdhesionEnergyParameter>" << mEndoLumenAdhesionEnergyParameter << "</EndoLumenAdhesionEnergyParameter> \n";

    *rParamsFile << "\t\t\t<EndoBoundaryAdhesionEnergyParameter>" << mEndoBoundaryAdhesionEnergyParameter << "</EndoBoundaryAdhesionEnergyParameter> \n";

    *rParamsFile << "\t\t\t<LumenBoundaryAdhesionEnergyParameter>" << mLumenBoundaryAdhesionEnergyParameter << "</LumenBoundaryAdhesionEnergyParameter> \n";

    *rParamsFile << "\t\t\t<EpiBoundaryAdhesionEnergyParameter>" << mEpiBoundaryAdhesionEnergyParameter << "</EpiBoundaryAdhesionEnergyParameter> \n";

    // Call method on direct parent class
    NagaiHondaForce<DIM>::OutputForceParameters(rParamsFile);
}

// Explicit instantiation
template class DifferentialAdhesionForce<1>;
template class DifferentialAdhesionForce<2>;
template class DifferentialAdhesionForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DifferentialAdhesionForce)
