#include "DifferentialAdhesionForce.hpp"
#include "CellEndo.hpp"
#include "CellLumen.hpp"

template<unsigned DIM>
DifferentialAdhesionForce<DIM>::DifferentialAdhesionForce()
    : NagaiHondaForce<DIM>(),
      mEndoEndoAdhesionEnergyParameter(1.0),
      mLumenLumenAdhesionEnergyParameter(1.0),
      mEndoEpiAdhesionEnergyParameter(1.0),
      mLumenEpiAdhesionEnergyParameter(1.0),
      mLumenEndoAdhesionEnergyParameter(1.0),
      mEndoBoundaryAdhesionEnergyParameter(1.0),
      mLumenBoundaryAdhesionEnergyParameter(1.0)
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
        else
        {
            // This cell is not labelled
            return this->GetNagaiHondaCellBoundaryAdhesionEnergyParameter();
        }
    }
    else
    {
        // Work out the number of labelled cells: 0,1 or 2
        unsigned num_endo_cells = 0;
        unsigned num_lumen_cells = 0;
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
        else if (num_endo_cells == 1 && num_lumen_cells == 1)
        {
            // One cell is labelled "endo" and the other "lumen"
            return this->GetLumenEndoAdhesionEnergyParameter();
        }
        else if (num_endo_cells == 1 && num_lumen_cells == 0)
        {
            // One cell is labelled "endo" and the other is not labelled
            return this->GetEndoEpiAdhesionEnergyParameter();
        }
        else if (num_lumen_cells == 1 && num_endo_cells == 0)
        {
            // One cell is labelled "lumen"
            return this->GetLumenEpiAdhesionEnergyParameter();
        }
        else
        {
            // Neither cell is labelled
            assert(num_endo_cells == 0);
            return this->GetNagaiHondaCellCellAdhesionEnergyParameter();
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
double DifferentialAdhesionForce<DIM>::GetEndoEpiAdhesionEnergyParameter()
{
    return mEndoEpiAdhesionEnergyParameter;
}

template<unsigned DIM>
double DifferentialAdhesionForce<DIM>::GetLumenEpiAdhesionEnergyParameter()
{
    return mLumenEpiAdhesionEnergyParameter;
}

template<unsigned DIM>
double DifferentialAdhesionForce<DIM>::GetLumenEndoAdhesionEnergyParameter()
{
    return mLumenEndoAdhesionEnergyParameter;
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
void DifferentialAdhesionForce<DIM>::SetEndoEpiAdhesionEnergyParameter(double endoEpiAdhesionEnergyParameter)
{
    mEndoEpiAdhesionEnergyParameter = endoEpiAdhesionEnergyParameter;
}

template<unsigned DIM>
void DifferentialAdhesionForce<DIM>::SetLumenEpiAdhesionEnergyParameter(double lumenEpiAdhesionEnergyParameter)
{
    mLumenEpiAdhesionEnergyParameter = lumenEpiAdhesionEnergyParameter;
}

template<unsigned DIM>
void DifferentialAdhesionForce<DIM>::SetLumenEndoAdhesionEnergyParameter(double lumenEndoAdhesionEnergyParameter)
{
    mLumenEndoAdhesionEnergyParameter = lumenEndoAdhesionEnergyParameter;
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
void DifferentialAdhesionForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    // Output member variables
    *rParamsFile << "\t\t\t<EndoEndoAdhesionEnergyParameter>" << mEndoEndoAdhesionEnergyParameter << "</EndoEndoAdhesionEnergyParameter> \n";

    *rParamsFile << "\t\t\t<LumenLumenAdhesionEnergyParameter>" << mLumenLumenAdhesionEnergyParameter << "</LumenLumenAdhesionEnergyParameter> \n";

    *rParamsFile << "\t\t\t<EndoEpiAdhesionEnergyParameter>" << mEndoEpiAdhesionEnergyParameter << "</EndoEpiAdhesionEnergyParameter> \n";

    *rParamsFile << "\t\t\t<LumenEpiAdhesionEnergyParameter>" << mLumenEpiAdhesionEnergyParameter << "</LumenEpiAdhesionEnergyParameter> \n";

    *rParamsFile << "\t\t\t<LumenEndoAdhesionEnergyParameter>" << mLumenEndoAdhesionEnergyParameter << "</LumenEndoAdhesionEnergyParameter> \n";

    *rParamsFile << "\t\t\t<EndoBoundaryAdhesionEnergyParameter>" << mEndoBoundaryAdhesionEnergyParameter << "</EndoBoundaryAdhesionEnergyParameter> \n";

    *rParamsFile << "\t\t\t<LumenBoundaryAdhesionEnergyParameter>" << mLumenBoundaryAdhesionEnergyParameter << "</LumenBoundaryAdhesionEnergyParameter> \n";

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
