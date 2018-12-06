#include "CellTypeWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "CellEndo.hpp"
#include "CellLumen.hpp"
#include "CellEpi.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellTypeWriter<ELEMENT_DIM, SPACE_DIM>::CellTypeWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("results.viztypes")
{
    this->mVtkCellDataName = "Cellular type";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellTypeWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double cell_type = 0.0;
    if (pCell->HasCellProperty<CellEndo>())
    {
        CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<CellEndo>();
        boost::shared_ptr<CellEndo> p_celltype = boost::static_pointer_cast<CellEndo>(collection.GetProperty());
        cell_type = p_celltype->GetColour();
    }
    else if (pCell->HasCellProperty<CellLumen>())
    {
        CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<CellLumen>();
        boost::shared_ptr<CellLumen> p_celltype = boost::static_pointer_cast<CellLumen>(collection.GetProperty());
        cell_type = p_celltype->GetColour();
    }
    else if (pCell->HasCellProperty<CellEpi>())
    {
        CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<CellEpi>();
        boost::shared_ptr<CellEpi> p_celltype = boost::static_pointer_cast<CellEpi>(collection.GetProperty());
        cell_type = p_celltype->GetColour();
    }
    return cell_type;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellTypeWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    unsigned cell_type = 0;
    if (pCell->HasCellProperty<CellEndo>())
    {
        CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<CellEndo>();
        boost::shared_ptr<CellEndo> p_celltype = boost::static_pointer_cast<CellEndo>(collection.GetProperty());
        cell_type = p_celltype->GetColour();
    }
    else if (pCell->HasCellProperty<CellLumen>())
    {
        CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<CellLumen>();
        boost::shared_ptr<CellLumen> p_celltype = boost::static_pointer_cast<CellLumen>(collection.GetProperty());
        cell_type = p_celltype->GetColour();
    }
    else if (pCell->HasCellProperty<CellEpi>())
    {
        CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<CellEpi>();
        boost::shared_ptr<CellEpi> p_celltype = boost::static_pointer_cast<CellEpi>(collection.GetProperty());
        cell_type = p_celltype->GetColour();
    }

    *this->mpOutStream << " " << cell_type;

    unsigned location_index = pCellPopulation->GetLocationIndexUsingCell(pCell);
    *this->mpOutStream << " " << location_index;

    c_vector<double, SPACE_DIM> coords = pCellPopulation->GetLocationOfCellCentre(pCell);
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        *this->mpOutStream << " " << coords[i];
    }
}

// Explicit instantiation
template class CellTypeWriter<1,1>;
template class CellTypeWriter<1,2>;
template class CellTypeWriter<2,2>;
template class CellTypeWriter<1,3>;
template class CellTypeWriter<2,3>;
template class CellTypeWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellTypeWriter)
