#include "CellPosWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "CellPeriph.hpp"
#include "CellCore.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellPosWriter<ELEMENT_DIM, SPACE_DIM>::CellPosWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("results.viztypes")
{
    this->mVtkCellDataName = "Cell position";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellPosWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double cell_pos = 0.0;
    if (pCell->HasCellProperty<CellPeriph>())
    {
        CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<CellPeriph>();
        boost::shared_ptr<CellPeriph> p_celltype = boost::static_pointer_cast<CellPeriph>(collection.GetProperty());
        cell_pos = p_celltype->GetColour();
    }
    else if (pCell->HasCellProperty<CellCore>())
    {
        CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<CellCore>();
        boost::shared_ptr<CellCore> p_celltype = boost::static_pointer_cast<CellCore>(collection.GetProperty());
        cell_pos = p_celltype->GetColour();
    }

    return cell_pos;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPosWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    unsigned cell_pos = 0.0;
    if (pCell->HasCellProperty<CellPeriph>())
    {
        CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<CellPeriph>();
        boost::shared_ptr<CellPeriph> p_celltype = boost::static_pointer_cast<CellPeriph>(collection.GetProperty());
        cell_pos = p_celltype->GetColour();
    }
    else if (pCell->HasCellProperty<CellCore>())
    {
        CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<CellCore>();
        boost::shared_ptr<CellCore> p_celltype = boost::static_pointer_cast<CellCore>(collection.GetProperty());
        cell_pos = p_celltype->GetColour();
    }


    *this->mpOutStream << " " << cell_pos;

    unsigned location_index = pCellPopulation->GetLocationIndexUsingCell(pCell);
    *this->mpOutStream << " " << location_index;

    c_vector<double, SPACE_DIM> coords = pCellPopulation->GetLocationOfCellCentre(pCell);
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        *this->mpOutStream << " " << coords[i];
    }
}

// Explicit instantiation
template class CellPosWriter<1,1>;
template class CellPosWriter<1,2>;
template class CellPosWriter<2,2>;
template class CellPosWriter<1,3>;
template class CellPosWriter<2,3>;
template class CellPosWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellPosWriter)
