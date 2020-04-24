#include "CellAllTypeWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "CellEndo.hpp"
#include "CellLumen.hpp"
#include "CellPeriph.hpp"
#include "CellCore.hpp"
#include "CellStalk.hpp"
#include "CellTip.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellAllTypeWriter<ELEMENT_DIM, SPACE_DIM>::CellAllTypeWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("results.viztypes")
{
    this->mVtkCellDataName = "Cellular All type";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellAllTypeWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double cell_type = 0.0;
    if (pCell->HasCellProperty<CellLumen>())
    {
        CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<CellLumen>();
        boost::shared_ptr<CellLumen> p_celltype = boost::static_pointer_cast<CellLumen>(collection.GetProperty());
        cell_type = p_celltype->GetColour();
    }
    else if (pCell->HasCellProperty<CellCore>())
    {
        CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<CellCore>();
        boost::shared_ptr<CellCore> p_celltype = boost::static_pointer_cast<CellCore>(collection.GetProperty());
        cell_type = p_celltype->GetColour();
    }
    else if (pCell->HasCellProperty<CellPeriph>())
    {
        CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<CellPeriph>();
        boost::shared_ptr<CellPeriph> p_celltype = boost::static_pointer_cast<CellPeriph>(collection.GetProperty());
        cell_type = p_celltype->GetColour();
    }
    if (pCell->HasCellProperty<CellTip>())
    {
        CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<CellTip>();
        boost::shared_ptr<CellTip> p_celltype = boost::static_pointer_cast<CellTip>(collection.GetProperty());
        cell_type = p_celltype->GetColour();
    }
    else if (pCell->HasCellProperty<CellEndo>())
    {
        CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<CellEndo>();
        boost::shared_ptr<CellEndo> p_celltype = boost::static_pointer_cast<CellEndo>(collection.GetProperty());
        cell_type = p_celltype->GetColour();
    }

    return cell_type;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellAllTypeWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    unsigned cell_type = 0;
if (pCell->HasCellProperty<CellLumen>())
    {
        CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<CellLumen>();
        boost::shared_ptr<CellLumen> p_celltype = boost::static_pointer_cast<CellLumen>(collection.GetProperty());
        cell_type = p_celltype->GetColour();
    }
    else if (pCell->HasCellProperty<CellCore>())
    {
        CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<CellCore>();
        boost::shared_ptr<CellCore> p_celltype = boost::static_pointer_cast<CellCore>(collection.GetProperty());
        cell_type = p_celltype->GetColour();
    }
    else if (pCell->HasCellProperty<CellPeriph>())
    {
        CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<CellPeriph>();
        boost::shared_ptr<CellPeriph> p_celltype = boost::static_pointer_cast<CellPeriph>(collection.GetProperty());
        cell_type = p_celltype->GetColour();
    }
    if (pCell->HasCellProperty<CellTip>())
    {
        CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<CellTip>();
        boost::shared_ptr<CellTip> p_celltype = boost::static_pointer_cast<CellTip>(collection.GetProperty());
        cell_type = p_celltype->GetColour();
    }
    else if (pCell->HasCellProperty<CellEndo>())
    {
        CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<CellEndo>();
        boost::shared_ptr<CellEndo> p_celltype = boost::static_pointer_cast<CellEndo>(collection.GetProperty());
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
template class CellAllTypeWriter<1,1>;
template class CellAllTypeWriter<1,2>;
template class CellAllTypeWriter<2,2>;
template class CellAllTypeWriter<1,3>;
template class CellAllTypeWriter<2,3>;
template class CellAllTypeWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellAllTypeWriter)
