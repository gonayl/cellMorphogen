#include "CellPolarityWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "CellEpi.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellPolarityWriter<ELEMENT_DIM, SPACE_DIM>::CellPolarityWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("results.viztypes")
{
    this->mVtkCellDataName = "Cellular polarisastion X";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellPolarityWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double vecPolaX = 0;
    if (pCell->HasCellProperty<CellEpi>())
    {
          vecPolaX = pCell->GetCellData()->GetItem("vecPolaX");
          vecPolaY = pCell->GetCellData()->GetItem("vecPolaY");
          //std::cout << vecPolaX << std::endl;
    }


    return vecPolaX;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPolarityWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double vecPolaX = 0;
    if (pCell->HasCellProperty<CellEpi>())
    {
        //CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<CellEpi>();
        //boost::shared_ptr<CellEpi> p_celltype = boost::static_pointer_cast<CellEpi>(collection.GetProperty());
        vecPolaX = pCell->GetCellData()->GetItem("vecPolaX");
        vecPolaY = pCell->GetCellData()->GetItem("vecPolaY");
    }

    unsigned location_index = pCellPopulation->GetLocationIndexUsingCell(pCell);
    *this->mpOutStream << " " << location_index;

    *this->mpOutStream << " " << vecPolaX << " " << vecPolaY ;

}

// Explicit instantiation
template class CellPolarityWriter<1,1>;
template class CellPolarityWriter<1,2>;
template class CellPolarityWriter<2,2>;
template class CellPolarityWriter<1,3>;
template class CellPolarityWriter<2,3>;
template class CellPolarityWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellPolarityWriter)
