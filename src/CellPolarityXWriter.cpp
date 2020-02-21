#include "CellPolarityXWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "CellEpi.hpp"


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellPolarityXWriter<ELEMENT_DIM, SPACE_DIM>::CellPolarityXWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("results.viztypes")
{
    this->mVtkCellDataName = "Cellular polarisastion X";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellPolarityXWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double vecPolaX = 0;
    if (pCell->HasCellProperty<CellEpi>())
    {
          vecPolaX = pCell->GetCellData()->GetItem("vecPolaX");
          //std::cout << vecPolaX << std::endl;
    }


    return vecPolaX;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPolarityXWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double vecPolaX = 0;
    if (pCell->HasCellProperty<CellEpi>())
    {
        //CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<CellEpi>();
        //boost::shared_ptr<CellEpi> p_celltype = boost::static_pointer_cast<CellEpi>(collection.GetProperty());
        vecPolaX = pCell->GetCellData()->GetItem("vecPolaX");
    }


    *this->mpOutStream << " " << vecPolaX;

    unsigned location_index = pCellPopulation->GetLocationIndexUsingCell(pCell);
    *this->mpOutStream << " " << location_index;

    c_vector<double, SPACE_DIM> coords = pCellPopulation->GetLocationOfCellCentre(pCell);
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        *this->mpOutStream << " " << coords[i];
    }
}

// Explicit instantiation
template class CellPolarityXWriter<1,1>;
template class CellPolarityXWriter<1,2>;
template class CellPolarityXWriter<2,2>;
template class CellPolarityXWriter<1,3>;
template class CellPolarityXWriter<2,3>;
template class CellPolarityXWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellPolarityXWriter)
