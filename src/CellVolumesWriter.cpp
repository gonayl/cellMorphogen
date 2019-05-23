#include "CellPerimetersWriter.hpp"
#include "AbstractCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellPerimetersWriter<ELEMENT_DIM, SPACE_DIM>::CellPerimetersWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("cellperimeters.dat")
{
    this->mVtkCellDataName = "Cell perimeters";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellPerimetersWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double perimeter = pCellPopulation->GetperimeterOfCell(pCell);
    return perimeter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellPerimetersWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    unsigned location_index = pCellPopulation->GetLocationIndexUsingCell(pCell);
    unsigned cell_id = pCell->GetCellId();
    c_vector<double, SPACE_DIM> centre_location = pCellPopulation->GetLocationOfCellCentre(pCell);
    double perimeter = pCellPopulation->GetElementCorrespondingToCell(pCell)->GetSurfaceAreaOfElement();

    if (perimeter < DBL_MAX)   // Only write cells with finite perimeter (avoids a case for boundary cells in MeshBasedCellPopulation)
    {
        *this->mpOutStream << location_index << " " << cell_id << " ";
        for (unsigned i=0; i<SPACE_DIM; i++)
        {
            *this->mpOutStream << centre_location[i] << " ";
        }

        *this->mpOutStream << perimeter << " ";
    }
}

// Explicit instantiation
template class CellPerimetersWriter<1,1>;
template class CellPerimetersWriter<1,2>;
template class CellPerimetersWriter<2,2>;
template class CellPerimetersWriter<1,3>;
template class CellPerimetersWriter<2,3>;
template class CellPerimetersWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellPerimetersWriter)
