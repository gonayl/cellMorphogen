#include "CellElongationWriter.hpp"
#include "AbstractCellPopulation.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellElongationWriter<ELEMENT_DIM, SPACE_DIM>::CellElongationWriter()
    : AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>("cellareas.dat")
{
    this->mVtkCellDataName = "Cell elongation";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
double CellElongationWriter<ELEMENT_DIM, SPACE_DIM>::GetCellDataForVtkOutput(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    double volume = pCellPopulation->GetVolumeOfCell(pCell);
    return volume;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CellElongationWriter<ELEMENT_DIM, SPACE_DIM>::VisitCell(CellPtr pCell, AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    unsigned location_index = pCellPopulation->GetLocationIndexUsingCell(pCell);
    unsigned cell_id = pCell->GetCellId();
    c_vector<double, SPACE_DIM> centre_location = pCellPopulation->GetLocationOfCellCentre(pCell);
    double volume = pCellPopulation->GetVolumeOfCell(pCell);

    if (volume < DBL_MAX)   // Only write cells with finite volume (avoids a case for boundary cells in MeshBasedCellPopulation)
    {
        *this->mpOutStream << location_index << " " << cell_id << " ";
        for (unsigned i=0; i<SPACE_DIM; i++)
        {
            *this->mpOutStream << centre_location[i] << " ";
        }

        *this->mpOutStream << volume << " ";
    }
}

// Explicit instantiation
template class CellElongationWriter<1,1>;
template class CellElongationWriter<1,2>;
template class CellElongationWriter<2,2>;
template class CellElongationWriter<1,3>;
template class CellElongationWriter<2,3>;
template class CellElongationWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellElongationWriter)
