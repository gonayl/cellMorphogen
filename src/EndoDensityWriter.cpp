

#include "EndoDensityWriter.hpp"
#include "AbstractCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "ApoptoticCellProperty.hpp"

#include "CellLumen.hpp"
#include "CellEndo.hpp"
#include "CellEpi.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
EndoDensityWriter<ELEMENT_DIM, SPACE_DIM>::EndoDensityWriter()
    : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>("EndoDensity.dat")
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void EndoDensityWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("EndoDensityWriter cannot be used with a MeshBasedCellPopulation");

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void EndoDensityWriter<ELEMENT_DIM, SPACE_DIM>::Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("EndoDensityWriter cannot be used with a CaBasedCellPopulation");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void EndoDensityWriter<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("EndoDensityWriter cannot be used with a NodeBasedCellPopulation");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void EndoDensityWriter<ELEMENT_DIM, SPACE_DIM>::Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("EndoDensityWriter cannot be used with a PottsBasedCellPopulation");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void EndoDensityWriter<ELEMENT_DIM, SPACE_DIM>::Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{

  double total_area = 0.0;
  double endo_area = 0.0;
  double lumen_area = 0.0;
  double epi_area = 0.0;

  for (typename AbstractCellPopulation<SPACE_DIM>::Iterator cell_iter = pCellPopulation->Begin();
       cell_iter != pCellPopulation->End();
       ++cell_iter)
  {
    CellPtr pCell = *cell_iter;

    //VertexElement<SPACE_DIM>* p_element = p_cell_population->GetElementCorrespondingToCell(*cell_iter);
    unsigned p_index = pCellPopulation->GetLocationIndexUsingCell(pCell);
    double cell_size = pCellPopulation->rGetMesh().GetVolumeOfElement(p_index);



    if (pCell->HasCellProperty<CellEndo>())
    {
      endo_area = endo_area + cell_size;
    }
    else if (pCell->HasCellProperty<CellLumen>())
    {
      lumen_area = lumen_area + cell_size;
    }
    else if (pCell->HasCellProperty<CellEpi>())
    {
      epi_area = epi_area + cell_size;
    }
    total_area = total_area + cell_size;
  }
  *this->mpOutStream << total_area << " " << endo_area << " " << lumen_area << " " << epi_area;
}

// Explicit instantiation
template class EndoDensityWriter<1,1>;
template class EndoDensityWriter<1,2>;
template class EndoDensityWriter<2,2>;
template class EndoDensityWriter<1,3>;
template class EndoDensityWriter<2,3>;
template class EndoDensityWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(EndoDensityWriter)
