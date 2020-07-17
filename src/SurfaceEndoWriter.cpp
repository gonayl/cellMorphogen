

#include "SurfaceEndoWriter.hpp"
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
SurfaceEndoWriter<ELEMENT_DIM, SPACE_DIM>::SurfaceEndoWriter()
    : AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>("SurfaceEndo.dat")
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void SurfaceEndoWriter<ELEMENT_DIM, SPACE_DIM>::Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("SurfaceEndoWriter cannot be used with a MeshBasedCellPopulation");

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void SurfaceEndoWriter<ELEMENT_DIM, SPACE_DIM>::Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("SurfaceEndoWriter cannot be used with a CaBasedCellPopulation");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void SurfaceEndoWriter<ELEMENT_DIM, SPACE_DIM>::Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("SurfaceEndoWriter cannot be used with a NodeBasedCellPopulation");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void SurfaceEndoWriter<ELEMENT_DIM, SPACE_DIM>::Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{
    EXCEPTION("SurfaceEndoWriter cannot be used with a PottsBasedCellPopulation");
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void SurfaceEndoWriter<ELEMENT_DIM, SPACE_DIM>::Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation)
{


  //std::string surface_tot = "";

  for (typename AbstractCellPopulation<SPACE_DIM>::Iterator cell_iter = pCellPopulation->Begin();
       cell_iter != pCellPopulation->End();
       ++cell_iter)
  {
    CellPtr pCell = *cell_iter;

    //VertexElement<SPACE_DIM>* p_element = p_cell_population->GetElementCorrespondingToCell(*cell_iter);




    if (pCell->HasCellProperty<CellEndo>())
    {
      unsigned p_index = pCellPopulation->GetLocationIndexUsingCell(pCell);
      double cell_size = pCellPopulation->rGetMesh().GetVolumeOfElement(p_index);
      //surface_tot = surface_tot << " " << cell_size;
      *this->mpOutStream << cell_size<< " ";
    }

  }
  //*this->mpOutStream << surface_tot;
}

// Explicit instantiation
template class SurfaceEndoWriter<1,1>;
template class SurfaceEndoWriter<1,2>;
template class SurfaceEndoWriter<2,2>;
template class SurfaceEndoWriter<1,3>;
template class SurfaceEndoWriter<2,3>;
template class SurfaceEndoWriter<3,3>;

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(SurfaceEndoWriter)