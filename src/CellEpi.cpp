#include "CellEpi.hpp"

CellEpi::CellEpi(unsigned colour)
    : AbstractCellProperty(),
      mColour(colour)
{
}

CellEpi::~CellEpi()
{
}

unsigned CellEpi::GetColour() const
{
    return mColour;
}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(CellEpi)
