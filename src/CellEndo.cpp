#include "CellEndo.hpp"

CellEndo::CellEndo(unsigned colour)
    : AbstractCellProperty(),
      mColour(colour)
{
}

CellEndo::~CellEndo()
{
}

unsigned CellEndo::GetColour() const
{
    return mColour;
}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(CellEndo)
