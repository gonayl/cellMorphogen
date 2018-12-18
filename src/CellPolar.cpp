#include "CellPolar.hpp"

CellPolar::CellPolar(unsigned colour)
    : AbstractCellProperty(),
      mColour(colour)
{
}

CellPolar::~CellPolar()
{
}

unsigned CellPolar::GetColour() const
{
    return mColour;
}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(CellPolar)
