#include "CellBoundary.hpp"

CellBoundary::CellBoundary(unsigned colour)
    : AbstractCellProperty(),
      mColour(colour)
{
}

CellBoundary::~CellBoundary()
{
}

unsigned CellBoundary::GetColour() const
{
    return mColour;
}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(CellBoundary)
