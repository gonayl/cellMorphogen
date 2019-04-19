#include "CellPeriph.hpp"

CellPeriph::CellPeriph(unsigned colour)
    : AbstractCellProperty(),
      mColour(colour)
{
}

CellPeriph::~CellPeriph()
{
}

unsigned CellPeriph::GetColour() const
{
    return mColour;
}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(CellPeriph)
