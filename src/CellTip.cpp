#include "CellTip.hpp"

CellTip::CellTip(unsigned colour)
    : AbstractCellProperty(),
      mColour(colour)
{
}

CellTip::~CellTip()
{
}

unsigned CellTip::GetColour() const
{
    return mColour;
}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(CellTip)
