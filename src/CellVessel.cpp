#include "CellVessel.hpp"

CellVessel::CellVessel(unsigned colour)
    : AbstractCellProperty(),
      mColour(colour)
{
}

CellVessel::~CellVessel()
{
}

unsigned CellVessel::GetColour() const
{
    return mColour;
}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(CellVessel)
