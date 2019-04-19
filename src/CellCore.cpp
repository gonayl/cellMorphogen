#include "CellCore.hpp"

CellCore::CellCore(unsigned colour)
    : AbstractCellProperty(),
      mColour(colour)
{
}

CellCore::~CellCore()
{
}

unsigned CellCore::GetColour() const
{
    return mColour;
}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(CellCore)
