#include "CellStalk.hpp"

CellStalk::CellStalk(unsigned colour)
    : AbstractCellProperty(),
      mColour(colour)
{
}

CellStalk::~CellStalk()
{
}

unsigned CellStalk::GetColour() const
{
    return mColour;
}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(CellStalk)
