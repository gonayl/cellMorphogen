#include "CellBase.hpp"

CellBase::CellBase(unsigned colour)
    : AbstractCellProperty(),
      mColour(colour)
{
}

CellBase::~CellBase()
{
}

unsigned CellBase::GetColour() const
{
    return mColour;
}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(CellBase)
