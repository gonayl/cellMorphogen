#include "CellMotile.hpp"

CellMotile::CellMotile(unsigned colour)
    : AbstractCellProperty(),
      mColour(colour)
{
}

CellMotile::~CellMotile()
{
}

unsigned CellMotile::GetColour() const
{
    return mColour;
}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(CellMotile)
