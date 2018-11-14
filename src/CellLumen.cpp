#include "CellLumen.hpp"

CellLumen::CellLumen(unsigned colour)
    : AbstractCellProperty(),
      mColour(colour)
{
}

CellLumen::~CellLumen()
{
}

unsigned CellLumen::GetColour() const
{
    return mColour;
}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(CellLumen)
