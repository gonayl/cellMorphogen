#ifndef CELLEPI_HPP_
#define CELLEPI_HPP_

#include <boost/shared_ptr.hpp>
#include "AbstractCellProperty.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

class CellEpi : public AbstractCellProperty
{
protected:

    unsigned mColour;

private:

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellProperty>(*this);
        archive & mColour;
    }

public:

    CellEpi(unsigned colour=1);

    virtual ~CellEpi();

    unsigned GetColour() const;
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CellEpi)

#endif /* CELLEPI_HPP_ */
