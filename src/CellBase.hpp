#ifndef CELLBASE_HPP_
#define CELLBASE_HPP_

#include <boost/shared_ptr.hpp>
#include "AbstractCellProperty.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

class CellBase : public AbstractCellProperty
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

    CellBase(unsigned colour=12);

    virtual ~CellBase();

    unsigned GetColour() const;
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CellBase)

#endif /* CellBase_HPP_ */
