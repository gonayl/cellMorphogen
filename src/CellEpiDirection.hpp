#ifndef CELLEPIDIR_HPP_
#define CELLEPIDIR_HPP_

#include <boost/shared_ptr.hpp>
#include "AbstractCellProperty.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

class CellEpiDirection : public AbstractCellProperty
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

    CellEpiDirection(unsigned colour=5);

    virtual ~CellEpiDirection();


    unsigned GetColour() const;
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CellEpiDirection)

#endif
