#ifndef CELLLUMEN_HPP_
#define CELLLUMEN_HPP_

#include <boost/shared_ptr.hpp>
#include "AbstractCellProperty.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

class CellLumen : public AbstractCellProperty
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

    CellLumen(unsigned colour=4);

    virtual ~CellLumen();


    unsigned GetColour() const;
};

#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(CellLumen)

#endif /* CELLLUMEN_HPP_ */
