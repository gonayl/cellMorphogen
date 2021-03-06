#ifndef CELLVESSEL_HPP_
#define CELLVESSEL_HPP_

#include <boost/shared_ptr.hpp>
#include "AbstractCellProperty.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * Cell label class.
 *
 * Each Cell owns a CellPropertyCollection, which may include a shared pointer
 * to an object of this type. When a Cell that is labelled divides, the daughter
 * cells are both labelled.
 *
 * The CellVessel object keeps track of the number of cells that have the label, as well
 * as what colour should be used by the visualizer to display cells with the label.
 */
class CellVessel : public AbstractCellProperty
{
protected:

    /**
     * Colour for use by visualizer.
     */
    unsigned mColour;

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellProperty>(*this);
        archive & mColour;
    }

public:

    /**
     * Constructor.
     *
     * @param colour  what colour cells with this label should be in the visualizer (defaults to 5)
     */
    CellVessel(unsigned colour=8);

    /**
     * Destructor.
     */
    virtual ~CellVessel();

    /**
     * @return #mColour.
     */
    unsigned GetColour() const;
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(CellVessel)

#endif /* CELLVESSEL_HPP_ */
