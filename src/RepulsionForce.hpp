#ifndef REPULSIONFORCE_HPP_
#define REPULSIONFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"

// Needed here to avoid serialization errors (on Boost<1.37)

/**
 * A morphogen driven cell force class.
 */
template<unsigned DIM>
class RepulsionForce : public AbstractForce<DIM>
{
private:
    double mTreshold;
    friend class boost::serialization::access;
    template<class Archive>

    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
	      archive & mTreshold;
    }

public:

    /**
     * Constructor.
     */
    RepulsionForce(double k=1.0);

    /**
     * Destructor.
     */
    ~RepulsionForce();

    /**
     * Overridden AddForceContribution() method.
     *
     * @param rCellPopulation a cell population object
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RepulsionForce)

#endif /*REPULSIONFORCE_HPP_*/
