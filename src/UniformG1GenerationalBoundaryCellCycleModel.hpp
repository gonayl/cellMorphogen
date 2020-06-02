#ifndef UNIFORMG1GENERATIONALBOUNDARYCELLCYCLEMODEL_HPP_
#define UNIFORMG1GENERATIONALBOUNDARYCELLCYCLEMODEL_HPP_

#include "AbstractSimpleGenerationalCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"

/**
 * A stochastic cell-cycle model employed by Meineke et al (2001) in their off-lattice
 * model of the intestinal crypt (doi:10.1046/j.0960-7722.2001.00216.x). Cells that
 * are proliferating are assigned G1 phase durations drawn from a uniform distribution,
 * which differs between stem and transit amplifying cells. All other cell-cycle phases
 * are held constant.
 */
class UniformG1GenerationalBoundaryCellCycleModel : public AbstractSimpleGenerationalCellCycleModel
{
    friend class TestSimpleCellCycleModels;

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell-cycle model and random number generator, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleGenerationalCellCycleModel>(*this);

        // Make sure the RandomNumberGenerator singleton gets saved too
        SerializableSingleton<RandomNumberGenerator>* p_wrapper = RandomNumberGenerator::Instance()->GetSerializationWrapper();
        archive & p_wrapper;
        archive & mCycleDuration;
    }

protected:

    double mCycleDuration ;

    void SetG1Duration();


    UniformG1GenerationalBoundaryCellCycleModel(const UniformG1GenerationalBoundaryCellCycleModel& rModel);

public:

    /**
     * Constructor - just a default, mBirthTime is set in the AbstractCellCycleModel class.
     * mG1Duration is set very high, it is set for the individual cells when InitialiseDaughterCell is called
     */
    UniformG1GenerationalBoundaryCellCycleModel();

    AbstractCellCycleModel* CreateCellCycleModel();

    virtual bool ReadyToDivide();

    void SetCycleDuration(double cycleduration) ;

    double GetCycleDuration() const ;

    /**
     * Overridden OutputCellCycleModelParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(UniformG1GenerationalBoundaryCellCycleModel)

#endif /*UNIFORMG1GENERATIONALBOUNDARYCELLCYCLEMODEL_HPP_*/
