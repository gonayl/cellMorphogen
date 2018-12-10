#ifndef STOCHASTICLUMENCELLCYCLEMODEL_HPP_
#define STOCHASTICLUMENCELLCYCLEMODEL_HPP_

#include "AbstractSimpleGenerationalCellCycleModel.hpp"

/**
 * Simple cell-cycle model where mature non-differentiated cells have a specified probability of
 * dividing per hour.
 *
 * The class includes two parameters: the first, mDivisionProbability, defines the probability
 * of dividing per hour; the second, mMinimumDivisionAge, defines a minimum age at which cells
 * may divide.
 */
class StochasticLumenCellCycleModel : public AbstractSimpleGenerationalCellCycleModel
{
private:

    friend class boost::serialization::access;

    /**
     * Boost Serialization method for archiving/checkpointing
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleGenerationalCellCycleModel>(*this);
        archive & mMinimumDivisionAge;
        archive & mMaxTransitGeneration;
    }

protected:

    /**
     * Minimum age of a cell at which it may divide.
     * Defaults to 1 hour.
     */
    double mMinimumDivisionAge;

    double mMaxTransitGeneration ;

    void SetMaxTransitGeneration() ;

    /**
     * Set the duration of G1 phase. This method is called on each cell at the
     * start of a simulation, and for both daughter cells immediately following
     * cell division.
     *
     * If the cell associated with this cell-cycle model has stem proliferative
     * type, then the G1 phase duration is drawn from the uniform distribution
     * U[14,18]. If the cell has transit proliferative type (semi-differentiated),
     * then the G1 phase duration is drawn from the uniform distribution U[4,6].
     * These two distributions, proposed by Meineke et al (doi:10.1046/j.0960-7722.2001.00216.x),
     * reflect indirect biological observations that stem cells cycle more
     * slowly than their progeny.
     *
     * If the cell is differentiated, then the G1 phase duration is set to DBL_MAX,
     * so that the cell will never reach the end of G1 phase.
     */
    void SetG1Duration();

    /**
     * Protected copy-constructor for use by CreateCellCycleModel.
     *
     * The only way for external code to create a copy of a cell cycle model
     * is by calling that method, to ensure that a model of the correct subclass
     * is created. This copy-constructor helps subclasses to ensure that all
     * member variables are correctly copied when this happens.
     *
     * This method is called by child classes to set member variables for a
     * daughter cell upon cell division. Note that the parent cell cycle model
     * will have had ResetForDivision() called just before CreateCellCycleModel()
     * is called, so performing an exact copy of the parent is suitable behaviour.
     * Any daughter-cell-specific initialisation can be done in InitialiseDaughterCell().
     *
     * @param rModel the cell cycle model to copy.
     */

    StochasticLumenCellCycleModel(const StochasticLumenCellCycleModel& rModel);

public:

    /**
     * Constructor.
     */
    StochasticLumenCellCycleModel();

    /**
     * Overridden ReadyToDivide() method.
     *
     * If the cell's age is greater than mMinimumDivisionAge, then we draw a uniform
     * random number r ~ U[0,1]. If r < mDivisionProbability*dt, where dt is the
     * simulation time step, then the cell is ready to divide and we return true.
     * Otherwise, the cell is not yet ready to divide and we return false.
     *
     * @return whether the cell is ready to divide.
     */
    virtual bool ReadyToDivide();


    /**
     * Overridden builder method to create new instances of
     * the cell-cycle model.
     *
     * @return new cell-cycle model
     */
    AbstractSimpleGenerationalCellCycleModel* CreateCellCycleModel();

    /**
     * Set the value of mMinimumDivisionAge.
     *
     * @param minimumDivisionAge the new value of mMinimumDivisionAge
     */
    void SetMinimumDivisionAge(double minimumDivisionAge);

    /**
     * Get mMinimumDivisionAge.
     *
     * @return mMinimumDivisionAge
     */
    double GetMinimumDivisionAge();

    /**
     * Overridden OutputCellCycleModelParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(StochasticLumenCellCycleModel)

#endif // STOCHASTICLUMENCELLCYCLEMODEL_HPP_
