#ifndef PERIMETERDEPENDENTCELLCYCLEMODEL_HPP_
#define PERIMETERDEPENDENTCELLCYCLEMODEL_HPP_

#include "AbstractCellCycleModel.hpp"

/**
 * Simple cell-cycle model where mature non-differentiated cells have a specified probability of
 * dividing per hour.
 *
 * The class includes two parameters: the first, mMaxStretch, defines the probability
 * of dividing per hour; the second, mMinimumDivisionAge, defines a minimum age at which cells
 * may divide.
 */
class PerimeterDependentCellCycleModel : public AbstractCellCycleModel
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
        archive & boost::serialization::base_object<AbstractCellCycleModel>(*this);
        archive & mMaxStretch;
        archive & mMinimumDivisionAge;
    }

protected:
    /**
     * Maximum stetching allowed.
     * Defaults to ??.
     */
    double mMaxStretch;

    /**
     * Minimum age of a cell at which it may divide.
     * Defaults to 1 hour.
     */
    double mMinimumDivisionAge;

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
    PerimeterDependentCellCycleModel(const PerimeterDependentCellCycleModel& rModel);

public:

    /**
     * Constructor.
     */
    PerimeterDependentCellCycleModel();

    /**
     * Overridden ReadyToDivide() method.
     *
     * If the cell's age is greater than mMinimumDivisionAge, then we draw a uniform
     * random number r ~ U[0,1]. If r < mMaxStretch*dt, where dt is the
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
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * Set the value of mMaxStretch.
     *
     * @param divisionProbability the new value of mMaxStretch
     */
    void SetMaxStretch(double divisionProbability);

    /**
     * Get mMaxStretch.
     *
     * @return mMaxStretch
     */
    double GetMaxStretch();

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
     * Overridden GetAverageTransitCellCycleTime() method.
     *
     * @return the average cell cycle time for cells of transit proliferative type
     */
    double GetAverageTransitCellCycleTime();

    /**
     * Overridden GetAverageStemCellCycleTime() method.
     *
     * @return the average cell cycle time for cells of stem proliferative type
     */
    double GetAverageStemCellCycleTime();

    /**
     * Overridden OutputCellCycleModelParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(PerimeterDependentCellCycleModel)

#endif // PERIMETERDEPENDENTCELLCYCLEMODEL_HPP_
