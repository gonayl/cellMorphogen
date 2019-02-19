#ifndef TARGETAREAMODIFIER_HPP_
#define TARGETAREAMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "AbstractTargetAreaModifier.hpp"

/**
 * A target area modifier class in which the target area of a cell grows linearly, up to
 * mReferenceTargetArea, over a prescribed duration.
 *
 * If used with a phase-based cell-cycle model (such as FixedG1GenerationalCellCycleModel),
 * the target area of a cell increases linearly from the value 0.5*mReferenceTargetArea
 * up to mReferenceTargetArea over the course of the cell's G1 phase.
 *
 * If used with a non-phase-based cell-cycle model, the target area of a cell increases
 * linearly from the value 0.5*mReferenceTargetArea up to mReferenceTargetArea while the
 * cell's age is less than mGrowthDuration.
 *
 * Here mReferenceTargetArea and mGrowthDuration are settable member variables. The default
 * value of mReferenceTargetArea is 1.0 and the default value of mGrowthDuration is DOUBLE_UNSET.
 *
 * Note that if mGrowthDuration is set by the user, then this value is are used to prescribe
 * target area growth as described earlier, regardless of whether a phase-based cell-cycle
 * model is present.
 */
template<unsigned DIM>
class TargetAreaModifier : public AbstractTargetAreaModifier<DIM>
{
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractTargetAreaModifier<DIM> >(*this);
        archive & mGrowthDuration;
    }

    /**
     * The duration that a cell's target area to increase from 0.5*mReferenceTargetArea
     * to mReferenceTargetArea at the start of its cell cycle. Defaults to DOUBLE_UNSET.
     * If this variable is set using SetGrowthDuration(), then it is used regardless of
     * whether a phase-based cell-cycle model is used.
     */
    double mGrowthDuration;

public:

    /**
     * Default constructor.
     */
    TargetAreaModifier();

    /**
     * Destructor.
     */
    virtual ~TargetAreaModifier();

    /**
     * Overridden UpdateTargetAreaOfCell() method.
     *
     * @param pCell pointer to the cell
     */
    virtual void UpdateTargetAreaOfCell(const CellPtr pCell);

    /**
     * @return #mGrowthDuration
     */
    double GetGrowthDuration();

    /**
     * Set #mGrowthDuration.
     *
     * @param growthDuration the new value of #mGrowthDuration
     */
    void SetGrowthDuration(double growthDuration);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(TargetAreaModifier)

#endif /*TARGETAREAMODIFIER_HPP_*/
