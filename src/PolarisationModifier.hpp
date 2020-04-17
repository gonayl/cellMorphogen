#ifndef POLARISATIONMODIFIER_HPP_
#define POLARISATIONMODIFIER_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractCellBasedSimulationModifier.hpp"

/**
 * A modifier class which Compute the size of the polarisation vector of epithelial cells
 */
template<unsigned DIM>
class PolarisationModifier : public AbstractCellBasedSimulationModifier<DIM,DIM>
{
private:

    double mVecPolarisationDecrease ;
    double mEpiEpiPolarisationParameter ;
    double mEndoEpiPolarisationParameter ;
    double mLumenEpiPolarisationParameter ;

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
        archive & boost::serialization::base_object<AbstractCellBasedSimulationModifier<DIM,DIM> >(*this);
        archive & mVecPolarisationDecrease;
        archive & mEpiEpiPolarisationParameter;
        archive & mEndoEpiPolarisationParameter;
        archive & mLumenEpiPolarisationParameter;
    }

public:

    /**
     * Default constructor.
     */
    PolarisationModifier();

    /**
     * Destructor.
     */
    virtual ~PolarisationModifier();

    /**
     * Overridden UpdateAtEndOfTimeStep() method.
     *
     * Specify what to do in the simulation at the end of each time step.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    /**
     * Overridden SetupSolve() method.
     *
     * Specify what to do in the simulation before the start of the time loop.
     *
     * @param rCellPopulation reference to the cell population
     * @param outputDirectory the output directory, relative to where Chaste output is stored
     */
    virtual void SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory);

    /**
     * Helper method to compute the perimeter of each cell in the population and store these in the CellData.
     *
     * @param rCellPopulation reference to the cell population
     */
    void UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation);

    double GetVecPolarisationDecrease();

    double GetEpiEpiPolarisationParameter();

    double GetEndoEpiPolarisationParameter();

    double GetLumenEpiPolarisationParameter();

    void SetVecPolarisationDecrease(double vecPolarisationDecrease);

    void SetEpiEpiPolarisationParameter(double epiEpiPolarisationParameter);

    void SetEndoEpiPolarisationParameter(double endoEpiPolarisationParameter);

    void SetLumenEpiPolarisationParameter(double lumenEpiPolarisationParameter);

    /**
     * Overridden OutputSimulationModifierParameters() method.
     * Output any simulation modifier parameters to file.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputSimulationModifierParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PolarisationModifier)

#endif /*PolarisationModifier_HPP_*/
