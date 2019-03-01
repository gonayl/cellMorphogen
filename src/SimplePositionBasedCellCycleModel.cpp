#include "SimplePositionBasedCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "CellEpi.hpp"
#include "CellEndo.hpp"
#include "CellLumen.hpp"

SimplePositionBasedCellCycleModel::SimplePositionBasedCellCycleModel()
{
}

void SimplePositionBasedCellCycleModel::UpdateCellCyclePhase()
{
    // mG1Duration is set when the cell-cycle model is given a cell

    if (mpCell->HasCellProperty<CellEpi>())
    {
        // Get cell's distance to boundary
        double dist2boundary = mpCell->GetCellData()->GetItem("mindistborder");

        AbstractSimplePhaseBasedCellCycleModel::UpdateCellCyclePhase();

        //if (mCurrentCellCyclePhase == G_ONE_PHASE)
        //{
            // Update G1 duration based on the position of the cell in the mass
            // double dt = SimulationTime::Instance()->GetTimeStep();
            RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

            mG1Duration += 10*p_gen->ranf()*dist2boundary;
            mTransitCellG1Duration += 10*p_gen->ranf()*dist2boundary; 

            std::cout << 10*p_gen->ranf()*dist2boundary << std::endl ;

        //}
    }
}

SimplePositionBasedCellCycleModel::SimplePositionBasedCellCycleModel(const SimplePositionBasedCellCycleModel& rModel)
   : AbstractSimplePhaseBasedCellCycleModel(rModel)
{
    /*
     * Initialize only those member variables defined in this class.
     *
     * The member variables mCurrentCellCyclePhase, mG1Duration,
     * mMinimumGapDuration, mStemCellG1Duration, mTransitCellG1Duration,
     * mSDuration, mG2Duration and mMDuration are initialized in the
     * AbstractPhaseBasedCellCycleModel constructor.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     *
     * Note that mG1Duration is (re)set as soon as InitialiseDaughterCell()
     * is called on the new cell-cycle model.
     */
}

AbstractCellCycleModel* SimplePositionBasedCellCycleModel::CreateCellCycleModel()
{
    return new SimplePositionBasedCellCycleModel(*this);
}

void SimplePositionBasedCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{

    // No new parameters to output, simply call method on direct parent class
    AbstractSimplePhaseBasedCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(SimplePositionBasedCellCycleModel)
