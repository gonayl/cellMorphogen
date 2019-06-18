#include "UniformG1GenerationalBoundaryCellCycleModel.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellCore.hpp"
#include "CellPeriph.hpp"
#include "CellLumen.hpp"

UniformG1GenerationalBoundaryCellCycleModel::UniformG1GenerationalBoundaryCellCycleModel()
{
}

UniformG1GenerationalBoundaryCellCycleModel::UniformG1GenerationalBoundaryCellCycleModel(const UniformG1GenerationalBoundaryCellCycleModel& rModel)
   : AbstractSimpleGenerationalCellCycleModel(rModel)
{
    /*
     * The member variables mGeneration and mMaxTransitGeneration are
     * initialized in the AbstractSimpleGenerationalCellCycleModel
     * constructor.
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

AbstractCellCycleModel* UniformG1GenerationalBoundaryCellCycleModel::CreateCellCycleModel()
{
    return new UniformG1GenerationalBoundaryCellCycleModel(*this);
}

void UniformG1GenerationalBoundaryCellCycleModel::SetG1Duration()
{
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    assert(mpCell != nullptr);

    if (mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
    {
        mG1Duration = GetStemCellG1Duration() + 4*p_gen->ranf(); // U[14,18] for default parameters (mStemCellG1Duration) according to Meineke
    }
    else if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
    {
        if (mpCell->HasCellProperty<CellCore>())
        {
          mG1Duration =  10 ; // ICI : MODIFIER DUREE PHASE G1 POUR CELLULES DE COEUR, tu peux remplacer GetTransitCellG1Duration() par un nombre
          mSDuration =  10 ;
          mG2Duration = 3 ;
          mMDuration =  3 ;
        }
        else if (mpCell->HasCellProperty<CellPeriph>())
        {
          mG1Duration = 5 ; // ICI : MODIFIER DUREE PHASE G1 POUR CELLULE PERIPH
          mSDuration =  5 ;
          mG2Duration = 1 ;
          mMDuration =  1 ;
        }
        else if (mpCell->HasCellProperty<CellLumen>())
        {
          mG1Duration = 2 ; // ICI : MODIFIER DUREE PHASE G1 POUR CELLULE PERIPH
          mSDuration =  2 ;
          mG2Duration = 1 ;
          mMDuration =  1 ;
        }
        else
        {
          mG1Duration = GetTransitCellG1Duration() + 3*p_gen->ranf();
        }
    }
    else if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
    {
        mG1Duration = DBL_MAX;
    }
    else
    {
        NEVER_REACHED;
    }
}

void UniformG1GenerationalBoundaryCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractSimpleGenerationalCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(UniformG1GenerationalBoundaryCellCycleModel)
