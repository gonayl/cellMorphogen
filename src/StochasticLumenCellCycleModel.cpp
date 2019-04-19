#include "StochasticLumenCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellLumen.hpp"
#include "CellEpi.hpp"
#include "CellEndo.hpp"

StochasticLumenCellCycleModel::StochasticLumenCellCycleModel()
    : AbstractSimpleGenerationalCellCycleModel(),
      mMinimumDivisionAge(1.0),
      mMaxTransitGeneration(1.0)
{
}

StochasticLumenCellCycleModel::StochasticLumenCellCycleModel(const StochasticLumenCellCycleModel& rModel)
   : AbstractSimpleGenerationalCellCycleModel(rModel),
     mMinimumDivisionAge(rModel.mMinimumDivisionAge),
     mMaxTransitGeneration(rModel.mMaxTransitGeneration)
{
    /*
     * Initialize only those member variables defined in this class.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     */
}

bool StochasticLumenCellCycleModel::ReadyToDivide()
{
    assert(mpCell != nullptr);

    if (!mReadyToDivide)
    {
      if (mpCell->HasCellProperty<CellLumen>())
      {
        mMinimumDivisionAge = 1.0 ;
      }
      else if (mpCell->HasCellProperty<CellEpi>())
      {
        mMinimumDivisionAge = 3.0 ;
      }
      else if (mpCell->HasCellProperty<CellEndo>())
      {
        mMinimumDivisionAge = 5.0 ;
      }
        if (GetAge() > mMinimumDivisionAge)
        {
            if (!(mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>()))
            {
                mReadyToDivide = true;
            }
        }
    }
    return mReadyToDivide;
}

AbstractSimpleGenerationalCellCycleModel* StochasticLumenCellCycleModel::CreateCellCycleModel()
{
    return new StochasticLumenCellCycleModel(*this);
}

void StochasticLumenCellCycleModel::SetG1Duration()
{
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
    // double weight = mpCell->GetCellData()->GetItem("mindistborder") ;
    double weight = 1.0 ;

    assert(mpCell != nullptr);

    if (mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
    {
        mG1Duration = GetStemCellG1Duration() + 4*p_gen->ranf(); // U[14,18] for default parameters (mStemCellG1Duration) according to Meineke
    }
    else if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
    {
      if (mpCell->HasCellProperty<CellLumen>())
      {
        mG1Duration = 0.5 * GetTransitCellG1Duration() + 1*p_gen->ranf();
      }
      else if (mpCell->HasCellProperty<CellEpi>())
      {
        mG1Duration = GetTransitCellG1Duration() + 2*p_gen->ranf() + weight ;
      }
      else if (mpCell->HasCellProperty<CellEndo>())
      {
        mG1Duration = GetTransitCellG1Duration() + 8*p_gen->ranf();
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

void StochasticLumenCellCycleModel::SetMaxTransitGeneration()
{
    assert(mpCell != nullptr);

    if (mpCell->HasCellProperty<CellLumen>())
      {
        mMaxTransitGeneration = 2.0;
      }
    else if (mpCell->HasCellProperty<CellEpi>())
      {
        mMaxTransitGeneration = 4.0;
      }
    else if (mpCell->HasCellProperty<CellEndo>())
      {
        mMaxTransitGeneration = 1.0;
      }
}

void StochasticLumenCellCycleModel::SetMinimumDivisionAge(double minimumDivisionAge)
{
    mMinimumDivisionAge = minimumDivisionAge;
}

double StochasticLumenCellCycleModel::GetMinimumDivisionAge()
{
    return mMinimumDivisionAge;
}

void StochasticLumenCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MinimumDivisionAge>" << mMinimumDivisionAge << "</MinimumDivisionAge>\n";

    // Call method on direct parent class
    AbstractSimpleGenerationalCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(StochasticLumenCellCycleModel)
