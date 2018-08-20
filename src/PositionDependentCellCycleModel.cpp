#include "PositionDependentCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CounterSingleton.hpp"
#include <cmath>

PositionDependentCellCycleModel::PositionDependentCellCycleModel()
    : AbstractCellCycleModel(),
      mMaxStretch(2.2),
      mMinimumDivisionAge(1.0)
{
}

PositionDependentCellCycleModel::PositionDependentCellCycleModel(const PositionDependentCellCycleModel& rModel)
   : AbstractCellCycleModel(rModel),
     mMaxStretch(rModel.mMaxStretch),
     mMinimumDivisionAge(rModel.mMinimumDivisionAge)
{
    /*
     * Initialize only those member variables defined in this class.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     */
}

bool PositionDependentCellCycleModel::ReadyToDivide()
{
    assert(mpCell != nullptr);

    if (!mReadyToDivide)
    {
        if (GetAge() > mMinimumDivisionAge)
        {
            // double dt = SimulationTime::Instance()->GetTimeStep();
            if (!(mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>()))
            {
                double cell_radius = mpCell->GetCellData()->GetItem("radius");
                double dist = mpCell->GetCellData()->GetItem("dist");
                double xpos = mpCell->GetCellData()->GetItem("xpos") ;
                double ypos = mpCell->GetCellData()->GetItem("ypos") ;

                if (dist > cell_radius)
                {
                    mReadyToDivide = true;
                    mpCell->GetCellData()->SetItem("xposini", xpos) ;
                    mpCell->GetCellData()->SetItem("yposini", ypos) ;
                    CounterSingleton::Instance()->IncrementCounter(); 
                }
            }
        }
    }
    return mReadyToDivide;
}

AbstractCellCycleModel* PositionDependentCellCycleModel::CreateCellCycleModel()
{
    return new PositionDependentCellCycleModel(*this);
}

void PositionDependentCellCycleModel::SetMaxStretch(double divisionProbability)
{
    mMaxStretch = divisionProbability;
}

double PositionDependentCellCycleModel::GetMaxStretch()
{
    return mMaxStretch;
}

void PositionDependentCellCycleModel::SetMinimumDivisionAge(double minimumDivisionAge)
{
    mMinimumDivisionAge = minimumDivisionAge;
}

double PositionDependentCellCycleModel::GetMinimumDivisionAge()
{
    return mMinimumDivisionAge;
}

double PositionDependentCellCycleModel::GetAverageTransitCellCycleTime()
{
    return 1.0/mMaxStretch;
}

double PositionDependentCellCycleModel::GetAverageStemCellCycleTime()
{
    return 1.0/mMaxStretch;
}

void PositionDependentCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MaxStretch>" << mMaxStretch << "</MaxStretch>\n";
    *rParamsFile << "\t\t\t<MinimumDivisionAge>" << mMinimumDivisionAge << "</MinimumDivisionAge>\n";

    // Call method on direct parent class
    AbstractCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(PositionDependentCellCycleModel)
