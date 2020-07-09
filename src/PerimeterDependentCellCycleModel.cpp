#include "PerimeterDependentCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellStalk.hpp"

PerimeterDependentCellCycleModel::PerimeterDependentCellCycleModel()
    : AbstractCellCycleModel(),
      mMaxStretch(3.2),
      mMinimumDivisionAge(0.2)
{
}

PerimeterDependentCellCycleModel::PerimeterDependentCellCycleModel(const PerimeterDependentCellCycleModel& rModel)
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

bool PerimeterDependentCellCycleModel::ReadyToDivide()
{
    assert(mpCell != nullptr);

    if (!mReadyToDivide)
    {
        if (GetAge() > mMinimumDivisionAge)
        {
            // double dt = SimulationTime::Instance()->GetTimeStep();
            if (!(mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>()))
            {
                double cell_elongation = mpCell->GetCellData()->GetItem("perimeter");
                //double have_tip_neighboor = mpCell->GetCellData()->GetItem("have_tip_neighboor");
              //  if (cell_elongation > mMaxStretch && have_tip_neighboor > 0)
                if (cell_elongation > mMaxStretch )
                {
                    mReadyToDivide = true;
                }
            }
        }
    }
    return mReadyToDivide;
}

AbstractCellCycleModel* PerimeterDependentCellCycleModel::CreateCellCycleModel()
{
    return new PerimeterDependentCellCycleModel(*this);
}

void PerimeterDependentCellCycleModel::SetMaxStretch(double divisionProbability)
{
    mMaxStretch = divisionProbability;
}

double PerimeterDependentCellCycleModel::GetMaxStretch()
{
    return mMaxStretch;
}

void PerimeterDependentCellCycleModel::SetMinimumDivisionAge(double minimumDivisionAge)
{
    mMinimumDivisionAge = minimumDivisionAge;
}

double PerimeterDependentCellCycleModel::GetMinimumDivisionAge()
{
    return mMinimumDivisionAge;
}

double PerimeterDependentCellCycleModel::GetAverageTransitCellCycleTime()
{
    return 1.0/mMaxStretch;
}

double PerimeterDependentCellCycleModel::GetAverageStemCellCycleTime()
{
    return 1.0/mMaxStretch;
}

void PerimeterDependentCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MaxStretch>" << mMaxStretch << "</MaxStretch>\n";
    *rParamsFile << "\t\t\t<MinimumDivisionAge>" << mMinimumDivisionAge << "</MinimumDivisionAge>\n";

    // Call method on direct parent class
    AbstractCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(PerimeterDependentCellCycleModel)
