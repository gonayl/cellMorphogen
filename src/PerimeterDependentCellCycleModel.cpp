#include "PerimeterDependentCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellStalk.hpp"
#include "CellPeriph.hpp"
#include "CellCore.hpp"
#include "CellVessel.hpp"
using namespace std ;

PerimeterDependentCellCycleModel::PerimeterDependentCellCycleModel()
    : AbstractCellCycleModel(),
      mMaxStretch(2.3),
      mMaxStretchPeriph(7.0),
      mMinimumDivisionAge(0.1)
{
}

PerimeterDependentCellCycleModel::PerimeterDependentCellCycleModel(const PerimeterDependentCellCycleModel& rModel)
   : AbstractCellCycleModel(rModel),
     mMaxStretch(rModel.mMaxStretch),
     mMaxStretchPeriph(rModel.mMaxStretchPeriph),
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

                double have_tip_neighboor = mpCell->GetCellData()->GetItem("have_tip_neighboor");
                double have_vessel_neighboor = mpCell->GetCellData()->GetItem("have_vessel_neighboor");
                double cell_elongation = mpCell->GetCellData()->GetItem("perimeter");
                double cell_maxmin = mpCell->GetCellData()->GetItem("maxmin");
                bool is_vessel = mpCell->template HasCellProperty<CellVessel>();

                if (cell_elongation > mMaxStretch && have_tip_neighboor > 0)
                //if (cell_elongation > 1.8 && cell_maxmin > mMaxStretchPeriph )
                {
                    mReadyToDivide = true;
                }

                else if (cell_elongation > 1.0 && cell_maxmin > 8.0 && have_tip_neighboor > 0)
                {
                    cout << cell_maxmin << endl ;
                    mReadyToDivide = true;
                }
                else if (cell_elongation > mMaxStretch && is_vessel)
                {
                    mReadyToDivide = true;
                }
                else if (cell_elongation > mMaxStretch && have_vessel_neighboor > 0)
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

void PerimeterDependentCellCycleModel::SetMaxStretchPeriph(double divisionProbability)
{
    mMaxStretchPeriph = divisionProbability;
}

double PerimeterDependentCellCycleModel::GetMaxStretch()
{
    return mMaxStretch;
}

double PerimeterDependentCellCycleModel::GetMaxStretchPeriph()
{
    return mMaxStretchPeriph;
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
    *rParamsFile << "\t\t\t<MaxStretchPeriph>" << mMaxStretchPeriph << "</MaxStretchPeriph>\n";
    *rParamsFile << "\t\t\t<MinimumDivisionAge>" << mMinimumDivisionAge << "</MinimumDivisionAge>\n";

    // Call method on direct parent class
    AbstractCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(PerimeterDependentCellCycleModel)
