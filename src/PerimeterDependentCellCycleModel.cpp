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
      mMinimumDivisionAge(0.001)
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
                // double have_vessel_neighboor = mpCell->GetCellData()->GetItem("have_vessel_neighboor");
                double cell_elongation = mpCell->GetCellData()->GetItem("perimeter");
                double cell_maxmin = mpCell->GetCellData()->GetItem("maxmin");
                double cell_min = mpCell->GetCellData()->GetItem("min");
                bool is_vessel = mpCell->template HasCellProperty<CellVessel>();

                if (cell_elongation > mMaxStretch && have_tip_neighboor > 0) // following stalk
                //if (cell_elongation > 1.8 && cell_maxmin > mMaxStretchPeriph )
                {
                    mReadyToDivide = true;
                }

                else if (cell_elongation > 1.0 && cell_maxmin > 8.0 && have_tip_neighboor > 0) // following stalk but with minmax
                {
                    cout << cell_maxmin << endl ;
                    mReadyToDivide = true;
                }

                else if (cell_elongation > 2.8 && is_vessel) // vessel
                {
                    mReadyToDivide = true;
                }

                else if (cell_elongation > 2.0 && cell_maxmin > 8.0 && is_vessel) // vessel
                {
                    mReadyToDivide = true;
                }

                else if (cell_elongation > 2.1 && cell_min < 0.11 && is_vessel) // if vessel nieghbour should be also vessel
                {
                    mReadyToDivide = true;
                }
                else if (cell_elongation > 2.3 && cell_min < 0.15 && is_vessel) // if vessel nieghbour should be also vessel
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
