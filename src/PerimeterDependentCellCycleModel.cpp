/*

Copyright (c) 2005-2017, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "PerimeterDependentCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

PerimeterDependentCellCycleModel::PerimeterDependentCellCycleModel()
    : AbstractCellCycleModel(),
      mMaxStretch(2.2),
      mMinimumDivisionAge(1.0)
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
                if (cell_elongation > mMaxStretch)
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
