#include "TargetAreaModifier.hpp"
#include "AbstractPhaseBasedCellCycleModel.hpp"
#include "ApoptoticCellProperty.hpp"
#include "CellLumen.hpp"
#include "CellEpi.hpp"
#include "CellEndo.hpp"


template<unsigned DIM>
TargetAreaModifier<DIM>::TargetAreaModifier()
    : AbstractTargetAreaModifier<DIM>(),
      mGrowthDuration(DOUBLE_UNSET)
{
}

template<unsigned DIM>
TargetAreaModifier<DIM>::~TargetAreaModifier()
{
}

template<unsigned DIM>
void TargetAreaModifier<DIM>::UpdateTargetAreaOfCell(CellPtr pCell)
{
    // Get target area A of a healthy cell in S, G2 or M phase
    double cell_target_area = this->mReferenceTargetArea;

    double growth_duration = mGrowthDuration;
    if (growth_duration == DOUBLE_UNSET)
    {
        if (dynamic_cast<AbstractPhaseBasedCellCycleModel*>(pCell->GetCellCycleModel()) == nullptr)
        {
            EXCEPTION("If SetGrowthDuration() has not been called, a subclass of AbstractPhaseBasedCellCycleModel must be used");
        }
        AbstractPhaseBasedCellCycleModel* p_model = static_cast<AbstractPhaseBasedCellCycleModel*>(pCell->GetCellCycleModel());

        growth_duration = p_model->GetG1Duration();

        // If the cell is differentiated then its G1 duration is infinite
        if (growth_duration == DBL_MAX)
        {
            // This is just for fixed cell-cycle models, need to work out how to find the g1 duration
            growth_duration = p_model->GetTransitCellG1Duration();
        }
    }

    if (pCell->HasCellProperty<ApoptoticCellProperty>())
    {
        // Age of cell when apoptosis begins
        if (pCell->GetStartOfApoptosisTime() - pCell->GetBirthTime() < growth_duration)
        {
            cell_target_area *= 0.5*(1 + (pCell->GetStartOfApoptosisTime() - pCell->GetBirthTime())/growth_duration);
        }

        // The target area of an apoptotic cell decreases linearly to zero
        double time_spent_apoptotic = SimulationTime::Instance()->GetTime() - pCell->GetStartOfApoptosisTime();

        cell_target_area *= 1.0 - 0.5/(pCell->GetApoptosisTime())*time_spent_apoptotic;
        if (cell_target_area < 0)
        {
            cell_target_area = 0;
        }
    }
    else if (pCell->HasCellProperty<CellLumen>())
    {
      double cell_age = pCell->GetAge();

      // The target area of a proliferating cell increases linearly from A/2 to A over the course of the prescribed duration
      if (cell_age < growth_duration)
      {
          //cell_target_area *= 0.5*(1 + cell_age/growth_duration);
          cell_target_area = 0.4*this->mReferenceTargetArea;
      }
      else
      {
          /**
           * At division, daughter cells inherit the cell data array from the mother cell.
           * Here, we assign the target area that we want daughter cells to have to cells
           * that we know to divide in this time step.
           *
           * \todo This is a little hack that we might want to clean up in the future.
           */
          if (pCell->ReadyToDivide())
          {
              cell_target_area = 0.3*this->mReferenceTargetArea;
          }
      }
    }
    else if (pCell->HasCellProperty<CellEpi>())
    {
      double cell_age = pCell->GetAge();

      // The target area of a proliferating cell increases linearly from A/2 to A over the course of the prescribed duration
      if (cell_age < growth_duration)
      {
          cell_target_area *= 0.5*(1 + cell_age/growth_duration);
          //cell_target_area = 0.5*this->mReferenceTargetArea;
      }
      else
      {
          /**
           * At division, daughter cells inherit the cell data array from the mother cell.
           * Here, we assign the target area that we want daughter cells to have to cells
           * that we know to divide in this time step.
           *
           * \todo This is a little hack that we might want to clean up in the future.
           */
          if (pCell->ReadyToDivide())
          {
              cell_target_area = 0.5*this->mReferenceTargetArea;
          }
      }
    }
    else if (pCell->HasCellProperty<CellEndo>())
    {
      double cell_age = pCell->GetAge();

      // The target area of a proliferating cell increases linearly from A/2 to A over the course of the prescribed duration
      if (cell_age < growth_duration)
      {
          cell_target_area *= 0.5*(0.5 + cell_age/growth_duration);
          //cell_target_area = 0.5*this->mReferenceTargetArea;
      }
      else
      {
          /**
           * At division, daughter cells inherit the cell data array from the mother cell.
           * Here, we assign the target area that we want daughter cells to have to cells
           * that we know to divide in this time step.
           *
           * \todo This is a little hack that we might want to clean up in the future.
           */
          if (pCell->ReadyToDivide())
          {
              cell_target_area = 0.5*this->mReferenceTargetArea;
          }
      }
    }
    else
    {
        double cell_age = pCell->GetAge();

        // The target area of a proliferating cell increases linearly from A/2 to A over the course of the prescribed duration
        if (cell_age < growth_duration)
        {
            //cell_target_area *= 0.5*(1 + cell_age/growth_duration);
            cell_target_area = 0.8*this->mReferenceTargetArea;
        }
        else
        {
            /**
             * At division, daughter cells inherit the cell data array from the mother cell.
             * Here, we assign the target area that we want daughter cells to have to cells
             * that we know to divide in this time step.
             *
             * \todo This is a little hack that we might want to clean up in the future.
             */
            if (pCell->ReadyToDivide())
            {
                cell_target_area = 0.4*this->mReferenceTargetArea;
            }
        }
    }

    // Set cell data
    pCell->GetCellData()->SetItem("target area", cell_target_area);
}

template<unsigned DIM>
double TargetAreaModifier<DIM>::GetGrowthDuration()
{
    return mGrowthDuration;
}

template<unsigned DIM>
void TargetAreaModifier<DIM>::SetGrowthDuration(double growthDuration)
{
    assert(growthDuration >= 0.0);
    mGrowthDuration = growthDuration;
}

template<unsigned DIM>
void TargetAreaModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<GrowthDuration>" << mGrowthDuration << "</GrowthDuration>\n";

    // Next, call method on direct parent class
    AbstractTargetAreaModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class TargetAreaModifier<1>;
template class TargetAreaModifier<2>;
template class TargetAreaModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(TargetAreaModifier)
