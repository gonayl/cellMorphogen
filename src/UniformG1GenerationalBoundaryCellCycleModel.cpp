#include "UniformG1GenerationalBoundaryCellCycleModel.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "CellCore.hpp"
#include "CellPeriph.hpp"
#include "SimulationParameters.hpp"
#include "CellLumen.hpp"
#include "CellEpi.hpp"
#include <math.h>

UniformG1GenerationalBoundaryCellCycleModel::UniformG1GenerationalBoundaryCellCycleModel()
: mCycleDuration(11)
{
}

UniformG1GenerationalBoundaryCellCycleModel::UniformG1GenerationalBoundaryCellCycleModel(const UniformG1GenerationalBoundaryCellCycleModel& rModel)
   : AbstractSimpleGenerationalCellCycleModel(rModel),
     mCycleDuration(rModel.mCycleDuration)
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

void UniformG1GenerationalBoundaryCellCycleModel::SetCycleDuration(double cycleduration)
{
    mCycleDuration = cycleduration;
}

double UniformG1GenerationalBoundaryCellCycleModel::GetCycleDuration() const
{
    return mCycleDuration;
}


bool UniformG1GenerationalBoundaryCellCycleModel::ReadyToDivide()
{
   assert(mpCell != NULL);
   //si Normal
   if (!mReadyToDivide)
   {
     if (!mpCell->HasCellProperty<CellLumen>())
     {
       UpdateCellCyclePhase();
       if ( (GetCurrentCellCyclePhase() != G_ZERO_PHASE) &&
            (GetAge() >=  GetMDuration() +  GetG1Duration() +  GetSDuration() + GetG2Duration()) )
       {
           mReadyToDivide = true;
       }
     }
   }


   //Si pour lumen
   if (!mReadyToDivide)
   {
     if (mpCell->HasCellProperty<CellEpi>() && mpCell->HasCellProperty<CellCore>())
     {
       double timeFromLastGen = mpCell->GetCellData()->GetItem("timeFromLastLumenGeneration");

       double xl = mpCell->GetCellData()->GetItem("vecPolaX");
       double yl = mpCell->GetCellData()->GetItem("vecPolaY");

       double norme = sqrt(xl*xl + yl*yl);
       if(norme >SimulationParameters::THRESHOLD_POLARISATION_EPI)
       {

         if(timeFromLastGen > SimulationParameters::TIME_BEETWEN_TWO_LUMEN_GENERATION)
         {
           if (GetAge() > SimulationParameters::AGE_DIV_LUMEN_MIN)
           {
             if(mpCell->GetCellData()->GetItem("lumenNearby") == 0)
             {
               //on ajoute un tag pour savoir que c'est pas une vraie division
               mpCell->GetCellData()->SetItem("hadALumenDivision",1);

               mReadyToDivide = true;
             }

           }
         }
       }
     }
   }


   return mReadyToDivide;
}




void UniformG1GenerationalBoundaryCellCycleModel::SetG1Duration()
{
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
    double a_thyr = 1.0035 ;
    double b_thyr = -0.747 ;
    double simulation_time = 48.0 ;
    double x_thyr = SimulationTime::Instance()->GetTime()/simulation_time ;

    assert(mpCell != nullptr);

    if (mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
    {
        mG1Duration = GetStemCellG1Duration() + 4*p_gen->ranf(); // U[14,18] for default parameters (mStemCellG1Duration) according to Meineke
    }
    else if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
    {
        if (mpCell->HasCellProperty<CellCore>())
        {
          // double CoreCycleDuration = 32;
          mG1Duration =  13/(a_thyr*exp(b_thyr*x_thyr)) ; // ICI : MODIFIER DUREE PHASE G1 POUR CELLULES DE COEUR, tu peux remplacer GetTransitCellG1Duration() par un nombre
          mSDuration =  13/(a_thyr*exp(b_thyr*x_thyr)) ;
          mG2Duration = 4/(a_thyr*exp(b_thyr*x_thyr)) ;
          mMDuration =  2/(a_thyr*exp(b_thyr*x_thyr)) ;
          mMaxTransitGenerations  = 6 ;
        }
        else if (mpCell->HasCellProperty<CellPeriph>())
        {
          // double BorderCycleDuration = 11 ;
          mG1Duration = 4/((a_thyr*exp(b_thyr*x_thyr))) ; // ICI : MODIFIER DUREE PHASE G1 POUR CELLULE PERIPH
          mSDuration =  4/(a_thyr*exp(b_thyr*x_thyr)) ;
          mG2Duration = 2/(a_thyr*exp(b_thyr*x_thyr)) ;
          mMDuration =  1/(a_thyr*exp(b_thyr*x_thyr)) ;
          mMaxTransitGenerations  = 6 ;
        }
        else if (mpCell->HasCellProperty<CellLumen>())
        {
          mG1Duration = 1/(a_thyr*exp(b_thyr*x_thyr)) ; // ICI : MODIFIER DUREE PHASE G1 POUR CELLULE LUMEN
          mSDuration =  1/(a_thyr*exp(b_thyr*x_thyr)) ;
          mG2Duration = 1/(a_thyr*exp(b_thyr*x_thyr)) ;
          mMDuration =  1/(a_thyr*exp(b_thyr*x_thyr)) ;
          mMaxTransitGenerations  = 20 ;
        }
        else
        {
          double G1Duration = (mCycleDuration*6/15)/(a_thyr*exp(b_thyr*x_thyr))    ;
          double SDuration = (mCycleDuration*6/15)/(a_thyr*exp(b_thyr*x_thyr)) ;
          double G2Duration = (mCycleDuration*2/15)/(a_thyr*exp(b_thyr*x_thyr))  ;
          double MDuration = (mCycleDuration*1/15)/(a_thyr*exp(b_thyr*x_thyr))  ;
          mG1Duration = G1Duration ;
          mSDuration =  SDuration ;
          mG2Duration = G2Duration ;
          mMDuration =  MDuration ;
          mMaxTransitGenerations = 6 ;
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
