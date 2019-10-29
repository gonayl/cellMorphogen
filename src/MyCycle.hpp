#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"

#include "SmartPointers.hpp"
#include "Exception.hpp"

#include "AbstractSimpleGenerationalCellCycleModel.hpp"

#include "CheckReadyToDivideAndPhaseIsUpdated.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "WildTypeCellMutationState.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "OffLatticeSimulation.hpp"
#include "StemCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

class MyCycle : public AbstractSimpleGenerationalCellCycleModel
{
private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleGenerationalCellCycleModel>(*this);
        RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
        archive & *p_gen;
        archive & p_gen;
    }

    void SetG1Duration()
    {
        /*assert(mpCell != NULL);

        double uniform_random_number = RandomNumberGenerator::Instance()->ranf();

        if (mpCell->GetCellProliferativeType()->IsType<StemCellProliferativeType>())
        {
            mG1Duration = -log(uniform_random_number)*GetStemCellG1Duration();
        }
        else if (mpCell->GetCellProliferativeType()->IsType<TransitCellProliferativeType>())
        {
            mG1Duration = -log(uniform_random_number)*GetTransitCellG1Duration();
        }
        else if (mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>())
        {
            mG1Duration = DBL_MAX;
        }
        else
        {
            NEVER_REACHED;
        }*/
        mG1Duration = 999999;
    }

public:

    MyCycle()
    {}

    AbstractCellCycleModel* CreateCellCycleModel()
    {
        MyCycle* p_model = new MyCycle();

        p_model->SetBirthTime(mBirthTime);
        p_model->SetMinimumGapDuration(mMinimumGapDuration);
        p_model->SetStemCellG1Duration(mStemCellG1Duration);
        p_model->SetTransitCellG1Duration(mTransitCellG1Duration);
        p_model->SetSDuration(mSDuration);
        p_model->SetG2Duration(mG2Duration);
        p_model->SetMDuration(mMDuration);
        p_model->SetGeneration(mGeneration);
        p_model->SetMaxTransitGenerations(mMaxTransitGenerations);

        return p_model;
    }
};
