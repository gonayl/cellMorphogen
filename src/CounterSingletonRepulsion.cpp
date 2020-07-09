#include <cassert>
#include <cmath>
#include "CounterSingletonRepulsion.hpp"

/** Pointer to the single instance */
CounterSingletonRepulsion* CounterSingletonRepulsion::mpInstance = nullptr;

CounterSingletonRepulsion* CounterSingletonRepulsion::Instance()
{
    if (mpInstance == nullptr)
    {
        mpInstance = new CounterSingletonRepulsion;
        std::atexit(Destroy);
    }
    return mpInstance;
}

unsigned CounterSingletonRepulsion::GetCount() const
{
    return mCount;
}

void CounterSingletonRepulsion::IncrementCounter()
{
    mCount++;
}

void CounterSingletonRepulsion::ResetCountToZero()
{
    mCount = 0;
}

void CounterSingletonRepulsion::Destroy()
{
    if (mpInstance)
    {
        delete mpInstance;
        mpInstance = nullptr;
    }
}

CounterSingletonRepulsion::CounterSingletonRepulsion()
    : mCount(0u)
{
    // Make sure there's only one instance - enforces correct serialization
    assert(mpInstance == nullptr);
}
