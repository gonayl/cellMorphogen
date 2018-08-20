#include <cassert>
#include <cmath>
#include "CounterSingleton.hpp"

/** Pointer to the single instance */
CounterSingleton* CounterSingleton::mpInstance = nullptr;

CounterSingleton* CounterSingleton::Instance()
{
    if (mpInstance == nullptr)
    {
        mpInstance = new CounterSingleton;
        std::atexit(Destroy);
    }
    return mpInstance;
}

unsigned CounterSingleton::GetCount() const
{
    return mCount;
}

void CounterSingleton::IncrementCounter()
{
    mCount++;
}

void CounterSingleton::ResetCountToZero()
{
    mCount = 0;
}

void CounterSingleton::Destroy()
{
    if (mpInstance)
    {
        delete mpInstance;
        mpInstance = nullptr;
    }
}

CounterSingleton::CounterSingleton()
    : mCount(0u)
{
    // Make sure there's only one instance - enforces correct serialization
    assert(mpInstance == nullptr);
}
