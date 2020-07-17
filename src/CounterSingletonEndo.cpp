#include <cassert>
#include <cmath>
#include "CounterSingletonEndo.hpp"

/** Pointer to the single instance */
CounterSingletonEndo* CounterSingletonEndo::mpInstance = nullptr;

CounterSingletonEndo* CounterSingletonEndo::Instance()
{
    if (mpInstance == nullptr)
    {
        mpInstance = new CounterSingletonEndo;
        std::atexit(Destroy);
    }
    return mpInstance;
}

unsigned CounterSingletonEndo::GetCount() const
{
    return mCount;
}

void CounterSingletonEndo::IncrementCounter()
{
    mCount++;
}

void CounterSingletonEndo::ResetCountToZero()
{
    mCount = 0;
}

void CounterSingletonEndo::Destroy()
{
    if (mpInstance)
    {
        delete mpInstance;
        mpInstance = nullptr;
    }
}

CounterSingletonEndo::CounterSingletonEndo()
    : mCount(0u)
{
    // Make sure there's only one instance - enforces correct serialization
    assert(mpInstance == nullptr);
}
