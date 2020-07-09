#ifndef COUNTERSINGLETONREPULSION_HPP_
#define COUNTERSINGLETONREPULSION_HPP_

#include <boost/shared_ptr.hpp>
#include "ChasteSerialization.hpp"
#include <boost/serialization/shared_ptr.hpp>
#include "SerializableSingleton.hpp"

/**
 * Class using the singleton pattern to keep count of some event in a globally
 * consistent time.
 */
class CounterSingletonRepulsion : public SerializableSingleton<CounterSingletonRepulsion>
{
public:

    /**
     * @return a pointer to the counter object.
     * The first time this is called the counter object is created.
     */
    static CounterSingletonRepulsion* Instance();

    /**
     * @return the count
     */
    unsigned GetCount() const;

    /**
     * Increment the counter by one.
     */
    void IncrementCounter();

    /**
     * Reset the counter to zero.
     */
    void ResetCountToZero();

    /**
     * Destroy the CounterSingletonRepulsion instance. The next call to Instance() will
     * create a new instance, on which ResetCountToZero() must be called again
     * to reset the counter.
     *
     * This method *must* be called before program exit, to avoid a memory
     * leak.
     */
    static void Destroy();

protected:

    /**
     * Default constructor. Sets the counter to zero.
     */
    CounterSingletonRepulsion();

private:

    /**
     * A pointer to the singleton instance of this class.
     */
    static CounterSingletonRepulsion* mpInstance;

    /**
     * Stores the count.
     */
    unsigned mCount;

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Serialization of a CounterSingletonRepulsion object must be done with care.
     * Do not serialize this singleton directly.  Instead, serialize
     * the object returned by GetSerializationWrapper().
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mCount;
    }
};

#endif /*COUNTERSINGLETONREPULSION_HPP_*/
