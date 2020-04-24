#ifndef SURFACELUMENWRITER_HPP_
#define SURFACELUMENWRITER_HPP_

#include "AbstractCellPopulationWriter.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * A class written using the visitor pattern for writing cell population volume data to file.
 * Used by MeshBasedCellPopulation.
 *
 * The output file is called SurfaceLumens.dat by default.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class SurfaceLumenWriter : public AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>
{
private:
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:

    /**
     * Default constructor.
     */
    SurfaceLumenWriter();

    /**
     * Visit the population and write the areas (or volume, in 3D) occupied by the
     * entire cell population and by apoptotic cells. Any ghost nodes present in
     * the population are neglected from these areas.
     *
     * Outputs a line of space-separated values of the form:
     * [total area] [apoptotic area]
     *
     * This line is appended to the output written by AbstractCellBasedWriter, which is a single
     * value [present simulation time], followed by a tab.
     *
     * @param pCellPopulation a pointer to the MeshBasedCellPopulation to visit.
     */
    virtual void Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write the areas (or volume, in 3D) occupied by the
     * entire cell population and by apoptotic cells. Any ghost nodes present in
     * the population are neglected from these areas.
     *
     * Outputs a line of space-separated values of the form:
     * [total area] [apoptotic area]
     *
     * This line is appended to the output written by AbstractCellBasedWriter, which is a single
     * value [present simulation time], followed by a tab.
     *
     * @param pCellPopulation a pointer to the CaBasedCellPopulation to visit.
     */
    virtual void Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write the areas (or volume, in 3D) occupied by the
     * entire cell population and by apoptotic cells. Any ghost nodes present in
     * the population are neglected from these areas.
     *
     * Outputs a line of space-separated values of the form:
     * [total area] [apoptotic area]
     *
     * This line is appended to the output written by AbstractCellBasedWriter, which is a single
     * value [present simulation time], followed by a tab.
     *
     * @param pCellPopulation a pointer to the NodeBasedCellPopulation to visit.
     */
    virtual void Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write the areas (or volume, in 3D) occupied by the
     * entire cell population and by apoptotic cells. Any ghost nodes present in
     * the population are neglected from these areas.
     *
     * Outputs a line of space-separated values of the form:
     * [total area] [apoptotic area]
     *
     * This line is appended to the output written by AbstractCellBasedWriter, which is a single
     * value [present simulation time], followed by a tab.
     *
     * @param pCellPopulation a pointer to the PottsBasedCellPopulation to visit.
     */
    virtual void Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write the areas (or volume, in 3D) occupied by the
     * entire cell population and by apoptotic cells. Any ghost nodes present in
     * the population are neglected from these areas.
     *
     * Outputs a line of space-separated values of the form:
     * [total area] [apoptotic area]
     *
     * This line is appended to the output written by AbstractCellBasedWriter, which is a single
     * value [present simulation time], followed by a tab.
     *
     * @param pCellPopulation a pointer to the VertexBasedCellPopulation to visit.
     */
    virtual void Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation);
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_ALL_DIMS(SurfaceLumenWriter)

#endif /*SurfaceLumenWRITER_HPP_*/
