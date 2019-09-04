#ifndef CALIBRATIONERRORWRITER_HPP_
#define CALIBRATIONERRORWRITER_HPP_

#include "AbstractCellPopulationWriter.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * A class written using the visitor pattern for writing population elements from a cell population to file.
 *
 * The output file is called results.vizelements by default.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CalibrationErrorWriter : public AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>
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
    CalibrationErrorWriter();

    /**
     * Visit the MeshBasedCellPopulation and write the index of each Element.
     *
     * Outputs a line of space-separated values of the form:
     * ...[element index]...
     *
     * This line is appended to the output written by AbstractCellBasedWriter, which is a single
     * value [present simulation time], followed by a tab.
     *
     * @param pCellPopulation a pointer to the MeshBasedCellPopulation to visit.
     */
    virtual void Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write the data.
     *
     * This is an empty dummy function, since this class is not defined for use with a CaBasedCellPopulation.
     *
     * @param pCellPopulation a pointer to the CaBasedCellPopulation to visit.
     */
    virtual void Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write the data.
     *
     * This is an empty dummy function, since this class is not defined for use with a NodeBasedCellPopulation.
     *
     * @param pCellPopulation a pointer to the NodeBasedCellPopulation to visit.
     */
    virtual void Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     * Visit the PottsBasedCellPopulation and write data for each Element.
     *
     * Outputs a line of space-separated values of the form:
     * ...[num nodes in element] [node 0 index] [node 1 index] [node 2 index]...
     *
     * where [node 0 index] denotes the global index of the Node that is contained
     * in the Element with local index 0, and so on.
     *
     * This line is appended to the output written by AbstractCellBasedWriter, which is a single
     * value [present simulation time], followed by a tab.
     *
     * @param pCellPopulation a pointer to the PottsBasedCellPopulation to visit.
     */
    virtual void Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     * Visit the VertexBasedCellPopulation and write data for each Element.
     *
     * Outputs a line of space-separated values of the form:
     * ...[num nodes in element] [node 0 index] [node 1 index] [node 2 index]...
     *
     * where [node 0 index] denotes the global index of the Node that is contained
     * in the Element with local index 0, and so on.
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
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CalibrationErrorWriter)

#endif /*CALIBRATIONERRORWRITER_HPP_*/
