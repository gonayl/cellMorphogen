#ifndef DIFFERENTIALADHESIONFORCE_HPP_
#define DIFFERENTIALADHESIONFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "NagaiHondaForce.hpp"

#include <iostream>

/**
 * A force class for use in vertex-based simulations, based on a model
 * model proposed by T. Nagai and H. Honda ("A dynamic cell model for
 * the formation of epithelial tissues", Philosophical Magazine Part B
 * 81:699-719) to include differential adhesion between normal and
 * labelled cells. To include differential adhesion we override the
 * GetAdhesionParameter() method.
 *
 * Each of the model parameter member variables are rescaled such that
 * mDampingConstantNormal takes the default value 1, whereas Nagai and
 * Honda (who denote the parameter by nu) take the value 0.01.
 */
template<unsigned DIM>
class DifferentialAdhesionForce  : public NagaiHondaForce<DIM>
{
private:

    /**
     * Cell-cell adhesion energy parameter for two endo cells.
     * Has units of kg (cell size at equilibrium rest length)^2 s^-2.
     * Takes the default value 1.0.
     */
    double mEndoEndoAdhesionEnergyParameter;

    /**
     * Cell-cell adhesion energy parameter for two lumen cells.
     * Has units of kg (cell size at equilibrium rest length)^2 s^-2.
     * Takes the default value 1.0.
     */
    double mLumenLumenAdhesionEnergyParameter;

    /**
     * Cell-cell adhesion energy parameter for two lumen cells.
     * Has units of kg (cell size at equilibrium rest length)^2 s^-2.
     * Takes the default value 1.0.
     */
    double mEpiEpiAdhesionEnergyParameter;

    /**
     * Cell-cell adhesion energy parameter for two epi cells.
     * Has units of kg (cell size at equilibrium rest length)^2 s^-2.
     * Takes the default value 1.0.
     */
    double mCoreCoreAdhesionEnergyParameter;

    /**
     * Cell-cell adhesion energy parameter for two epi cells.
     * Has units of kg (cell size at equilibrium rest length)^2 s^-2.
     * Takes the default value 1.0.
     */
    double mCorePeriphAdhesionEnergyParameter;

    /**
     * Cell-cell adhesion energy parameter for two epi cells.
     * Has units of kg (cell size at equilibrium rest length)^2 s^-2.
     * Takes the default value 1.0.
     */
    double mPeriphPeriphAdhesionEnergyParameter;

    /**
     * Cell-cell adhesion energy parameter for endothelial and epithelial cells.
     * Has has units of kg (cell size at equilibrium rest length)^2 s^-2.
     * Takes the default value 1.0.
     */
    double mEndoEpiAdhesionEnergyParameter;

    /**
     * Cell-cell adhesion energy parameter for endothelial and epithelial cells.
     * Has has units of kg (cell size at equilibrium rest length)^2 s^-2.
     * Takes the default value 1.0.
     */
    double mEndoPeriphAdhesionEnergyParameter;

    /**
     * Cell-cell adhesion energy parameter for lumen and epithelial cells.
     * Has has units of kg (cell size at equilibrium rest length)^2 s^-2.
     * Takes the default value 1.0.
     */
    double mEpiLumenAdhesionEnergyParameter;

    /**
     * Cell-cell adhesion energy parameter for lumen and endothelial cells.
     * Has has units of kg (cell size at equilibrium rest length)^2 s^-2.
     * Takes the default value 1.0.
     */
    double mEndoLumenAdhesionEnergyParameter;

    /**
     * Cell-cell adhesion energy parameter for lumen and endothelial cells.
     * Has has units of kg (cell size at equilibrium rest length)^2 s^-2.
     * Takes the default value 1.0.
     */
    double mStalkStalkAdhesionEnergyParameter;

    /**
     * Cell-cell adhesion energy parameter for lumen and endothelial cells.
     * Has has units of kg (cell size at equilibrium rest length)^2 s^-2.
     * Takes the default value 1.0.
     */
    double mStalkTipAdhesionEnergyParameter;

    /**
     * Cell-cell adhesion energy parameter for lumen and endothelial cells.
     * Has has units of kg (cell size at equilibrium rest length)^2 s^-2.
     * Takes the default value 1.0.
     */
    double mTipTipAdhesionEnergyParameter;

    /**
     * Cell-boundary adhesion energy parameter for endothelial cells.
     * Has units of kg (cell size at equilibrium rest length)^2 s^-2.
     * Takes the default value 1.0.
     */
    double mEndoBoundaryAdhesionEnergyParameter;

    /**
     * Cell-boundary adhesion energy parameter for lumen cells.
     * Has units of kg (cell size at equilibrium rest length)^2 s^-2.
     * Takes the default value 1.0.
     */
    double mLumenBoundaryAdhesionEnergyParameter;

    /**
     * Cell-boundary adhesion energy parameter for epithelial cells.
     * Has units of kg (cell size at equilibrium rest length)^2 s^-2.
     * Takes the default value 1.0.
     */
    double mEpiBoundaryAdhesionEnergyParameter;

    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<NagaiHondaForce<DIM> >(*this);
        archive & mEndoEndoAdhesionEnergyParameter;
        archive & mLumenLumenAdhesionEnergyParameter;
        archive & mEpiEpiAdhesionEnergyParameter;
        archive & mCoreCoreAdhesionEnergyParameter;
        archive & mCorePeriphAdhesionEnergyParameter;
        archive & mPeriphPeriphAdhesionEnergyParameter;
        archive & mEndoEpiAdhesionEnergyParameter;
        archive & mEndoPeriphAdhesionEnergyParameter;
        archive & mEpiLumenAdhesionEnergyParameter;
        archive & mEndoLumenAdhesionEnergyParameter;
        archive & mStalkStalkAdhesionEnergyParameter;
        archive & mStalkTipAdhesionEnergyParameter;
        archive & mTipTipAdhesionEnergyParameter;
        archive & mEndoBoundaryAdhesionEnergyParameter;
        archive & mLumenBoundaryAdhesionEnergyParameter;
        archive & mEpiBoundaryAdhesionEnergyParameter;
    }

public:

    /**
     * Constructor.
     */
    DifferentialAdhesionForce();

    /**
     * Destructor.
     */
    virtual ~DifferentialAdhesionForce();

    /**
     * Overridden GetAdhesionParameter() method.
     *
     * Get the adhesion parameter for the edge between two given nodes. Depends
     * on the type of cells attached to the elements.
     *
     * @param pNodeA one node
     * @param pNodeB the other node
     * @param rVertexCellPopulation reference to the cell population
     *
     * @return the adhesion parameter for this edge.
     */

    virtual double GetAdhesionParameter(Node<DIM>* pNodeA, Node<DIM>* pNodeB, VertexBasedCellPopulation<DIM>& rVertexCellPopulation);


    double GetEndoEndoAdhesionEnergyParameter();

    double GetLumenLumenAdhesionEnergyParameter();

    double GetEpiEpiAdhesionEnergyParameter();

    double GetCoreCoreAdhesionEnergyParameter();

    double GetCorePeriphAdhesionEnergyParameter();

    double GetPeriphPeriphAdhesionEnergyParameter();

    double GetEndoEpiAdhesionEnergyParameter();

    double GetEndoPeriphAdhesionEnergyParameter();

    double GetEpiLumenAdhesionEnergyParameter();

    double GetEndoLumenAdhesionEnergyParameter();

    double GetStalkStalkAdhesionEnergyParameter();

    double GetStalkTipAdhesionEnergyParameter();

    double GetTipTipAdhesionEnergyParameter();

    double GetEndoBoundaryAdhesionEnergyParameter();

    double GetLumenBoundaryAdhesionEnergyParameter();

    double GetEpiBoundaryAdhesionEnergyParameter();


    void SetEndoEndoAdhesionEnergyParameter(double endoEndoAdhesionEnergyParameter);

    void SetLumenLumenAdhesionEnergyParameter(double lumenLumenAdhesionEnergyParameter);

    void SetEpiEpiAdhesionEnergyParameter(double epiEpiAdhesionEnergyParameter);

    void SetCoreCoreAdhesionEnergyParameter(double coreCoreAdhesionEnergyParameter);

    void SetCorePeriphAdhesionEnergyParameter(double corePeriphAdhesionEnergyParameter);

    void SetPeriphPeriphAdhesionEnergyParameter(double periphPeriphAdhesionEnergyParameter);

    void SetEndoEpiAdhesionEnergyParameter(double endoEpiAdhesionEnergyParameter);

    void SetEndoPeriphAdhesionEnergyParameter(double endoPeriphAdhesionEnergyParameter);

    void SetEpiLumenAdhesionEnergyParameter(double epiLumenAdhesionEnergyParameter);

    void SetEndoLumenAdhesionEnergyParameter(double endoLumenAdhesionEnergyParameter);

    void SetStalkStalkAdhesionEnergyParameter(double stalkStalkAdhesionEnergyParameter);

    void SetStalkTipAdhesionEnergyParameter(double stalkTipAdhesionEnergyParameter);

    void SetTipTipAdhesionEnergyParameter(double tipTipAdhesionEnergyParameter);

    void SetEndoBoundaryAdhesionEnergyParameter(double endoBoundaryAdhesionEnergyParameter);

    void SetLumenBoundaryAdhesionEnergyParameter(double lumenBoundaryAdhesionEnergyParameter);

    void SetEpiBoundaryAdhesionEnergyParameter(double epiBoundaryAdhesionEnergyParameter);


    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(DifferentialAdhesionForce)

#endif /*DIFFERENTIALADHESIONFORCE_HPP_*/
