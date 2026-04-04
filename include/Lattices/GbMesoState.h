//
// Created by Nikhil Chandra Admal on 5/26/24.
//

#ifndef OILAB_GBMESOSTATE_H
#define OILAB_GBMESOSTATE_H

#include "Gb.h"
#include "GbContinuum.h"
#include "LatticeCore.h"
#include "../Math/PeriodicFunction.h"

namespace oILAB {

/*! Class template that defines a GB mesostate.
 *
 */
template <int dim> class GbMesoState : public GbContinuum<dim> {
    using VectorDimD = LatticeCore<dim>::VectorDimD;
    using XuPairs = std::deque<std::pair<VectorDimD,VectorDimD>>;

  /*!
   * \brief Returns the cartesian coordinates of the CSL vectors that define a
   * mesostate's GB.
   * @param mesoStateCslVectors - a vector (size = \p dim) of CSL vectors that
   * define the box of the mesostate.
   * @return Cartesian coordinates of the \p dim-1 grain boundary CSL vectors
   */
  // ensure that the input is of the right dimension
  static Eigen::Matrix<double, dim, dim - 1> getMesoStateGbDomain(const std::vector<LatticeVector<dim>> &mesoStateCslVectors);

  /*!
   * \brief Returns the nodes \f$\textbf x\f$ and their displacements \f$\textbf u\f$ of the faceted boundaries of the
   * grains that form the GB.
   * @param gb - grain boundary
   * @param mesoStateCslVectors - the box vectors of the mesostate
   * @param engagedTsPairs - a deque of engaged (t,s) pairs (translation vector \p t, shift vector \p s).
   * \p t is a DSCL vector while \p s is expressed in Cartesian coordinates
   * @return A pair of deques for the boundaries of the two undeformed grains that are deformed and glued to form the GB.
   * Each deque contains `(x,u)` pairs.
   */
  static std::pair<XuPairs,XuPairs> getFacetedSurfaces(const Gb<dim> &gb,
      const std::vector<LatticeVector<dim>> &mesoStateCslVectors,
      const std::deque<std::pair<LatticeVector<dim>, VectorDimD>> &engagedTsPairs);

public:
  /*!
   * Grain boundary
   */
  const Gb<dim> &gb;

  /*!
   * Grain boundary tilt axis
   */
  const ReciprocalLatticeVector<dim> &axis;

  /*!
   * A vector (size = \p dim) of CSL vectors that define the box of the
   * mesostate. The second and third vectors should be parallel to the grain
   * boundary, while the third vector should be out of the grain boundary plane.
   */
  const std::vector<LatticeVector<dim>> &mesoStateCslVectors;

  /*!
   * @param engagedTsPairs a deque of pairs (translation vector \f$\textbf t\f$, shift
   * vector \f$\textbf s\f$) that defines a mesostate. Translating lattice
    * \f$\mathcal A\f$ by \f$\textbf t/2\f$ and lattice
     * \f$\mathcal B\f$ by \f$-\textbf t/2\f$ results in a CSL shift of \f$\textbf
   * s\f$.
   */
  const std::deque<std::pair<LatticeVector<dim>, VectorDimD>> engagedTsPairs;

  explicit GbMesoState(
      const Gb<dim> &gb,
      const ReciprocalLatticeVector<dim> &axis,
      const std::deque<std::pair<LatticeVector<dim>, VectorDimD>>& engagedTsPairs,
      const std::vector<LatticeVector<dim>> &mesoStateCslVectors);

  /*!
   * \brief Calculate the energy of a mesostate using lammps
   * @return Energy of the mesostate
   */
  // std::pair<double,double> densityEnergy() const;
  std::tuple<double, double>
  densityEnergy(const std::string &lmpLocation,
                const std::string &potentialName) const;

  /*! This function outputs/prints a grain boundary mesostate
   * @param filename name of the file to be written to
   */
  typename std::enable_if<dim == 3, void>::type
  box(const std::string &filename) const;
    };
    } // namespace oILAB

#include "GbMesoStateImplementation.h"
#endif //OILAB_GBMESOSTATE_H
