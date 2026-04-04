//
// Created by Nikhil Chandra Admal on 5/16/24.
//

#ifndef OILAB_GBPLASTICITY_H
#define OILAB_GBPLASTICITY_H

#include "LatticeCore.h"
#include "LatticeVector.h"
#include "Eigen/Dense"

namespace oILAB {

    template<int dim>
    class GbContinuum {
        using VectorDimD= LatticeCore<dim>::VectorDimD;
        using XuPairs = std::deque<std::pair<VectorDimD,VectorDimD>>;
    public:
        /*! Describes the normally flat GB. */
        const Eigen::Matrix<double,dim,dim-1> gbDomain;

        /*! This is a pair of XuPairs, where each XuPairs is a deque of nodal
         * position-displacement pairs that describe the faceted boundaries and
         * their displacements of respective grains.
         */
        const std::pair<XuPairs,XuPairs> xuPairsOfFacetedSurfaces;

        /*! Describes the outward unit normal to grainA in a flat GB. This is
         * used to infer which side of the faceted boundaries grainA occupies.
         */
        const VectorDimD normalGrainA;

        GbContinuum(const Eigen::Matrix<double, dim,dim-1>& domain,
                    const std::pair<XuPairs,XuPairs>& xuPairsOfFacetedSurfaces,
                    const VectorDimD& normalGrainA,
                    const bool& verbosity=false);

        /*! Returns the displacement at position \f$\mathbf{x}\f$ in the region occupied
         * by grain \f$i\f$. */
        VectorDimD displacement(const VectorDimD& x, const int& i) const;

        /*! Tells if the point \f$\mathbf{x}\f$ belongs to undeformed Grain A. */
        bool inGrainA(const VectorDimD& x) const;

        /*! Tells if the point \f$\mathbf{x}\f$ belongs to undeformed Grain B. */
        bool inGrainB(const VectorDimD& x) const;
    };

} // namespace oILAB

#include "GbContinuumImplementation.h"
#endif //OILAB_GBPLASTICITY_H
