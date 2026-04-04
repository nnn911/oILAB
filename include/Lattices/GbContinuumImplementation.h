//
// Created by Nikhil Chandra Admal on 5/27/24.
//

#ifndef OILAB_GBCONTINUUMIMPLEMENTATION_H
#define OILAB_GBCONTINUUMIMPLEMENTATION_H
#include "GbMesoState.h"

namespace oILAB {

    template <int dim>
    GbContinuum<dim>::GbContinuum(const Eigen::Matrix<double, dim, dim - 1> &domain,
                                  const std::pair<XuPairs,XuPairs>& xuPairsOfFacetedSurfaces,
                                  const VectorDimD& normalGrainA,
                                  const bool &verbosity) :
    /*init*/ gbDomain(domain),
    /*init*/ xuPairsOfFacetedSurfaces(xuPairsOfFacetedSurfaces),
    /*init*/ normalGrainA(normalGrainA){
        // you will have to assert that xuPairsOfFacetedSurfaces is consistent: if you
        // displace the nodes of the two surfaces they must match.
        if (verbosity) {
            std::cout << "-------------------------------------------------------------"
                         "-----------------"
                      << std::endl;
            std::cout << std::endl;
        }
    }

    template<int dim>
    GbContinuum<dim>::VectorDimD GbContinuum<dim>::displacement(const VectorDimD& x,
                                                                const int& i) const
    {
        VectorDimD u;
        u.setZero();

        /*
        We also need to be careful about the following: sometimes coincident points other
        than those engaged are mistakenly engaged, i.e. form coincidence, after deformation.
        In such cases we delete those unwanted coincidences. If done successfully, the number of
        coincidence points/density should exactly match with that expected from the mesostate signature
        */
        return u;
    }

    template<int dim>
    bool GbContinuum<dim>::inGrainA(const VectorDimD &x) const {
        // check if x belongs to the grain 1, i.e. the appropriate side of the first
        // faceted surface. You will have to use the variable normalGrainA. It is sometimes
        // possible that x belongs to both grainA region and grainB region
        return false;
    }
    template<int dim>
    bool GbContinuum<dim>::inGrainB(const VectorDimD &x) const {
        // check if x belongs to the grain 2, i.e. the appropriate side of the second
        // faceted surface. You will have to use the variable normalGrainA. It is sometimes
        // possible that x belongs to both grainA region and grainB region.
        return false;
    }
} // namespace oILAB

#endif