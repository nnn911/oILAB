//
// Created by Nikhil Chandra Admal on 5/27/24.
//

#ifndef OILAB_GBMESOSTATEIMPLEMENTATION_H
#define OILAB_GBMESOSTATEIMPLEMENTATION_H

#include "../IO/Lammps.h"
#include <iostream>
#include <set>
#include "GbMesoState.h"
#include "OrderedTuplet.h"

namespace oILAB {
template <int dim>
GbMesoState<dim>::GbMesoState(
    const Gb<dim> &gb,
    const ReciprocalLatticeVector<dim> &axis,
    const std::deque<std::pair<LatticeVector<dim>, VectorDimD>>& engagedTsPairs,
    const std::vector<LatticeVector<dim>> &mesoStateCslVectors) try:
    /*init*/ GbContinuum<dim>(getMesoStateGbDomain(mesoStateCslVectors),
                              getFacetedSurfaces(gb,mesoStateCslVectors,engagedTsPairs),
                              gb.nA.cartesian(),
                              true),
    /*init*/ gb(gb),
    /*init*/ axis(axis),
    /*init*/ mesoStateCslVectors(mesoStateCslVectors),
    /*init*/ engagedTsPairs(engagedTsPairs) {}
    catch(std::runtime_error& e)
    {
        throw;
    }


    template<int dim>
    Eigen::Matrix<double, dim,dim-1> GbMesoState<dim>::getMesoStateGbDomain(const std::vector<LatticeVector<dim>>& mesoStateCslVectors)
    {
        Eigen::Matrix<double, dim,dim-1> mesoStateGbDomain;
        for(int i=1; i<dim; ++i)
            mesoStateGbDomain.col(i-1)= mesoStateCslVectors[i].cartesian();
        return mesoStateGbDomain;
    }

    template<int dim>
    std::pair<typename GbMesoState<dim>::XuPairs, typename GbMesoState<dim>::XuPairs>
    GbMesoState<dim>::getFacetedSurfaces(const Gb<dim>& gb,
                                 const std::vector<LatticeVector<dim>>& mesoStateCslVectors,
                                 const std::deque<std::pair<LatticeVector<dim>,VectorDimD>>& engagedTsPairs)
    {
        XuPairs xuPairsA, xuPairsB;

        std::vector<LatticeVector<dim>> bicrystalBoxVectors(mesoStateCslVectors);
        bicrystalBoxVectors[0]= 2*mesoStateCslVectors[0];
        VectorDimD shift;
        shift << -0.5-FLT_EPSILON,-FLT_EPSILON,-FLT_EPSILON;

        std::set<OrderedTuplet<dim>> xAIntegerCoordsSet, xBIntegerCoordsSet;
        for(const auto& [t,s] : engagedTsPairs)
        {
            VectorDimD u= t.cartesian()/2;
            VectorDimD xA= s-u;
            VectorDimD xB= s+u;

            // modulo tempx w.r.t the bicrystal box
            LatticeVector<dim>::modulo(xA,bicrystalBoxVectors,shift);
            LatticeVector<dim>::modulo(xB,bicrystalBoxVectors,shift);

            OrderedTuplet<dim> xAIntegerCoords, xBIntegerCoords;
            try {
                xAIntegerCoords<< gb.bc.A.latticeVector(xA);
            }
            catch(std::runtime_error& e) {
                std::cout << e.what() << std::endl;
                std::cout << "xA not a lattice vector of A" << std::endl;
                exit(0);
            }
            try {
                xBIntegerCoords << gb.bc.B.latticeVector(xB);
            }
            catch(std::runtime_error& e) {
                std::cout << e.what() << std::endl;
                std::cout << "xB not a lattice vector of B" << std::endl;
                exit(0);
            }
            const bool insertedA= xAIntegerCoordsSet.insert(xAIntegerCoords).second;
            const bool insertedB= xBIntegerCoordsSet.insert(xBIntegerCoords).second;
            if(!insertedA || !insertedB)
                throw std::runtime_error("Clash in constraints.");

            xuPairsA.emplace_back(xA,u);
            xuPairsB.emplace_back(xB,-u);
        }
        return std::make_pair(xuPairsA,xuPairsB);
    }

    /*-------------------------------------*/
    template<int dim>
    std::tuple<double,double> GbMesoState<dim>::densityEnergy(const std::string& lmpLocation,
                                                             const std::string& potentialName) const
    {
        box("temp" + std::to_string(omp_get_thread_num()));
        std::pair<double,double> densityEnergyPair= energy(lmpLocation,
                                                           "temp" + std::to_string(omp_get_thread_num()) + "_reference1.txt",
                                                           potentialName);


        return {densityEnergyPair.first,densityEnergyPair.second};

    }


    template<int dim>
    typename std::enable_if<dim==3,void>::type
    GbMesoState<dim>::box(const std::string& name) const
    {
        const auto& config= gb.bc.box(mesoStateCslVectors,0);
        std::vector<LatticeVector<3>> boxVectors;
        boxVectors.push_back(this->mesoStateCslVectors[0]);
        boxVectors.push_back(this->mesoStateCslVectors[1]);
        boxVectors.push_back(this->mesoStateCslVectors[2]);

        std::vector<VectorDimD> referenceConfigA, deformedConfigA;
        std::vector<VectorDimD> referenceConfigB, deformedConfigB;

        for (const auto &latticeVector: config) {
            VectorDimD x;
            if (&(latticeVector.lattice) == &(gb.bc.A) && this->inGrainA(latticeVector.cartesian())) {
                x= latticeVector.cartesian() + this->displacement(latticeVector.cartesian(),1);
                referenceConfigA.push_back(latticeVector.cartesian());
                deformedConfigA.push_back(x);
            }
            else if (&(latticeVector.lattice) == &(gb.bc.B) && this->inGrainB(latticeVector.cartesian())) {
                x= latticeVector.cartesian() + this->displacement(latticeVector.cartesian(),2);
                referenceConfigB.push_back(latticeVector.cartesian());
                deformedConfigB.push_back(x);
            }
        }

        int nAtoms= referenceConfigA.size()+referenceConfigB.size();

        std::string referenceFile= name + "_reference0.txt";
        std::string deformedFile= name + "_reference1.txt";

        std::ofstream reference, deformed;
        reference.open(referenceFile);
        deformed.open(deformedFile);
        if (!reference || !deformed) std::cerr << "Unable to open files";
        reference << nAtoms << std::endl; deformed << nAtoms << std::endl;
        reference << "Lattice=\" "; deformed << "Lattice=\" ";

        reference << std::setprecision(15) << (2*boxVectors[0].cartesian()).transpose() << " ";
        deformed << std::setprecision(15) << (2*boxVectors[0].cartesian()).transpose() << " ";
        reference << std::setprecision(15) << (boxVectors[1].cartesian()).transpose() << " ";
        deformed << std::setprecision(15) << (boxVectors[1].cartesian()).transpose() << " ";
        reference << std::setprecision(15) << (boxVectors[2].cartesian()).transpose();
        deformed << std::setprecision(15) << (boxVectors[2].cartesian()).transpose();
        reference << "\" Properties=atom_types:I:1:pos:R:3:radius:R:1 PBC=\"F T T\" origin=\" "; deformed << "\" Properties=atom_types:I:1:pos:R:3:radius:R:1 PBC=\" F T T\" origin=\" ";
        reference << std::setprecision(15) << (-1 * boxVectors[0].cartesian()).transpose() << "\"" << std::endl;
        deformed << std::setprecision(15) << (-1 * boxVectors[0].cartesian()).transpose() << "\"" << std::endl;

        for(const auto& position : referenceConfigA)
            reference << 1 << " " << std::setprecision(15) << position.transpose() << "  " << 0.05 << std::endl;
        for(const auto& position : referenceConfigB)
            reference << 2 << " " << std::setprecision(15) << position.transpose() << "  " << 0.05 << std::endl;
        for(const auto& position : deformedConfigA)
            deformed << 1 << " " << std::setprecision(15) << position.transpose() << "  " << 0.05 << std::endl;
        for(const auto& position : deformedConfigB)
            deformed << 2 << " " << std::setprecision(15) << position.transpose() << "  " << 0.05 << std::endl;

        reference.close();
        deformed.close();
    }

 } // namespace oILAB

#endif
