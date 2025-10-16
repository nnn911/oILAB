//
// Created by Nikhil Chandra Admal on 8/14/24.
//

#ifndef OILAB_EVOLUTIONALGORITHMIMPLEMENTATION_H
#define OILAB_EVOLUTIONALGORITHMIMPLEMENTATION_H

#include "../Utilities/randomInteger.h"
#include "../Lattices/OrderedTuplet.h"
#include "CanonicalTP.h"
#include "LandauWangTP.h"
#include "../Lattices/GbMesoState.h"

namespace oILAB {
template <typename StateType, typename SystemType,
          typename TransitionProbabilityType>
EvolutionAlgorithm<StateType, SystemType,
                   TransitionProbabilityType>::EvolutionAlgorithm()
    : transitionProbability(static_cast<TransitionProbabilityType &>(*this)) {}

template <typename StateType, typename SystemType,
          typename TransitionProbabilityType>
bool EvolutionAlgorithm<StateType, SystemType, TransitionProbabilityType>::
    acceptMove(
        const std::pair<StateType, SystemType> &proposedStateSystem,
        const std::pair<StateType, SystemType> &currentStateSystem) const {
  double probability = transitionProbability.probability(proposedStateSystem,
                                                         currentStateSystem);
  if (random<double>(0.0, 1.0) <= probability)
    return true;
  else
    return false;
    }

    } // namespace oILAB
#endif