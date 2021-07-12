// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file FemtoDreamDetaDphiStar.h
/// \brief Definition of the FemtoDreamDetaDphiStar Container for math calculations math calculations of CPR quantities
/// \author Valentina Mantovani Sarti, TU München, valentina.mantovani-sarti@ph.tum.de

#ifndef ANALYSIS_TASKS_PWGCF_FEMTODREAM_FEMTODREAMDETADPHISTAR_H_
#define ANALYSIS_TASKS_PWGCF_FEMTODREAM_FEMTODREAMDETADPHISTAR_H_

#include "Math/Vector4D.h"
#include "Math/Boost.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include <numbers>

#include <iostream>

namespace o2::analysis
{
namespace femtoDream
{

/// \class FemtoDreamDetaDphiStar
/// \brief Container for math calculations of CPR quantities

class FemtoDreamDetaDphiStar
{
 public:
  static const float PiGreek = TMath::Pi();

  template <typename T>
  static void PhiAtRadiiTPC(const T& part, const float magfield, std::vector<float> tmpRadiiTPC, std::vector<float> tmpVec)
  {
    //ISSUE: We do not have the magnetic field in the Collisions femtotable;
    //since in the PairTask we have the col, we can just define in the process:
    //float magfield = col.Bfield();
    float phi0 = part->phi();
    float charge = (part->sign()); //so we will need a Join table with FemtodreamParticles and FEMTODEBUGPARTS
    float pt = part->pt();
    for (size_t i = 0; i < tmpRadiiTPC.size(); i++) {
      tmpVec.push_back(phi0 - std::asin(0.3 * charge * 0.1 * magfield * tmpRadiiTPC.at(i) * 0.01 / (2. * pt)));
    }
  }

  static float AveragePhiStar(std::vector<float> tmpVec1, std::vector<float> tmpVec2)
  {
    const int num = tmpVec1.size();
    float dPhiAvg = 0;
    for (size_t i = 0; i < num; i++) {
      float dphi = tmpVec1.at(i) - tmpVec2.at(i);
      if (dphi > PiGreek) {
        dphi += -2 * PiGreek;
      } else if (dphi < - PiGreek) {
        dphi += 2 * PiGreek;
      }
      dphi = TVector2::Phi_mpi_pi(dphi);
      dPhiAvg += dphi;
    }
    return (dPhiAvg / (float)num);
  }

  template <typename T>
  static float DeltaEta(const T& part1, const T& part2) {
    return (part1->eta() - part2->eta());
  }


 };

} /* namespace femtoDream */
} /* namespace o2::analysis */

#endif /* ANALYSIS_TASKS_PWGCF_FEMTODREAM_FEMTODREAMDETADPHISTAR_H_ */