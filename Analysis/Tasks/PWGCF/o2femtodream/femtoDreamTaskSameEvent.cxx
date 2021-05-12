// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file femtoDreamTaskSameEvent.cxx
/// \brief Analysis task for particle pairing in the same event
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de

#include "o2femtodream/FemtoDerived.h"
#include "o2femtodream/FemtoDreamContainer.h"

#include "o2femtodream/FemtoDreamTrackSelection.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"

#include "TDatabasePDG.h"

using namespace o2;
using namespace o2::analysis::femtoDream;
using namespace o2::framework;

struct femtoDreamTaskSameEvent {

  FemtoDreamTrackSelection trackCuts;
  Configurable<FemtoDreamTrackSelection> CfgTrackSel{"FemtoDreamTrackSelection", {1, 0.5f, 4.05f, 0.8f, {{120}}, 0.83f, 70, 1, 0.1f, 0.2f, 3.f, 0.75f, o2::track::PID::Proton}, FemtoDreamTrackSelection::getCutHelp()};

  O2_DEFINE_CONFIGURABLE(CfgPDGCodePartOne, int, 2212, "PDG Code of particle one");
  O2_DEFINE_CONFIGURABLE(CfgPDGCodePartTwo, int, 2212, "PDG Code of particle two");

  /// Histograms
  FemtoDreamContainer sameEventCont;
  HistogramRegistry resultRegistry{"Correlations", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    sameEventCont.init(&resultRegistry);
    sameEventCont.setMasses(TDatabasePDG::Instance()->GetParticle(CfgPDGCodePartOne)->Mass(),
                            TDatabasePDG::Instance()->GetParticle(CfgPDGCodePartTwo)->Mass());

    trackCuts = CfgTrackSel;
    trackCuts.init(); //&qaRegistry);
  }

  void process(aod::FemtoDreamCollision const& col, aod::FemtoDreamParticles const& parts)
  {
    for (auto& [p1, p2] : combinations(parts, parts)) {
      sameEventCont.setPair(&resultRegistry, p1, p2);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<femtoDreamTaskSameEvent>(cfgc)};
  return workflow;
}
