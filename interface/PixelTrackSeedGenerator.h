#ifndef _PixelTrackSeedGenerator_h_
#define _PixelTrackSeedGenerator_h_

#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Common/interface/EDProduct.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "MagneticField/Engine/interface/MagneticField.h"

#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include <vector>
using namespace std;

class PixelTrackSeedGenerator
{
  public:
    PixelTrackSeedGenerator(const edm::EventSetup& es);
    ~PixelTrackSeedGenerator();
    TrajectorySeed seed(const reco::Track& recTrack);

  private:
    double getRadius(const TrackingRecHit& recHit);
    void sortRecHits(vector<pair<double,int> >& recHits);

    const MagneticField* theMagField;
    const Propagator* thePropagator;
    const TrackerGeometry* theTracker;
    const TransientTrackingRecHitBuilder* theTTRecHitBuilder;
};

#endif
