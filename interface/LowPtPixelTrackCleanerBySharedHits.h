#ifndef _LowPtPixelTrackCleanerBySharedHits_h_
#define _LowPtPixelTrackCleanerBySharedHits_h_

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoPixelVertexing/PixelTrackFitting/interface/TracksWithHits.h"
#include "RecoPixelVertexing/PixelTrackFitting/interface/PixelTrackCleaner.h"

#include <utility>
#include <vector>

using namespace std;
using namespace pixeltrackfitting;

class LowPtPixelTrackCleanerBySharedHits : public PixelTrackCleaner
{
  public:
    LowPtPixelTrackCleanerBySharedHits(const edm::ParameterSet& ps);
    virtual ~LowPtPixelTrackCleanerBySharedHits();

    virtual TracksWithRecHits cleanTracks
     (const TracksWithRecHits & tracksWithRecHits);

  private:
    int getLayer(const DetId & id);
    bool hasCommonDetUnit (vector<const TrackingRecHit *> recHitsA,
                           vector<const TrackingRecHit *> recHitsB,
                           vector<DetId> detIds);
    bool hasCommonLayer (vector<const TrackingRecHit *> recHitsA,
                         vector<const TrackingRecHit *> recHitsB,
                         vector<int> detLayers);

};

#endif

