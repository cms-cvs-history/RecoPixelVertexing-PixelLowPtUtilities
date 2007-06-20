#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "RecoPixelVertexing/PixelLowPtUtilities/interface/EventPlotter.h"

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByHits.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"

#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"

#include "TROOT.h"
#include "TFile.h"
#include "TNtuple.h"

#include <fstream>
using namespace std;

/*****************************************************************************/
class LowPtTrackAnalyzer : public edm::EDAnalyzer
{
 public:
   explicit LowPtTrackAnalyzer(const edm::ParameterSet& pset);
   ~LowPtTrackAnalyzer();
   virtual void beginJob(const edm::EventSetup& es);
   virtual void analyze(const edm::Event& ev, const edm::EventSetup& es);
   virtual void endJob();

 private:
   int layerFromDetid(const DetId& detId);
   int getNumberOfSimHits(const TrackingParticle& simTrack);
   int getDetLayerId(const PSimHit& simHit);
   int getNumberOfPixelHits(const TrackingParticle& simTrack);
   void checkSimTracks(edm::Handle<TrackingParticleCollection>& simCollection,
                       reco::SimToRecoCollection& q);
   void checkRecTracks(edm::Handle<reco::TrackCollection>& recCollection,
                       reco::RecoToSimCollection& p);

   const TrackerGeometry* theTracker;
   const TrackAssociatorByHits* theAssociatorByHits;

   TFile * resultFile; 
   TNtuple * trackSim;
   TNtuple * trackRec;

   bool plotEvent;
   bool zipFiles;
};

/*****************************************************************************/
LowPtTrackAnalyzer::LowPtTrackAnalyzer(const edm::ParameterSet& pset)
{
  plotEvent = pset.getParameter<bool>("plotEvent");
  zipFiles  = pset.getParameter<bool>("zipFiles");
}

/*****************************************************************************/
LowPtTrackAnalyzer::~LowPtTrackAnalyzer()
{
}

/*****************************************************************************/
void LowPtTrackAnalyzer::beginJob(const edm::EventSetup& es)
{
  // Get tracker geometry
  edm::ESHandle<TrackerGeometry> tracker;
  es.get<TrackerDigiGeometryRecord>().get(tracker);
  theTracker = tracker.product();
 
  // Get aassociator
  edm::ESHandle<TrackAssociatorBase> theHitsAssociator;
  es.get<TrackAssociatorRecord>().get("TrackAssociatorByHits",
                                    theHitsAssociator);
  theAssociatorByHits =
   (const TrackAssociatorByHits*)theHitsAssociator.product();
 
  // Root
  resultFile = new TFile("result.root","recreate");
  resultFile->cd();
   
  trackSim = new TNtuple("trackSim","trackSim",
    "id:eta:pt:d0:nhit:nrec");
  trackRec = new TNtuple("trackRec","trackRec",
    "qr:etar:ptr:d0r:nsim:ids:parids:etas:pts:d0s");

}

/*****************************************************************************/
void LowPtTrackAnalyzer::endJob()
{
  resultFile->cd();
  trackSim->Write();
  trackRec->Write();
  resultFile->Close();
}

/*****************************************************************************/
int LowPtTrackAnalyzer::layerFromDetid(const DetId& detId)
{
  int layerNumber=0;
  unsigned int subdetId = static_cast<unsigned int>(detId.subdetId());
  if ( subdetId == StripSubdetector::TIB)
    {
      TIBDetId tibid(detId.rawId());
      layerNumber = tibid.layer();
    }
  else if ( subdetId ==  StripSubdetector::TOB )
    {
      TOBDetId tobid(detId.rawId());
      layerNumber = tobid.layer();
    }
  else if ( subdetId ==  StripSubdetector::TID)
    {
      TIDDetId tidid(detId.rawId());
      layerNumber = tidid.wheel();
    }
  else if ( subdetId ==  StripSubdetector::TEC )
    {
      TECDetId tecid(detId.rawId());
      layerNumber = tecid.wheel();
    }
  else if ( subdetId ==  PixelSubdetector::PixelBarrel )
    {
      PXBDetId pxbid(detId.rawId());
      layerNumber = pxbid.layer();
    }
  else if ( subdetId ==  PixelSubdetector::PixelEndcap )
    {
      PXFDetId pxfid(detId.rawId());
      layerNumber = pxfid.disk();
    }
  else
    edm::LogVerbatim("TrackValidator") << "Unknown subdetid: " << subdetId;

  return layerNumber;
}

/*****************************************************************************/
int LowPtTrackAnalyzer::getNumberOfSimHits(const TrackingParticle& simTrack)
{
  int oldlay = 0; int newlay = 0;
  int olddet = 0; int newdet = 0;

  int nhit = 0;

  for(std::vector<PSimHit>::const_iterator
       simHit = simTrack.pSimHit_begin();
       simHit!= simTrack.pSimHit_end(); simHit++)
  {
    const DetId detId = DetId(simHit->detUnitId());
    oldlay = newlay; newlay = layerFromDetid(detId);
    olddet = newdet; newdet = detId.subdetId();
    if(oldlay != newlay || (oldlay == newlay && olddet != newdet) ) nhit++;
  }

  return nhit;
}

/*****************************************************************************/
int LowPtTrackAnalyzer::getDetLayerId(const PSimHit& simHit)
{
  int layerId;

  DetId id = DetId(simHit.detUnitId());
  LocalPoint lpos = simHit.localPosition();
  GlobalPoint gpos = theTracker->idToDetUnit(id)->toGlobal(lpos);

  if(theTracker->idToDetUnit(id)->subDetector() ==
      GeomDetEnumerators::PixelBarrel)
  { // barrel
    if(gpos.perp2() < 6 * 6) layerId = 0;
    else
    {
      if(gpos.perp2() < 9 * 9) layerId = 1;
                          else layerId = 2;
    }
  }
  else
  { // endcap
    if(fabsf(gpos.z()) < 40) layerId = 3;
                        else layerId = 4;
  }

  return layerId;
}

/*****************************************************************************/
int LowPtTrackAnalyzer::getNumberOfPixelHits(const TrackingParticle&
simTrack)
{
  // How many pixel hits?
  const int nLayers = 5;
  vector<bool> filled(nLayers,false);

  int numberOfPixelHits = 0;

  for(std::vector<PSimHit>::const_iterator simHit = simTrack.pSimHit_begin();
                                           simHit!= simTrack.pSimHit_end();
                                           simHit++)
  {
    DetId id = DetId(simHit->detUnitId());

    if(theTracker->idToDetUnit(id)->subDetector() ==
       GeomDetEnumerators::PixelBarrel ||
       theTracker->idToDetUnit(id)->subDetector() ==
       GeomDetEnumerators::PixelEndcap)
    {
      filled[getDetLayerId(*simHit)] = true;
      numberOfPixelHits++;
    }
  }
  
  // Count the number of filled pixel layers
  int fLayers = 0;
  for(int i=0; i<nLayers; i++)
    if(filled[i] == true) fLayers++;
  
  return numberOfPixelHits * (fLayers >= 3 ? 1 : -1);
}

/*****************************************************************************/
void LowPtTrackAnalyzer::checkSimTracks
  (edm::Handle<TrackingParticleCollection>& simCollection,
   reco::SimToRecoCollection& q)
{
  for(TrackingParticleCollection::size_type i=0;
               i < simCollection.product()->size(); ++i)
  {
    const TrackingParticleRef simTrack(simCollection, i);

    vector<float> result;

    // sim
    result.push_back(simTrack->pdgId());
    result.push_back(simTrack->eta());
    result.push_back(simTrack->pt());
    result.push_back(simTrack->vertex().Rho());
    result.push_back(getNumberOfPixelHits(*simTrack));

    // rec
    int nRec = 0;
    try
    {
      vector<pair<reco::TrackRef, double> > recTracks = q[simTrack];

      for(vector<pair<reco::TrackRef,double> >::const_iterator
            it = recTracks.begin(); it != recTracks.end(); ++it)
      {
        reco::TrackRef recTrack = it->first;
        int nShared =
         (int)(it->second * getNumberOfSimHits(*simTrack) + 0.5);

        if(nShared >= 3) nRec++;
      }
    }
    catch (cms::Exception& event)
    { }

    result.push_back(nRec);

    // fill
    trackSim->Fill(&result[0]);
  }
}

/*****************************************************************************/
void LowPtTrackAnalyzer::checkRecTracks
  (edm::Handle<reco::TrackCollection>& recCollection,
   reco::RecoToSimCollection& p)
{
  for(reco::TrackCollection::size_type i=0;
          i < recCollection.product()->size(); ++i)
  {
    reco::TrackRef recTrack(recCollection, i);

    vector<float> result;

    // rec
    result.push_back(recTrack->charge());
    result.push_back(recTrack->eta());
    result.push_back(recTrack->pt());
    result.push_back(recTrack->d0());

    // sim 
    TrackingParticleRef matchedSimTrack;
    bool match = false;

    int nSim = 0;
    try
    {
      vector<pair<TrackingParticleRef, double> > simTracks = p[recTrack];

      for(vector<pair<TrackingParticleRef, double> >::const_iterator
            it = simTracks.begin(); it != simTracks.end(); ++it)
      {
        TrackingParticleRef simTrack = it->first;
        float fraction = it->second;

        // If all hits are shared 
        if(fraction == 1.)
        { matchedSimTrack = simTrack; match = true; nSim = 1 ; break; }

        // If some, but not all hits are shared, increase nsim
        if(fraction > 0. && fraction < 1.) nSim++;
      }
    }
    catch (cms::Exception& event)
    { }
    result.push_back(nSim);

    if(match)
    {
      int parentId;
  
      if(matchedSimTrack->parentVertex()->nSourceTracks() == 1)
      {
        TrackingVertex::tp_iterator iv =
          matchedSimTrack->parentVertex()->sourceTracks_begin();
        parentId = (*iv)->pdgId();
      }
      else parentId = 0;

      result.push_back(matchedSimTrack->pdgId());
      result.push_back(parentId);
      result.push_back(matchedSimTrack->eta());
      result.push_back(matchedSimTrack->pt());
      result.push_back(matchedSimTrack->vertex().Rho());
    }

    // fill
    trackRec->Fill(&result[0]);
  }
}

/*****************************************************************************/
void LowPtTrackAnalyzer::analyze
  (const edm::Event& ev, const edm::EventSetup& es)
{
  // Get simulated
  edm::Handle<TrackingParticleCollection> simCollection;
  ev.getByType(simCollection);

  // Get reconstructed
  edm::Handle<reco::TrackCollection> recCollection;
  ev.getByLabel("ctfTripletTracks",  recCollection);

  // Associators
  reco::SimToRecoCollection simToReco =
    theAssociatorByHits->associateSimToReco(recCollection, simCollection,&ev);
  reco::RecoToSimCollection recoToSim =
    theAssociatorByHits->associateRecoToSim(recCollection, simCollection,&ev);

  // Analyze
  checkSimTracks(simCollection, simToReco);
  checkRecTracks(recCollection, recoToSim);

  // Plot event
  if(plotEvent)
  {
    EventPlotter theEventPlotter(es);
    theEventPlotter.printEvent(ev);

    if(zipFiles)
    {
      system("zip -q -m event.zip *.m");
      cerr << "[LowPtTrackAnalyzer] event plotted (event.zip)" << endl;
    }
    else
      cerr << "[LowPtTrackAnalyzer] event plotted (*.m)" << endl;

    cerr << "[LowPtTrackAnalyzer] done " << ev.id()          << endl;
    cerr << "----------------------------------------------" << endl;

    while(getchar() == 0);
  }
  else
  {
    cerr << "[LowPtTrackAnalyzer] done " << ev.id()          << endl;
    cerr << "----------------------------------------------" << endl;
  }
}

DEFINE_FWK_MODULE(LowPtTrackAnalyzer);
