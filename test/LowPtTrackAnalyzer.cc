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

//#include "RecoPixelVertexing/PixelLowPtUtilities/interface/EventPlotter.h"

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

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

/*
class TransientTrackFromFTSFactory {
 public:

    reco::TransientTrack build (const FreeTrajectoryState & fts) const;
    reco::TransientTrack build (const FreeTrajectoryState & fts,
        const edm::ESHandle<GlobalTrackingGeometry>& trackingGeometry);
};
*/

#include "RecoVertex/KalmanVertexFit/interface/SingleTrackVertexConstraint.h"

#include "TROOT.h"
#include "TFile.h"
#include "TNtuple.h"

#include <fstream>
using namespace std;
using namespace reco;

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
   void checkSimTracks (edm::Handle<TrackingParticleCollection>& simCollection,
                       reco::SimToRecoCollection& q);

   pair<float,float> refitWithVertex(const reco::Track & recTrack,
                                     const reco::VertexCollection* vertices);
   int getParticleId(edm::RefToBase<reco::Track>& recTrack, int & ptype);
   void checkRecTracks(edm::Handle<edm::View<reco::Track> >& recCollection,
                       const reco::VertexCollection* vertices,
                       reco::RecoToSimCollection& p);

   const TrackerGeometry * theTracker;
   const TrackAssociatorByHits * theAssociatorByHits;
   const TransientTrackBuilder * theTTBuilder;
   TrackerHitAssociator * theHitAssociator;

   TFile * resultFile; 
   TNtuple * trackSim;
   TNtuple * trackRec;
   TNtuple * eventInfo;

   vector<string> trackCollectionLabels;
   string resultFileLabel;
   bool plotEvent, zipFiles;
   int proc, ntrk,nvtx;
};

/*****************************************************************************/
LowPtTrackAnalyzer::LowPtTrackAnalyzer(const edm::ParameterSet& pset)
{
  trackCollectionLabels = pset.getParameter<vector<string> >("trackCollection");
  resultFileLabel       = pset.getParameter<string>("resultFile");

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
 
  // Get associator
  edm::ESHandle<TrackAssociatorBase> theHitsAssociator;
  es.get<TrackAssociatorRecord>().get("TrackAssociatorByHits",
                                    theHitsAssociator);
  theAssociatorByHits =
   (const TrackAssociatorByHits*)theHitsAssociator.product();

  // Get transient track builder
  edm::ESHandle<TransientTrackBuilder> builder;
  es.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
  theTTBuilder = builder.product();
 
  // Root
  resultFile = new TFile(resultFileLabel.c_str(),"recreate");
  resultFile->cd();
   
  trackSim = new TNtuple("trackSim","trackSim",
    "proc:ntrk:nvtx:ids:prim:parids:etas:pts:rhos:nhits:nrec:qr:etar:ptr:d0r");

  trackRec = new TNtuple("trackRec","trackRec",
    "proc:ntrk:nvtx:qr:etar:ptr:chir:ptv:chiv:nhitr:d0r:nsim:ids:ptype:prim:parids:etas:pts:rhos");

  eventInfo = new TNtuple("eventInfo","eventInfo",
    "proc:ntrkr:nvtxr");
//    "proc:ntrks:ntrkr:nvtxr");
}

/*****************************************************************************/
void LowPtTrackAnalyzer::endJob()
{
  resultFile->cd();
  trackSim->Write();
  trackRec->Write();
  eventInfo->Write();
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
int LowPtTrackAnalyzer::getNumberOfSimHits(const TrackingParticle& constSimTrack)
{
  int oldlay = 0; int newlay = 0;
  int olddet = 0; int newdet = 0;

  int nhit = 0;

  // to use trackerPSimHit methods TP cannot be const
  TrackingParticle simTrack = *(const_cast<TrackingParticle*>(&constSimTrack));
  
  for(std::vector<PSimHit>::const_iterator
       simHit = simTrack.trackerPSimHit_begin();
       simHit!= simTrack.trackerPSimHit_end(); simHit++)
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
int LowPtTrackAnalyzer::getNumberOfPixelHits
  (const TrackingParticle& simTrack)
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

    if(simTrack->charge() != 0) 
    {
    vector<float> result;

    //
    result.push_back(proc);
    result.push_back(ntrk);
    result.push_back(nvtx);

    // sim
    result.push_back(simTrack->pdgId());               // ids
    result.push_back(simTrack->parentVertex()->position().T()); // ?
                                                       // prim
    result.push_back(simTrack->pdgId());               // parids
    result.push_back(simTrack->eta());                 // etas
    result.push_back(simTrack->pt());                  // pts
    result.push_back(simTrack->vertex().Rho());        // rhos

    result.push_back(getNumberOfPixelHits(*simTrack)); // nhits >= 3 accepted

    // rec
    edm::RefToBase<reco::Track> matchedRecTrack;
    int nRec = 0;

    try
    {
      vector<pair<edm::RefToBase<reco::Track>, double> > recTracks = q[simTrack];

      for(vector<pair<edm::RefToBase<reco::Track>,double> >::const_iterator
            it = recTracks.begin(); it != recTracks.end(); ++it)
      {
        edm::RefToBase<reco::Track> recTrack = it->first;
        int nShared =
         (int)(it->second * getNumberOfSimHits(*simTrack) + 0.5);

        if( (recTrack->found() >= 3 &&
               nShared >  recTrack->found() / 2 && nShared >= 3) ||
            (recTrack->found() == 2 &&
               nShared == recTrack->found()) )
        { matchedRecTrack = recTrack; nRec++; }
      }
    }
    catch (cms::Exception& event)
    { }

    result.push_back(nRec); // nrec

    if(nRec > 0)
    {
      result.push_back(matchedRecTrack->charge()); // qr
      result.push_back(matchedRecTrack->eta());    // etar
      result.push_back(matchedRecTrack->pt());     // ptr
      result.push_back(matchedRecTrack->d0());     // d0r
    }
    else
    {
      result.push_back(-9999);
      result.push_back(-9999);
      result.push_back(-9999);
      result.push_back(-9999);
    }


    // fill
    trackSim->Fill(&result[0]);
    }
  }
}

/*****************************************************************************/
pair<float,float> LowPtTrackAnalyzer::refitWithVertex
  (const reco::Track & recTrack,
   const reco::VertexCollection* vertices)
{
  TransientTrack theTransientTrack = theTTBuilder->build(recTrack);

  // If there are vertices found
  if(vertices->size() > 0)
  {
    float dzmin = -1.;
    const reco::Vertex * closestVertex = 0;

    // Look for the closest vertex in z
    for(reco::VertexCollection::const_iterator
        vertex = vertices->begin(); vertex!= vertices->end(); vertex++)
    {
      float dz = fabs(recTrack.vertex().z() - vertex->position().z());
      if(vertex == vertices->begin() || dz < dzmin)
      { dzmin = dz ; closestVertex = &(*vertex); }
    }
  
    // Get vertex position and error matrix
    GlobalPoint vertexPosition(closestVertex->position().x(),
                               closestVertex->position().y(),
                               closestVertex->position().z());

    float beamSize = 15e-4; // 15 um
    GlobalError vertexError(beamSize*beamSize, 0,
                            beamSize*beamSize, 0,
                            0,closestVertex->covariance(2,2));

    // Refit track with vertex constraint
    SingleTrackVertexConstraint stvc;
    pair<TransientTrack, float> result =
      stvc.constrain(theTransientTrack, vertexPosition, vertexError);

    return pair<float,float>(result.first.impactPointTSCP().pt(),
                             result.second);
  }
  else
    return pair<float,float>(recTrack.pt(), -9999);
}

/*****************************************************************************/
int LowPtTrackAnalyzer::getParticleId(edm::RefToBase<reco::Track>& recTrack, int & ptype)
{
  int pid = 0;
  ptype = 0;
  double tmin = 0.;

  for(trackingRecHit_iterator recHit = recTrack->recHitsBegin();
                              recHit!= recTrack->recHitsEnd(); recHit++)
  {
    vector<PSimHit> simHits = theHitAssociator->associateHit(**recHit);

    for(vector<PSimHit>::const_iterator simHit = simHits.begin(); 
                                        simHit!= simHits.end(); simHit++)
      if(simHit == simHits.begin() || simHit->tof() < tmin )
      {
        pid   = simHit->particleType();
        ptype = simHit->processType();
        tmin  = simHit->tof();
      }
  }  

  return pid;
}

/*****************************************************************************/
void LowPtTrackAnalyzer::checkRecTracks
  (edm::Handle<edm::View<reco::Track> >& recCollection,
   const reco::VertexCollection* vertices,
   reco::RecoToSimCollection& p)
{
  for(edm::View<reco::Track> ::size_type i=0;
          i < recCollection.product()->size(); ++i)
  {
    edm::RefToBase<reco::Track> recTrack(recCollection, i);

    vector<float> result;

    //
    result.push_back(proc); // proc
    result.push_back(ntrk); // ntrk
    result.push_back(nvtx); // nvtx

    // rec
    result.push_back(recTrack->charge());                         // qr
    result.push_back(recTrack->eta());                            // etar
    result.push_back(recTrack->pt());                             // ptr
    result.push_back(recTrack->chi2());                           // chir
    result.push_back(refitWithVertex(*recTrack,vertices).first);  // ptv
    result.push_back(refitWithVertex(*recTrack,vertices).second); // chiv
    result.push_back(recTrack->numberOfValidHits());              // nhitr
    result.push_back(recTrack->d0());                             // d0r

    // sim 
    TrackingParticleRef matchedSimTrack;
    int nSim = 0;

    try
    {
      vector<pair<TrackingParticleRef, double> > simTracks = p[recTrack];

      for(vector<pair<TrackingParticleRef, double> >::const_iterator
            it = simTracks.begin(); it != simTracks.end(); ++it)
      {
        TrackingParticleRef simTrack = it->first;
        float fraction = it->second;

        // If more than half is shared
        if(fraction > 0.5)
        { matchedSimTrack = simTrack; nSim++; }
      }
    }
    catch (cms::Exception& event)
    { }

    result.push_back(nSim); // nsim

    if(nSim > 0)
    {
      int parentId;
      float T;

      int ptype;
      int ids = getParticleId(recTrack, ptype);

      if(matchedSimTrack->parentVertex()->nSourceTracks() == 0)
      {
        // track is primary, has no parent
        // recTrack can be a true primary, or an untracked daughter
        if(ptype == 2) parentId = 0;                        // primary
                  else parentId = matchedSimTrack->pdgId(); // hadronic, decay
      }
      else
      {
        // track is not primary, has a parent
        TrackingVertex::tp_iterator iv =
          matchedSimTrack->parentVertex()->sourceTracks_begin();
        parentId = (*iv)->pdgId();
      }

      T = matchedSimTrack->parentVertex()->position().T(); // ?

      result.push_back(ids);                              // ids
      result.push_back(ptype);                            // ptype
      result.push_back(T);                                // prim
      result.push_back(parentId);                         // parids
      result.push_back(matchedSimTrack->eta());           // etas
      result.push_back(matchedSimTrack->pt());            // pts
      result.push_back(matchedSimTrack->vertex().Rho());  // rhos
    }
    else 
    {
      result.push_back(-9999); // ids
      result.push_back(-9999); // ptype
      result.push_back(-9999); // prim
      result.push_back(-9999); // parids
      result.push_back(-9999); // etas
      result.push_back(-9999); // pts
      result.push_back(-9999); // rhos
    }

    // fill
    trackRec->Fill(&result[0]);
  }
}

/*****************************************************************************/
void LowPtTrackAnalyzer::analyze
  (const edm::Event& ev, const edm::EventSetup& es)
{
  // Get associator
  theHitAssociator = new TrackerHitAssociator::TrackerHitAssociator(ev);

  // Get generated
  edm::Handle<edm::HepMCProduct> hepEv;
  ev.getByType(hepEv);
  proc = hepEv->GetEvent()->signal_process_id();
  cerr << "[LowPtTrackAnalyzer] process = " << proc << endl;

  // Get simulated
  edm::Handle<TrackingParticleCollection> simCollection;
//  ev.getByLabel("trackingtruthprod",simCollection);
  ev.getByType(simCollection);

  // Get reconstructed
  edm::Handle<edm::View<reco::Track> >  recCollection;
  ev.getByLabel(trackCollectionLabels[0], recCollection); // !!

  cerr << "[LowPtTrackAnalyzer] recTracks = "
       << recCollection.product()->size() << endl;

  // Get vertices
  edm::Handle<reco::VertexCollection> vertexCollection;
  ev.getByLabel("pixelVertices",vertexCollection);
  const reco::VertexCollection * vertices = vertexCollection.product();

  ntrk = recCollection.product()->size();
  nvtx = vertexCollection.product()->size();

  {
  vector<float> result;
  result.push_back(proc); // proc
  result.push_back(ntrk); // ntrkr
  result.push_back(nvtx); // nvtxr

  eventInfo->Fill(&result[0]);
  }

  // Associators
  reco::SimToRecoCollection simToReco =
    theAssociatorByHits->associateSimToReco(recCollection, simCollection,&ev);
  reco::RecoToSimCollection recoToSim =
    theAssociatorByHits->associateRecoToSim(recCollection, simCollection,&ev);

  // Analyze
  checkSimTracks(simCollection,           simToReco);
  checkRecTracks(recCollection, vertices, recoToSim);

  cerr << "[LowPtTrackAnalyzer] done, " << ev.id()          << endl;
  cerr << "----------------------------------------------" << endl;

  // Plot event
/*
  if(plotEvent)
  {
    EventPlotter theEventPlotter(es, trackCollectionLabels);
    theEventPlotter.printEvent(ev);

    if(zipFiles)
    {
      system("zip -q -m event.zip *.m");
      cerr << "[LowPtTrackAnalyzer] event plotted (event.zip)" << endl;
    }
    else
      cerr << "[LowPtTrackAnalyzer] event plotted (*.m)" << endl;

    cerr << "[LowPtTrackAnalyzer] done, " << ev.id()          << endl;
    cerr << "----------------------------------------------" << endl;

    while(getchar() == 0);
  }
  else
  {
    cerr << "[LowPtTrackAnalyzer] done, " << ev.id()          << endl;
    cerr << "----------------------------------------------" << endl;
  }
*/

  delete theHitAssociator;
}

DEFINE_FWK_MODULE(LowPtTrackAnalyzer);
