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
   void checkSimTracks(edm::Handle<TrackingParticleCollection>& simCollection);
   void checkRecTracks(edm::Handle<reco::TrackCollection>& recCollection);

   TFile * resultFile; 
   TNtuple * trackSim;
   TNtuple * trackRec;
};

/*****************************************************************************/
LowPtTrackAnalyzer::LowPtTrackAnalyzer(const edm::ParameterSet& pset)
{
}

/*****************************************************************************/
LowPtTrackAnalyzer::~LowPtTrackAnalyzer()
{
}

/*****************************************************************************/
void LowPtTrackAnalyzer::beginJob(const edm::EventSetup& es)
{
  resultFile = new TFile("result.root","recreate");
  resultFile->cd();

  trackSim = new TNtuple("trackSim","trackSim", "id:eta:pt:d0");
  trackRec = new TNtuple("trackRec","trackRec",  "q:eta:pt:d0");
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
void LowPtTrackAnalyzer::checkSimTracks
  (edm::Handle<TrackingParticleCollection>& simCollection)
{
  for(TrackingParticleCollection::size_type i=0;
                 i<simCollection.product()->size(); ++i)
  {
    TrackingParticleRef simTrack(simCollection, i);

    vector<float> result;

    // sim
    result.push_back(simTrack->pdgId());
    result.push_back(simTrack->eta());
    result.push_back(simTrack->pt());
    result.push_back(simTrack->vertex().Rho());

    // fill
    trackSim->Fill(&result[0]);
  }
}

/*****************************************************************************/
void LowPtTrackAnalyzer::checkRecTracks
  (edm::Handle<reco::TrackCollection>& recCollection)
{
  for(reco::TrackCollection::size_type i=0;
            i<recCollection.product()->size(); ++i)
  {
    reco::TrackRef recTrack(recCollection, i);

    vector<float> result;

    // rec
    result.push_back(recTrack->charge());
    result.push_back(recTrack->eta());
    result.push_back(recTrack->pt());
    result.push_back(recTrack->d0());

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
  ev.getByType(recCollection);

  // Check simulated
  checkSimTracks(simCollection);

  // Check reconstructed
  checkRecTracks(recCollection);
}

DEFINE_FWK_MODULE(LowPtTrackAnalyzer);
