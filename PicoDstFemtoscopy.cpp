// This is needed for calling standalone classes
#define _VANILLA_ROOT_
// C++ headers
#include <iostream>
#include <map>
#include <vector>
#include <list>
// ROOT headers
#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TMath.h>
#include <TVector3.h>
#include <Math/Vector4D.h>
#include <TLorentzVector.h>
#include <TGraph.h>
#include <TStopwatch.h>
// PicoDst headers
#include "StPicoDstReader.h"
#include "StPicoDst.h"
#include "StPicoEvent.h"
#include "StPicoTrack.h"
#include "StPicoBTofPidTraits.h"
#include "StRefMultCorr.h"

#include "constants.C"
Int_t harmonic = 2;                   // set harmonic for eta-sub event plane, Q-Cumulants, and scalar product method
Int_t debug = 1;
// #include "StUtilities.C"
using std::cout;
using std::cerr;
using std::endl;
using std::vector;

struct myPart
{
  myPart(TLorentzVector fourMom, UInt_t topoMap0, UInt_t topoMap1, Int_t nHits, StPicoPhysicalHelix hel)
   : m4Mom(fourMom),mTopoMap0(topoMap0),mTopoMap1(topoMap1),mNhits(nHits),mHelix(hel) {}
 public:
  TLorentzVector      m4Mom;
  UInt_t              mTopoMap0;
  UInt_t              mTopoMap1;
  Int_t               mNhits;
  StPicoPhysicalHelix mHelix;
//  private:
};

Double_t Getkt(const TVector3 &p1, const TVector3 &p2) {
  // return TMath::Sqrt(TMath::Power(p1.Px()+p2.Px(),2.) + TMath::Power(p1.Py()+p2.Py(),2.))/2.;
  return (p1+p2).Perp()/2.;
}

Double_t GetSpittingLevel(const myPart &tr1, const myPart &tr2) {
  // splittingLevel = OL/Nhits
  // Nhits - total number of hits of both tracks
  // The occupancy level (OL) of the pad rows
  // is the sum over 45 rows of A_i, where 
  // A = -1 if both tracks have hits on i-th row
  // A = 0  if non of the tracks have hits on i-th row
  // A = 1  if only 1 of 2 tracks has hit on i-th row

  Int_t numerator = 0;
  for (Int_t i=8; i<32; i++) {
    if      ( tr1.mTopoMap0>>i&1U &&  tr2.mTopoMap0>>i&1U) numerator--;
    else if (!tr1.mTopoMap0>>i&1U && !tr2.mTopoMap0>>i&1U) continue;
    else numerator++;
  }
  for (Int_t i=0; i<21; i++) {
    if      ( tr1.mTopoMap1>>i&1U &&  tr2.mTopoMap1>>i&1U) numerator--;
    else if (!tr1.mTopoMap1>>i&1U && !tr2.mTopoMap1>>i&1U) continue;
    else numerator++;
  }
  return (Double_t) numerator/(tr1.mNhits + tr2.mNhits);
}

Double_t GetAverageDistance(const myPart &tr1, const myPart &tr2) {
  Double_t tpc_radius;
  StPicoPhysicalHelix helix1 = tr1.mHelix;
  StPicoPhysicalHelix helix2 = tr2.mHelix;
  Double_t aveDist = 0.;
  for (UInt_t i=0;i<11;i++) {
    tpc_radius = 50. + i*15.;
    pair<Double_t,Double_t> pathLengthPi1 = helix1.pathLength(tpc_radius);
    pair<Double_t,Double_t> pathLengthPi2 = helix2.pathLength(tpc_radius);
    Double_t pathLPi1 = pathLengthPi1.first > 0 ? pathLengthPi1.first : pathLengthPi1.second;
    Double_t pathLPi2 = pathLengthPi2.first > 0 ? pathLengthPi2.first : pathLengthPi2.second;
    aveDist += (helix1.at(pathLPi1) - helix2.at(pathLPi2)).Mag();
  } // for (UInt_t i=0;i<11;i++)
  return aveDist/11.;
}

Bool_t isGoodPID(StPicoTrack *const &track, StPicoBTofPidTraits *const &trait)
{
  // Check if track has TOF signal
  if (!track->isTofTrack()) return false;
  if (!trait) { std::cout << "Warning: There's no BTofPidTrait # " 
                          << track->bTofPidTraitsIndex() << " for track # " 
                          << track->id() << std::endl; return false; }
  // if (track->massSqr() < cutMass2Min) return false;
  return true;
}

Int_t  GetPID(StPicoTrack *const &track, StPicoBTofPidTraits *const &trait);

void PicoDstFemtoscopy(const Char_t *inFile = "", const Char_t *outFile = "", Float_t colEnergy = 27.)
{
  cout << "Hi! Lets do some physics, Master!" << endl;
  TStopwatch timer;
  timer.Start();
  const Bool_t bFXT = kTRUE;
  const Float_t pion_mass = 0.13957061;
  const UInt_t nkt = 4;
  const Double_t ktBin[nkt+1] = {0.2, 0.4, 0.6, 0.8, 1.0};
  const UInt_t nSL = 3;
  const Double_t SL[nSL] = {0.8, 0.6, 0.4};
  const UInt_t nAD = 3;
  const Double_t AD[nAD] = {4,6,8}; // in centimeter

  // ============================================================= //
  // ================= Mix event configuration =================== //
  // ============================================================= //
  const int nBuffer = 5;  // length of the mixing buffer
  const int nMyCent = 4;  // 0-10, 10-30, 30-50, 50-80
  const int nVtxZ = 1;    // for FXT
  std::list< vector< myPart > > mixBuf[nMyCent][nVtxZ];
  for (int ic(0);ic<nMyCent;ic++) {
    for (int iv(0);iv<nVtxZ;iv++) {
      mixBuf[ic][iv].clear(); // just to be sure the initialized lists are empty
    } // for (int iv(0);iv<nVtxZ;iv++)
  } // for (int ic(0);ic<nMyCent;ic++)


  vector<myPart> pion;
  // Set up output file
  TFile *fo = new TFile(outFile, "recreate");
  TH1D *hQinvNum[nMyCent];
  TH1D *hQinvDen[nMyCent];
  TH1D *hQinvNum_kt[nMyCent][nkt];
  TH1D *hQinvDen_kt[nMyCent][nkt];

  for (UInt_t icent=0;icent<nMyCent;icent++) {
    hQinvNum[icent] = new TH1D(Form("hQinvNum_%i",icent), Form(";q;Entries"), 100, 0., 1.);
    hQinvDen[icent] = new TH1D(Form("hQinvDen_%i",icent), Form(";q;Entries"), 100, 0., 1.);
    for (UInt_t ikt=0;ikt<nkt;ikt++) {
      hQinvNum_kt[icent][ikt] = new TH1D(Form("hQinvNum_kt_%i_%i",icent,ikt), Form(";q;Entries"), 100, 0., 1.);
      hQinvDen_kt[icent][ikt] = new TH1D(Form("hQinvDen_kt_%i_%i",icent,ikt), Form(";q;Entries"), 100, 0., 1.);
    }
  } // for (Int_t icent=0;icent<nMyCent;icent++)

  TH2D *hSLNum[nMyCent];
  TH2D *hSLDen[nMyCent];
  TH2D *hSLNum_kt[nMyCent][nkt];
  TH2D *hSLDen_kt[nMyCent][nkt];
  for (UInt_t icent=0;icent<nMyCent;icent++) {
    hSLNum[icent] = new TH2D(Form("hSLNum_%i",icent), Form(";q;splittingLevel;Entries"), 100, 0., 1.,300,-0.5,1.0);
    hSLDen[icent] = new TH2D(Form("hSLDen_%i",icent), Form(";q;splittingLevel;Entries"), 100, 0., 1.,300,-0.5,1.0);
    for (UInt_t ikt=0;ikt<nkt;ikt++) {
      hSLNum_kt[icent][ikt] = new TH2D(Form("hSLNum_kt_%i_%i",icent,ikt), Form(";q;splittingLevel;Entries"), 100, 0., 1.,300,-0.5,1.0);
      hSLDen_kt[icent][ikt] = new TH2D(Form("hSLDen_kt_%i_%i",icent,ikt), Form(";q;splittingLevel;Entries"), 100, 0., 1.,300,-0.5,1.0);
    } // for (UInt_t ikt=0;ikt<nkt;ikt++)
  } // for (UInt_t icent=0;icent<nMyCent;icent++)

  TH2D *hADNum[nMyCent];
  TH2D *hADDen[nMyCent];
  TH2D *hADNum_kt[nMyCent][nkt];
  TH2D *hADDen_kt[nMyCent][nkt];
  TProfile *pADNum[nMyCent];
  TProfile *pADDen[nMyCent];
  for (UInt_t icent=0;icent<nMyCent;icent++) {
    hADNum[icent] = new TH2D(    Form("hADNum_%i",icent), Form(";q;Average distance, cm;Entries"), 100, 0., 1.,1000,0.,10.0);
    hADDen[icent] = new TH2D(    Form("hADDen_%i",icent), Form(";q;Average distance, cm;Entries"), 100, 0., 1.,1000,0.,10.0);
    pADNum[icent] = new TProfile(Form("pADNum_%i",icent), Form(";q;Average distance, cm"),         100, 0., 1.             );
    pADDen[icent] = new TProfile(Form("pADDen_%i",icent), Form(";q;Average distance, cm"),         100, 0., 1.             );
    for (UInt_t ikt=0;ikt<nkt;ikt++) {
      hADNum_kt[icent][ikt] = new TH2D(Form("hADNum_kt_%i_%i",icent,ikt), Form(";q;Average distance, cm;Entries"), 100, 0., 1.,1000,0.,10.0);
      hADDen_kt[icent][ikt] = new TH2D(Form("hADDen_kt_%i_%i",icent,ikt), Form(";q;Average distance, cm;Entries"), 100, 0., 1.,1000,0.,10.0);
    } // for (UInt_t ikt=0;ikt<nkt;ikt++)

  } // for (Int_t icent=0;icent<nMyCent;icent++)

  TH1D *hQinvNum_SL[nMyCent][nSL];
  TH1D *hQinvDen_SL[nMyCent][nSL];
  for (UInt_t icent=0;icent<nMyCent;icent++) {
    for (UInt_t iSL=0;iSL<nSL;iSL++) {
      hQinvNum_SL[icent][iSL] = new TH1D(Form("hQinvNum_SL_%i_%i",icent,iSL), Form(";q;Entries"), 100, 0., 1.);
      hQinvDen_SL[icent][iSL] = new TH1D(Form("hQinvDen_SL_%i_%i",icent,iSL), Form(";q;Entries"), 100, 0., 1.);
    } // for (UInt_t iSL=0;iSL<nSL;iSL++)
  } // for (UInt_t icent=0;icent<nMyCent;icent++)


  TH1D *hQinvNum_AD[nMyCent][nAD];
  TH1D *hQinvDen_AD[nMyCent][nAD];
  for (UInt_t icent=0;icent<nMyCent;icent++) {
    for (UInt_t iAD=0;iAD<nAD;iAD++) {
      hQinvNum_AD[icent][iAD] = new TH1D(Form("hQinvNum_AD_%i_%i",icent,iAD), Form(";q;Entries"), 100, 0., 1.);
      hQinvDen_AD[icent][iAD] = new TH1D(Form("hQinvDen_AD_%i_%i",icent,iAD), Form(";q;Entries"), 100, 0., 1.);
    } // for (UInt_t iAD=0;iAD<nAD;iAD++)
  } // for (UInt_t icent=0;icent<nMyCent;icent++)

  TH1D *hQinvNum_SL_AD[nMyCent][nAD];
  TH1D *hQinvDen_SL_AD[nMyCent][nAD];
  for (UInt_t icent=0;icent<nMyCent;icent++) {
    for (UInt_t iAD=0;iAD<nAD;iAD++) {
      hQinvNum_SL_AD[icent][iAD] = new TH1D(Form("hQinvNum_SL_AD_%i_%i",icent,iAD), Form(";q;Entries"), 100, 0., 1.);
      hQinvDen_SL_AD[icent][iAD] = new TH1D(Form("hQinvDen_SL_AD_%i_%i",icent,iAD), Form(";q;Entries"), 100, 0., 1.);
    } // for (UInt_t iAD=0;iAD<nAD;iAD++)
  } // for (UInt_t icent=0;icent<nMyCent;icent++)

  TH1D *hDeltaRNum[nMyCent];
  TH1D *hDeltaRDen[nMyCent];
  TH1D *hDeltaRNum_kt[nMyCent][nkt];
  TH1D *hDeltaRDen_kt[nMyCent][nkt];
  TH1D *hDeltaRNum_SL[nMyCent][nSL];
  TH1D *hDeltaRDen_SL[nMyCent][nSL];
  for (UInt_t icent=0;icent<nMyCent;icent++) {
    hDeltaRNum[icent] = new TH1D(Form("hDeltaRNum_%i",icent), Form(";Average distance, cm;Entries"),         200, 0., 20.             );
    hDeltaRDen[icent] = new TH1D(Form("hDeltaRDen_%i",icent), Form(";Average distance, cm;Entries"),         200, 0., 20.             );
    for (UInt_t ikt=0;ikt<nkt;ikt++) {
      hDeltaRNum_kt[icent][ikt] = new TH1D(Form("hDeltaRNum_kt_%i_%i",icent,ikt), Form(";Average distance, cm;Entries"), 200, 0., 20.);
      hDeltaRDen_kt[icent][ikt] = new TH1D(Form("hDeltaRDen_kt_%i_%i",icent,ikt), Form(";Average distance, cm;Entries"), 200, 0., 20.);
    } // for (UInt_t ikt=0;ikt<nkt;ikt++)
    for (UInt_t iSL=0;iSL<nSL;iSL++) {
      hDeltaRNum_SL[icent][iSL] = new TH1D(Form("hDeltaRNum_SL_%i_%i",icent,iSL), Form(";Average distance, cm;Entries"),200, 0., 20.);
      hDeltaRDen_SL[icent][iSL] = new TH1D(Form("hDeltaRDen_SL_%i_%i",icent,iSL), Form(";Average distance, cm;Entries"),200, 0., 20.);
    } // for (UInt_t iSL=0;iSL<nSL;iSL++)

  } // for (Int_t icent=0;icent<nMyCent;icent++)

  const char *eventType[8] = {"RAW","Min. bias trigger","Bad-run rejection","Pile-up rejection","V_{z} cut","V_{r} cut","|V_{z}-V_{z}^{VPD}| => OFF","Centrality cut"};
  TH1D* hEvt = new TH1D("hEvt",";;Entries",8,0,8);
  hEvt->SetCanExtend(TH1::kAllAxes);
  TH1D *hVz;
  if (bFXT) hVz = new TH1D("hVz",";V_{z}, cm",100,190.,210.);
  else      hVz = new TH1D("hVz",";V_{z}, cm",200,-100.,100.);
  TH2D *hVr = new TH2D("hVr",";V_{x}, cm;V_{y}, cm",100,-5.,5.,100,-5.,5.);
  TH2D *hPidEdxpq = new TH2D("hPidEdxpq",";Rigidity p/q, GeV/c;dE/dx in TPC, (keV/cm)", 300, -5., 5., 300, 0., 100.);    
  TH2D *hPiM2pq = new TH2D("hPiM2pq",";Rigidity p/q, GeV/c;m^{2}, (GeV/c^{2})^{2}", 300, -5., 5., 300, -0.8, 1.8);    
  TH1D *hPiPtot_TPC_only = new TH1D("hPiPtot_TPC_only",";p, GeV/c;Entries",400, 0., 4.);
  TH1D *hPiPtot_TPC_TOF = new TH1D("hPiPtot_TPC_TOF",";p, GeV/c;Entries",400, 0., 4.);
  TH1D *hPiPtot = new TH1D("hPiPtot",";p, GeV/c;Entries",400, 0., 4.);
  TH2D *hPiEtaPt = new TH2D(Form("hPiEtaPt"), Form(";p_{T}, GeV/c;#eta;Entries"), 400, 0., 4., 400, -1.2, 1.2);

  int nTracks, pid, icent, ikt, starCent, iVtxZ, iFxtMult;
  Double_t kt, charge, qinv, splittingLevel, aveDist, ptot, pt, eta, phi;
  TLorentzVector lrtzVecPi;
  Float_t bField;

  StPicoDstReader* picoReader = new StPicoDstReader(inFile);
  picoReader->Init();
  cout << "Explicit read status for some branches" << endl;
  picoReader->SetStatus("*",0);
  picoReader->SetStatus("Event*", 1);
  picoReader->SetStatus("Track*", 1);
  picoReader->SetStatus("BTofHit*", 1);
  picoReader->SetStatus("BTofPidTraits*", 1);
  cout << "Status has been set" << endl;
  cout << "Now I know what to read, Master!" << endl;
  if( !picoReader->chain() ) cout << "No chain has been found." << endl;
  Long64_t eventsInTree = picoReader->tree()->GetEntries();
  cout << "eventsInTree: "  << eventsInTree << endl;
  Long64_t events2read = picoReader->chain()->GetEntries();
  cout << "Number of events to read: " << events2read << endl;
  cout << "Starting filling histograms for femtoscopy analysis..." << endl;
  // StRefMultCorr *refmultCorrUtil = new StRefMultCorr("refMult");
  // Loop over events
  for(Long64_t iEvent=0; iEvent<events2read; iEvent++) {
    if (iEvent % 1000 == 0) cout << "Woking on event #[" << (iEvent+1) << "/" << events2read << "]" << endl;
    Bool_t readEvent = picoReader->readPicoEvent(iEvent);
    if( !readEvent ) { cout << "Something went wrong, Master! Nothing to analyze..." << endl; break; }
    StPicoDst *dst = picoReader->picoDst();
    StPicoEvent *event = dst->event();
    if( !event ) { cout << "Something went wrong, Master! Event is hiding from me..." << endl; break; }
    hEvt->Fill(eventType[0], 1);
    for (auto mbtrig : MinBiasTrigger.at(colEnergy)) if (!event->isTrigger(mbtrig)) continue;
    hEvt->Fill(eventType[1], 1);
    // refmultCorrUtil->init(event -> runId());
    // if (refmultCorrUtil->isBadRun(event -> runId())) continue;
    std::vector<Int_t> Bad = BadRuns.at(colEnergy);
    if( std::find(Bad.begin(), Bad.end(), event -> runId()) != Bad.end() )  continue;
    hEvt->Fill(eventType[2], 1);
    TVector3 pVtx = event->primaryVertex();
    // if (refmultCorrUtil->isPileUpEvent(event->refMult(),event->nBTOFMatch(),pVtx.Z()))            continue; // pile-up remover
    hEvt->Fill(eventType[3], 1);
    if ( bFXT && (event->primaryVertex().Z() > 205 || event->primaryVertex().Z() < 195) ) continue;
    if (!bFXT && TMath::Abs(event->primaryVertex().Z()) > cutVtxZEnergy.at(colEnergy)) continue;
    hEvt->Fill(eventType[4], 1);
    if (!bFXT && event->primaryVertex().Perp() > cutVtxR)                                                  continue;
    if ( bFXT && TMath::Sqrt(TMath::Power(event->primaryVertex().Y()-VyShift,2.) + TMath::Power(event->primaryVertex().X()-VxShift,2.)) > cutVtxR) continue;
    hEvt->Fill(eventType[5], 1);
    // if (TMath::Abs(event->primaryVertex().Z() - event->vzVpd()) > cutVpdVz)                       continue;    
    hEvt->Fill(eventType[6], 1);
    // if (!bFXT) refmultCorrUtil->initEvent(event->refMult(),pVtx.Z(),event->ZDCx());
    // else refmultCorrUtil->initEvent(event->fxtMult(),pVtx.Z(),event->ZDCx());
    // starCent = refmultCorrUtil->getCentralityBin9();
    // if (starCent < 0) continue; // bad centrality
    hEvt->Fill(eventType[7], 1);
    // if      (starCent == 8 || starCent == 7) icent = 0; // 0-10%
    // else if (starCent == 6 || starCent == 5) icent = 1; // 10-30%
    // else if (starCent == 4 || starCent == 3) icent = 2; // 30-50%
    // else if (starCent == 2 || starCent == 1 || starCent == 0) icent = 3; // 50-80%
    // else icent = -1;
    // if (icent<0) continue;
    // if (icent >= (int) nMyCent || icent<0) {
    //   cerr << "Error: icent = " << icent << endl;
    //   exit(0);
    // }


    if (bFXT) {
      // A dummy way to get the fxtMult
      iFxtMult = 0;
      // Track analysis
      nTracks = dst->numberOfTracks();    
      for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {
        StPicoTrack *picoTrack = dst->track(iTrk);
        if (!picoTrack) continue; 
        if (!picoTrack->isPrimary()) continue;
        iFxtMult++;
      }

      // Centrality cuts                                              3.2 GeV
      if      ( iFxtMult >  195 )    starCent = -1; // upper cut      287
      else if ( iFxtMult >= 142 )    starCent = 8;  // 0-5            205 
      else if ( iFxtMult >= 119 )    starCent = 7;  // 5-10           170 
      else if ( iFxtMult >= 86  )    starCent = 3;  // 10-20          118 
      else if ( iFxtMult >= 60  )    starCent = 5;  // 20-30          80  
      else if ( iFxtMult >= 41  )    starCent = 7;  // 30-40          53  
      else if ( iFxtMult >= 26  )    starCent = 9;  // 40-50          33  
      else if ( iFxtMult >= 16  )    starCent = 11; // 50-60          20  
      else if ( iFxtMult >= 9   )    starCent = 13; // 60-70          11  
      else if ( iFxtMult >= 5   )    starCent = 15; // 70-80          6   
      else                           starCent = -1;

      if      (starCent == 8 || starCent == 7) icent = 0; // 0-10%
      else if (starCent == 6 || starCent == 5) icent = 1; // 10-30%
      else if (starCent == 4 || starCent == 3) icent = 2; // 30-50%
      else if (starCent == 2 || starCent == 1 || starCent == 0) icent = 3; // 50-80%
      else icent = -1;
      if (icent<0) continue;
      if (icent >= (int) nMyCent || icent<0) {
        cerr << "Error: icent = " << icent << endl;
        exit(0);
      }
    }


    iVtxZ = -1;
    if      (nVtxZ == 1) iVtxZ = 0;
    else if (nVtxZ == 2) iVtxZ = pVtx.Z() >= 0. ? 1 : 0;
    else if (pVtx.Z() == 70.) { iVtxZ = nVtxZ - 1; }
    else {
      Double_t dVertexZShift = pVtx.Z() + 70.;
      Double_t dBinWidth = 2.*70./nVtxZ;
      iVtxZ = (int) dVertexZShift/dBinWidth;
    }      
    nTracks = dst->numberOfTracks();
    hVz->Fill(pVtx.Z());
    hVr->Fill(pVtx.X(),pVtx.Y());
    pion.clear();
    bField = event->bField();
    for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {
      StPicoTrack *picoTrack = dst->track(iTrk);
      StPicoBTofPidTraits *trait = dst->btofPidTraits(picoTrack->bTofPidTraitsIndex());
      if (!picoTrack)                                          continue;
      if (!picoTrack->isPrimary())                             continue;
      if ( ( picoTrack -> dEdx() ) == 0. )                     continue;
      if (picoTrack->nHitsFit() < cutNhits)                    continue;
      // if ((Double_t)picoTrack->nHitsFit() / picoTrack->nHitsPoss() < cutNhitsRatio) continue;
      // if (picoTrack->nHitsPoss() < cutNhitsPoss)              continue;
      if (picoTrack->gDCA(pVtx).Mag() > cutDCA.at(colEnergy))  continue;
      // Kinematic track cut
      TVector3 mom = picoTrack->pMom();
      if (!bFXT && TMath::Abs(mom.Eta()) > cutEta)             continue;
      if ( bFXT && (mom.Eta() < cutEtaFTX || mom.Eta() > 0.) ) continue;
      if (picoTrack->pPt() < cutPtMin.at(colEnergy) )          continue;
      if (picoTrack->pPt() > cutPtMax)                         continue;
      if (picoTrack->pPtot() < cutPtotMin)                     continue;
      if (picoTrack->pPtot() > cutPMax)                        continue;

      charge = picoTrack->charge();
      pid = (charge>0) ? (GetPID(picoTrack, trait) + 1) : (GetPID(picoTrack, trait) + 5);
      if (pid == 1) {
        TVector3 mom = picoTrack->pMom();
        ptot = mom.Mag();
        pt = mom.Pt();
        eta = mom.Eta();
        phi = mom.Phi();
        if (!picoTrack->isTofTrack()) hPiPtot_TPC_only->Fill(ptot);
        else                          hPiPtot_TPC_TOF->Fill(ptot);
        hPiPtot->Fill(ptot);
        if (picoTrack->isTofTrack()) hPiM2pq->Fill(ptot, TMath::Power(ptot,2.)/(1./(1.-TMath::Power(trait->btofBeta(),2.))-1.));
        hPidEdxpq->Fill(ptot, picoTrack->dEdx());
        hPiEtaPt->Fill(mom.Perp(),mom.Eta());
        // lrtzVecPi.SetXYZM(mom.Px(), mom.Py(), mom.Pz(), pion_mass);
        lrtzVecPi.SetPtEtaPhiM(pt, eta, phi, pion_mass);
        pion.push_back(myPart(lrtzVecPi, picoTrack->topologyMap(0), picoTrack->topologyMap(1), picoTrack->nHits(),picoTrack->helix(bField)));
      } // if (pid == 1)
      pt = mom.Pt();
      eta = mom.Eta();
      phi = mom.Phi();
      lrtzVecPi.SetPtEtaPhiM(pt, eta, phi, pion_mass);
      pion.push_back(myPart(lrtzVecPi, picoTrack->topologyMap(0), picoTrack->topologyMap(1), picoTrack->nHits(),picoTrack->helix(bField)));
    } // for(Int_t iTrk=0; iTrk<nTracks; iTrk++)

    if (pion.size()>=2) {
      for (UInt_t iPion1(0); iPion1<pion.size(); iPion1++) {
        for (UInt_t iPion2(iPion1+1); iPion2<pion.size(); iPion2++) {
          qinv = - (pion.at(iPion1).m4Mom - pion.at(iPion2).m4Mom).M();
          // qinv = TMath::Sqrt( - TMath::Power(pion.at(iPion1).m4Mom.Energy() - pion.at(iPion2).m4Mom.Energy(), 2.) + (pion.at(iPion1).m4Mom.Vect() - pion.at(iPion2).m4Mom.Vect()).Mag2() );
          // cout << "qinv = " << qinv << "; qinv = " << (pion.at(iPion1).m4Mom.Vect() - pion.at(iPion2).m4Mom.Vect()).Mag() << endl;
          hQinvNum[icent]->Fill(qinv);
          splittingLevel = GetSpittingLevel(pion[iPion1],pion[iPion2]);
          hSLNum[icent]->Fill(qinv,splittingLevel);
          for (UInt_t iSL=0;iSL<nSL;iSL++) if (splittingLevel < SL[iSL]) hQinvNum_SL[icent][iSL]->Fill(qinv);
          aveDist = GetAverageDistance(pion[iPion1],pion[iPion2]);
          hADNum[icent]->Fill(qinv,aveDist);
          pADNum[icent]->Fill(qinv,aveDist);
          hDeltaRNum[icent]->Fill(aveDist);
          for (UInt_t iSL=0;iSL<nSL;iSL++) if (splittingLevel < SL[iSL]) hDeltaRNum_SL[icent][iSL]->Fill(aveDist);
          for (UInt_t iAD=0;iAD<nAD;iAD++) if (aveDist > AD[iAD]) {
            hQinvNum_AD[icent][iAD]->Fill(qinv);
            if (splittingLevel < SL[1]) hQinvNum_SL_AD[icent][iAD]->Fill(qinv);
          } // for (UInt_t iAD=0;iAD<nAD;iAD++) if (aveDist > AD[iAD])
          kt = Getkt(pion.at(iPion1).m4Mom.Vect(),pion.at(iPion2).m4Mom.Vect());
          ikt = -1;
          for (UInt_t j(0);j<nkt;j++) { if (kt>=ktBin[j] && kt<ktBin[j+1]) {ikt=j; break;} }
          if (ikt!=-1) {
            hSLNum_kt[icent][ikt]->Fill(qinv,splittingLevel);
            hADNum_kt[icent][ikt]->Fill(qinv,aveDist);
            hQinvNum_kt[icent][ikt]->Fill(qinv);
            hDeltaRNum_kt[icent][ikt]->Fill(aveDist);
          } // if (ikt!=-1)
        } // for (UInt_t iPion2(0); iPion2<pionMom.size(); iPion2++)
      } // for (UInt_t iPion1(0); iPion1<pionMom.size(); iPion1++)
    } // if (pionMom.size()>=2)

    // ============================================================= //
    // ========================= Mix event ========================= //
    // ============================================================= //
    if (pion.size()>=1) { // There must be at least 1 pion to do the mixing
      for (const vector<myPart> & previousEvent : mixBuf[icent][iVtxZ]) { // loop over the buffer
        for (auto &mixPion : previousEvent) { // loop over the pions from previous event
          for (auto &pionFromThisEvent : pion) { // loop over the pions from this event
            qinv = -(pionFromThisEvent.m4Mom-mixPion.m4Mom).M();
            hQinvDen[icent]->Fill(qinv);
            splittingLevel = GetSpittingLevel(pionFromThisEvent, mixPion);
            hSLDen[icent]->Fill(qinv,splittingLevel);
            for (UInt_t iSL=0;iSL<nSL;iSL++) if (splittingLevel < SL[iSL]) hQinvDen_SL[icent][iSL]->Fill(qinv);
            aveDist = GetAverageDistance(pionFromThisEvent,mixPion);
            hADDen[icent]->Fill(qinv,aveDist);
            pADDen[icent]->Fill(qinv,aveDist);
            hDeltaRDen[icent]->Fill(aveDist);
            for (UInt_t iSL=0;iSL<nSL;iSL++) if (splittingLevel < SL[iSL]) hDeltaRDen_SL[icent][iSL]->Fill(aveDist);
            for (UInt_t iAD=0;iAD<nAD;iAD++) if (aveDist > AD[iAD]) {
              hQinvDen_AD[icent][iAD]->Fill(qinv);
              if (splittingLevel < SL[1]) hQinvDen_SL_AD[icent][iAD]->Fill(qinv);
            } // for (UInt_t iAD=0;iAD<nAD;iAD++) if (aveDist > AD[iAD])
            kt = Getkt(pionFromThisEvent.m4Mom.Vect(),mixPion.m4Mom.Vect());
            ikt = -1;
            for (UInt_t j(0);j<nkt;j++) { if (kt>=ktBin[j] && kt<ktBin[j+1]) {ikt=j; break;} }
            if (ikt!=-1) {
              hSLDen_kt[icent][ikt]->Fill(qinv,splittingLevel);
              hADDen_kt[icent][ikt]->Fill(qinv,aveDist);
              hQinvDen_kt[icent][ikt]->Fill(qinv);
              hDeltaRDen_kt[icent][ikt]->Fill(aveDist);
            } // if (ikt!=-1)
          } // for (auto &pionFromThisEvent : pion)
        } // for (auto &mixPion : previousEvent)
      } // for (const vector<myPart> & previousEvent : mixBuf[icent][iVtxZ])
      mixBuf[icent][iVtxZ].push_front(pion);
      if (mixBuf[icent][iVtxZ].size() > nBuffer) mixBuf[icent][iVtxZ].pop_back();
    } // if (pion.size()>=1) 

  } //for(Long64_t iEvent=0; iEvent<events2read; iEvent++)
  picoReader->Finish();
  // Writing output
  fo->cd();
  fo->Write();
  fo->Close();
  cout << "I'm done with analysis. We'll have a Nobel Prize, Master!" << endl;
  timer.Stop();
  timer.Print();
}

int main(int argc, char **argv)
{
  TString iFileName, oFileName;
  Float_t centerOfMassEnergy = 0.;

  if (argc < 7)
  {
    cerr << "./PicoDstFemtoscopy -i INPUT -o OUTPUT -cme colEnergy" << endl;
    return 1;
  }
  for (Int_t i = 1; i < argc; i++)
  {
    if (std::string(argv[i]) != "-i" &&
        std::string(argv[i]) != "-o" &&
        std::string(argv[i]) != "-cme")
    {
      cerr << "\n[ERROR]: Unknown parameter " << i << ": " << argv[i] << endl;
      return 2;
    }
    else
    {
      if (std::string(argv[i]) == "-i" && i != argc - 1)
      {
        iFileName = argv[++i];
        continue;
      }
      if (std::string(argv[i]) == "-i" && i == argc - 1)
      {
        cerr << "\n[ERROR]: Input file name was not specified " << endl;
        return 3;
      }
      if (std::string(argv[i]) == "-o" && i != argc - 1)
      {
        oFileName = argv[++i];
        continue;
      }
      if (std::string(argv[i]) == "-o" && i == argc - 1)
      {
        cerr << "\n[ERROR]: Output file name was not specified " << endl;
        return 4;
      }
      if (std::string(argv[i]) == "-cme" && i != argc - 1)
      {
        centerOfMassEnergy = atof(argv[++i]);
        continue;
      }
      if (std::string(argv[i]) == "-cme" && i == argc - 1)
      {
        cerr << "\n[ERROR]: Collision energy was not specified " << endl;
        return 1;
      }
    }
  }
  PicoDstFemtoscopy(iFileName, oFileName, centerOfMassEnergy);

  return 0;
}


Int_t  GetPID(StPicoTrack *const &track, StPicoBTofPidTraits *const &trait)
{
  Int_t i_part = -1;
  // TPC-only
  if (!track->isTofTrack()) {
    
    if (track->pPtot() >= 0.15 && track->pPtot() < 0.6 && TMath::Abs(track->nSigmaPion())  < cutNsigPID && TMath::Abs(track->nSigmaElectron()) > 2 ) i_part = 0; // pion id
    if (track->pPtot() >= 0.15 && track->pPtot() < 0.5 && TMath::Abs(track->nSigmaKaon())  < cutNsigPID ) i_part = 1; // kaon id
    if (track->pPtot() >= 0.4 && track->pPtot() < 0.9 && TMath::Abs(track->nSigmaProton()) < cutNsigPID ) i_part = 2; // proton id
    if (track->pPtot() >= 0.6 && track->pPtot() < 0.7 && TMath::Abs(track->nSigmaPion())   < cutNsigPID  && TMath::Abs(track->nSigmaKaon()) > 2) i_part = 0; // pion id
    if (track->pPtot() >= 0.5 && track->pPtot() < 0.7 && TMath::Abs(track->nSigmaKaon())   < cutNsigPID  && TMath::Abs(track->nSigmaPion()) > 3) i_part = 1; // kaon id
    if (track->pPtot() >= 0.9 && track->pPtot() < 1.2 && TMath::Abs(track->nSigmaProton()) < cutNsigPID  && TMath::Abs(track->nSigmaPion()) > 3) i_part = 2; // proton id
  } // if (!track->isTofTrack())
  // TPC+TOF
  if (isGoodPID(track,trait)) {
    Double_t ptot = track->pPtot();
    Double_t m2 = TMath::Power(ptot,2.)/(1./(1.-TMath::Power(trait->btofBeta(),2.))-1.);
    // pion id
    if (ptot >= 0.2 && ptot < 1.0 &&
        TMath::Abs(track->nSigmaPion()) < cutNsigPID &&
        m2 >= -0.15 && m2 < 0.1) {
      i_part = 0;
    }
    // kaon id
    if (ptot >= 0.2 && ptot < 1.0 &&
        TMath::Abs(track->nSigmaKaon()) < cutNsigPID &&
        m2 >= 0.2 && m2 < 0.32) {
      i_part = 1;
    }
    // proton id
    if (ptot >= 0.4 && ptot < 3.4 &&
        TMath::Abs(track->nSigmaProton()) < cutNsigPID &&
        m2 >= 0.74 && m2 < 1.2) {
      i_part = 2;
    }
  } // if (isGoodPID(track))
  return i_part;
}