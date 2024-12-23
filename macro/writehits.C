{
    gROOT->Reset();
    TStopwatch timer;
    timer.Start();
    //   gROOT->LoadMacro("$VMCWORKDIR/gconfig/rootlogon.C");
    //   rootlogon();
    gSystem->Load("libSttMvdTracking");
    
    // MCpoints
    TFile filerun("FTF_7G0_1000_sim.root");
    TTree *treepnt = (TTree*) filerun.Get("pndsim");
    TClonesArray *pnt = new TClonesArray("PndSttPoint");
    treepnt->SetBranchAddress("STTPoint",&pnt);
    
    
    // Hits
    TFile filedigi("FTF_7G0_1000_digi.root");
    TTree *treedigi = (TTree*) filedigi.Get("pndsim");
    TClonesArray *digi = new TClonesArray("PndSttHit");
    treedigi->SetBranchAddress("STTHit",&digi);
    
    // HelixHits
    TFile filehelix("FTF_7G0_1000_reco.root");
    TTree *treereco = (TTree*) filehelix.Get("pndsim");
    TClonesArray *hh = new TClonesArray("PndSttHelixHit");
    treereco->SetBranchAddress("SttHelixHit",&hh);
    
    // Helix Tracks
    TClonesArray *track = new TClonesArray("PndTrackCand");
    treereco->SetBranchAddress("SttMvdGemTrackCand",&track);
    
    // event loop
    for (Int_t evt = 0; evt < 100 /*treereco->GetEntriesFast()*/; evt++) {
        
        treepnt->GetEntry(evt);
        treedigi->GetEntry(evt);
        treereco->GetEntry(evt);
        cout << "Event " << evt << ": MC Hits " << digi->GetEntriesFast() << endl;
 
        string name = "mct" + to_string(evt);
        ofstream fmct(name);
        for (Int_t m = 0; m < digi->GetEntriesFast(); m++) {
            PndSttPoint *mcpoint = (PndSttPoint*) pnt->At(m);
            if (mcpoint==0) continue;
            cout << m << ": x " << mcpoint->GetX() << " y " << mcpoint->GetY() << " z " << mcpoint->GetZ() << endl;
            fmct << mcpoint->GetX() << "   " << mcpoint->GetY() << "   " << mcpoint->GetZ() << endl;
        }
        fmct.close();
        
        string simname = "sim" + to_string(evt);
        ofstream fsim(simname);
        string diginame = "digi" + to_string(evt);
        ofstream fdigi(diginame);
        
        cout << "Event " << evt << ": num tracks " << track->GetEntriesFast() << endl;
        
        // tracks loop
        for (Int_t k = 0; k < track->GetEntriesFast(); k++) {
            
            PndTrackCand *stttrack = (PndTrackCand*) track->At(k);
            if(!stttrack) continue;
            
            cout << k << " Hits:" << stttrack->GetNHits() << endl;
            int mcid = stttrack->getMcTrackId();
            cout << "McTrackId " << mcid << endl;
            
            for (Int_t i = 0; i < stttrack->GetNHits(); i++) {
                PndTrackCandHit h = stttrack->GetSortedHit(i);
                Int_t hitindex = h.GetHitId();
                cout << " " << hitindex;
                PndSttHit* hit = (PndSttHit*) digi->At(hitindex);
                if (hit==0) continue;
                cout << " x " << hit->GetX() << " y " << hit->GetY() << " z " << hit->GetZ();
                fdigi << hit->GetX() << "   " << hit->GetY() << "   " << hit->GetZ() << endl;
                PndSttPoint *mcpoint = (PndSttPoint*) pnt->At(hit->GetRefIndex());
                if (mcpoint==0) continue;
                cout << "  MC: x " << mcpoint->GetX() << " y " << mcpoint->GetY() << " z " << mcpoint->GetZ() << endl;
                fsim << mcpoint->GetX() << "   " << mcpoint->GetY() << "   " << mcpoint->GetZ() << endl;
            }
            cout << endl;
            fsim.close();
            fdigi.close();
        }
    }
}
