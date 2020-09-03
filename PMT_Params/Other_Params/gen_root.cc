void gen_root()
{

    ifstream in;
    in.open("PmtData_copy1.csv");
    string line;
    int tmpcol, tmpid, tmpdyn, tmphqe;
    char tmpsn[20];
    double tmpgain, tmprsl, tmpdcr, tmptts, tmppde,  tmpttsss, tmphv, tmpamp, tmppvsv, tmpsvsn, tmprt, tmpft, tmpfwhm, tmptoff, tmpppp, tmpapp, tmpx, tmpy, tmpz;

    TFile* out = new TFile("PmtData_Lpmt.root", "recreate");
    TTree* tree = new TTree("PmtData_Lpmt", "lpmt parameters tree");
    tree->Branch("pmtID", &tmpid, "pmtID/I");
    tree->Branch("SN", &tmpsn, "SN[20]/C");
    tree->Branch("isDyn", &tmpdyn, "isDyn/I");
    tree->Branch("isHqe", &tmphqe, "isHqe/I");
    tree->Branch("gain", &tmpgain, "gain/D");
    tree->Branch("resolution", &tmprsl, "resolution/D");
    tree->Branch("pde",&tmppde, "pde/D");
    tree->Branch("darkRate", &tmpdcr, "darkRate/D");
    tree->Branch("tts", &tmptts, "tts/D");
    tree->Branch("tts_ss", &tmpttsss, "tss_ss/D");
    tree->Branch("hv", &tmphv, "hv/D");
    tree->Branch("amplitude", &tmpamp, "amplitude/D");
    tree->Branch("PvsV", &tmppvsv, "PvsV/D");
    tree->Branch("SvsN", &tmpsvsn, "SvsN/D");
    tree->Branch("riseTime", &tmprt, "riseTime/D");
    tree->Branch("fallTime", &tmpft, "fallTime/D");
    tree->Branch("fwhm", &tmpfwhm, "fwhm/D");
    tree->Branch("timeOffset", &tmptoff, "timeOffset/D");
    tree->Branch("prePulseProb", &tmpppp, "prePulseProb/D");
    tree->Branch("afterPulseProb", &tmpapp, "afterPulseProb/D");
    tree->Branch("pmtX", &tmpx, "pmtX/D");
    tree->Branch("pmtY", &tmpy, "pmtY/D");
    tree->Branch("pmtZ", &tmpz, "pmtZ/D");



    while(getline(in, line)) {
        istringstream ss(line);
        ss >> tmpcol >> tmpid >> tmpsn >> tmpdyn >> tmphqe >> tmpgain >> tmprsl >>tmpdcr >> tmppde >> tmptts >> tmpamp >> tmphv >> tmppvsv >> tmpsvsn >> tmprt >> tmpft >> tmpfwhm >> tmpttsss >> tmptoff >> tmpppp >> tmpapp >> tmpx >> tmpy >> tmpz;
        tree->Fill();
    }

    tree->Write();
    out->Close();


}
