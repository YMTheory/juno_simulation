void gen_root()
{

    ifstream in;
    in.open("PmtData_copy.csv");
    string line;
    int tmpcol, tmpid, tmpdyn, tmphqe;
    char tmpsn[20];
    double tmpgain, tmprsl, tmpdcr, tmptts, tmppde,  tmpttsss, tmphv, tmpamp, tmppvsv, tmpsvsn, tmprt, tmpft, tmpfwhm, tmptoff, tmpppp, tmpapp, tmpx, tmpy, tmpz;

    TFile* out = new TFile("PmtData_Lpmt.root", "recreate");
    TTree* tree = new TTree("PmtData_Lpmt", "lpmt parameters tree");
    tree->Branch("pmtID", &tmpid, "pmtID/I");
    tree->Branch("SN", &tmpsn, "SN[20]/C");
    tree->Branch("MCP_Hama", &tmpdyn, "MCP_Hama/I");
    tree->Branch("HiQE_MCP", &tmphqe, "HiQE_MCP/I");
    tree->Branch("Gain", &tmpgain, "Gain/D");
    tree->Branch("Resolution", &tmprsl, "Resolution/D");
    tree->Branch("PDE",&tmppde, "PDE/D");
    tree->Branch("DCR", &tmpdcr, "DCR/D");
    tree->Branch("TTS", &tmptts, "TTS/D");
    tree->Branch("TTS_SS", &tmpttsss, "TTS_SS/D");
    tree->Branch("HV", &tmphv, "HV/D");
    tree->Branch("Amplitude", &tmpamp, "Amplitude/D");
    tree->Branch("PvsV", &tmppvsv, "PvsV/D");
    tree->Branch("SvsN", &tmpsvsn, "SvsN/D");
    tree->Branch("RiseTime", &tmprt, "RiseTime/D");
    tree->Branch("FallTime", &tmpft, "FallTime/D");
    tree->Branch("FWHM", &tmpfwhm, "FWHM/D");
    tree->Branch("timeOffset", &tmptoff, "timeOffset/D");
    tree->Branch("prePulseProb", &tmpppp, "prePulseProb/D");
    tree->Branch("afterPulseProb", &tmpapp, "afterPulseProb/D");
    tree->Branch("pmtPosX", &tmpx, "pmtPosX/D");
    tree->Branch("pmtPosY", &tmpy, "pmtPosY/D");
    tree->Branch("pmtPosZ", &tmpz, "pmPostZ/D");



    while(getline(in, line)) {
        istringstream ss(line);
        ss >> tmpcol >> tmpid >> tmpsn >> tmpdyn >> tmphqe >> tmpgain >> tmprsl >>tmpdcr >> tmppde >> tmptts >> tmpamp >> tmphv >> tmppvsv >> tmpsvsn >> tmprt >> tmpft >> tmpfwhm >> tmpttsss >> tmptoff >> tmpppp >> tmpapp >> tmpx >> tmpy >> tmpz;
        tree->Fill();
    }

    tree->Write();
    out->Close();


}
