using namespace RooFit;                                                                                                           
using namespace std;

static const int nBins = 9;
float binBorders [nBins+1] = { 1, 2, 4.3, 6, 8.68, 10.09, 12.86, 14.18, 16, 19};

double PDGB0Mass = 5.27958;
double PDGJpsiMass = 3.096916;
double PDGPsiPrimeMass = 3.686109;

void Geteff(int year, int q2Bin, int binnum, int parity){
	  // year format: [6] for 2016
  //              [7] for 2017
  //              [8] for 2018
  //              [678] for total 3 years 
  // q2-bin format: [0-8] for one bin
    // flags to mark which q2 bins should be filled

cout << "year 201" << year << " q2Bin" << q2Bin << " binnum" << binnum << "parity" << parity << endl;
	bool runBin [nBins];
  string shortString [nBins];
  string longString  [nBins];
  for (int i=0; i<nBins; ++i) {
    runBin [i] = false;
    if ( q2Bin!=-1 && q2Bin!=i ) continue;
    runBin [i] = true;
    shortString [i] = Form("b%i",i);
    longString  [i] = Form("q2 bin %i",i);
  }
  // Load ntuples
  TChain* t_gen = new TChain();
  TChain* t_den = new TChain();
  TChain* t_num = new TChain();
  t_gen->Add("/eos/cms/store/user/fiorendi/p5prime/2016/skims/GEN_NoFilter/newphi/GEN_BFilter_B0MuMuKstar_p*.root/ntuple");
  if ( year==6 ) {
    // 2016
    t_den->Add("/eos/cms/store/user/fiorendi/p5prime/2016/skims/newphi/2016GEN_MC_LMNR.root/ntuple");
    t_num->Add("/eos/cms/store/user/fiorendi/p5prime/2016/skims/newphi/2016MC_LMNR.root/ntuple");
  } else if ( year==7 ) {
    // 2017
    t_den->Add("/eos/cms/store/user/fiorendi/p5prime/2017/skims/newphi/2017GEN_MC_LMNR.root/ntuple");
    t_num->Add("/eos/cms/store/user/fiorendi/p5prime/2017/skims/newphi/2017MC_LMNR.root/ntuple");
  } else if ( year==8 ) {
    // 2018
    t_den->Add("/eos/cms/store/user/fiorendi/p5prime/2018/skims/newphi/2018GEN_MC_LMNR.root/ntuple");
    t_num->Add("/eos/cms/store/user/fiorendi/p5prime/2018/skims/newphi/2018MC_LMNR.root/ntuple");
  } else if ( year==678 ) {
    // total 3 years
    t_den->Add("/eos/cms/store/user/fiorendi/p5prime/2016/skims/newphi/2016GEN_MC_LMNR.root/ntuple");
    t_den->Add("/eos/cms/store/user/fiorendi/p5prime/2017/skims/newphi/2017GEN_MC_LMNR.root/ntuple");
    t_den->Add("/eos/cms/store/user/fiorendi/p5prime/2018/skims/newphi/2018GEN_MC_LMNR.root/ntuple");
    t_num->Add("/eos/cms/store/user/fiorendi/p5prime/2016/skims/newphi/2016MC_LMNR.root/ntuple");
    t_num->Add("/eos/cms/store/user/fiorendi/p5prime/2017/skims/newphi/2017MC_LMNR.root/ntuple");
    t_num->Add("/eos/cms/store/user/fiorendi/p5prime/2018/skims/newphi/2018MC_LMNR.root/ntuple");
  }
  int genEntries = t_gen->GetEntries();
  int denEntries = t_den->GetEntries();
  int numEntries = t_num->GetEntries();


    // Import branches from ntuples:
  // angular variables
  double genCosThetaK, genCosThetaL, genPhi;
  double recoCosThetaK, recoCosThetaL, recoPhi;
  t_gen->SetBranchAddress( "cos_theta_k"     , &genCosThetaK  );
  t_gen->SetBranchAddress( "cos_theta_l"     , &genCosThetaL  );
  t_gen->SetBranchAddress( "phi_kst_mumu"    , &genPhi        );
  t_den->SetBranchAddress( "gen_cos_theta_k" , &genCosThetaK  );
  t_den->SetBranchAddress( "gen_cos_theta_l" , &genCosThetaL  );
  t_den->SetBranchAddress( "gen_phi_kst_mumu", &genPhi        );
  t_num->SetBranchAddress( "cos_theta_k"     , &recoCosThetaK );
  t_num->SetBranchAddress( "cos_theta_l"     , &recoCosThetaL );
  t_num->SetBranchAddress( "phi_kst_mumu"    , &recoPhi       );

  // variables for applying GEN-filter
  double genmupEta, genmumEta, genkstTrkpEta, genkstTrkmEta, genmupPt, genmumPt, genkstTrkpPt, genkstTrkmPt;
  t_gen->SetBranchAddress( "genmupEta", &genmupEta );
  t_gen->SetBranchAddress( "genmumEta", &genmumEta );
  t_gen->SetBranchAddress( "genkstTrkpEta", &genkstTrkpEta );
  t_gen->SetBranchAddress( "genkstTrkmEta", &genkstTrkmEta );
  t_gen->SetBranchAddress( "genmupPt", &genmupPt );
  t_gen->SetBranchAddress( "genmumPt", &genmumPt );
  t_gen->SetBranchAddress( "genkstTrkpPt", &genkstTrkpPt );
  t_gen->SetBranchAddress( "genkstTrkmPt", &genkstTrkmPt );

  // dimuon mass variables
  double genDimuMass2, recoDimuMass;
  t_gen->SetBranchAddress( "genQ2"   , &genDimuMass2 );
  t_den->SetBranchAddress( "genQ2"   , &genDimuMass2 );
  t_num->SetBranchAddress( "mumuMass", &recoDimuMass );


    // B0 mass variable
  double recoB0Mass;
  t_num->SetBranchAddress( "tagged_mass", &recoB0Mass );

  // B0-kinematic variables
  // double genB0pT, genB0eta;
  // double recoB0pT, recoB0eta;
  // t_gen->SetBranchAddress( "genbPt" , &genB0pT   );
  // t_gen->SetBranchAddress( "genbEta", &genB0eta  );
  // t_den->SetBranchAddress( "genbPt" , &genB0pT   );
  // t_den->SetBranchAddress( "genbEta", &genB0eta  );
  // t_num->SetBranchAddress( "bPt"    , &recoB0pT  );
  // t_num->SetBranchAddress( "bEta"   , &recoB0eta );

  // flavour tagging variables
  double genSignal, tagB0;
  t_num->SetBranchAddress( "genSignal", &genSignal );
  t_num->SetBranchAddress( "tagB0"    , &tagB0     );

  // event number for even/odd splitting
  double eventN_Dou;
  Long64_t eventN;
  t_gen->SetBranchAddress( "eventN", &eventN_Dou );
  t_den->SetBranchAddress( "eventN", &eventN     );
  t_num->SetBranchAddress( "eventN", &eventN     );

  // event pileup weight
  float PUweight = 1;
  t_den->SetBranchAddress( "weight", &PUweight );
  t_num->SetBranchAddress( "weight", &PUweight );

    // final state radiation flag
  double genSignHasFSR;
  t_gen->SetBranchAddress( "genSignHasFSR", &genSignHasFSR );

  string filestring = Form(parity==0 ? "1Deffb%iy%i_ev.root" : "1Deffb%iy%i_od.root",q2Bin,year);
  cout << "efficiency will be saved in " << filestring << endl;
  TFile *fout = new TFile(filestring.c_str(),"recreate");



  TH1D *f[9][9];
  string name[9] = {"f1s","m6s","m6c","f3","f4","f5","f7","f8","f9"};
  string kind[9] = {"genDen","genNum","den","CTreco","WTreco","CTeff","WTeff","CTweight","WTweight"};


  //int counter;
  double cos_theta_k,cos_theta_l,phi_kst_mumu;
  double f_1s,M_6s,M_6c,f_3,f_4,f_5,f_7,f_8,f_9;
  //Get binedges 

  //cout << "Get binedges....." << endl;
  vector<vector<double>> fi;

  for (int i=0;i<9;i++){
    fi.push_back(vector<double>());
  }

  double gennumentries=0;

  // Prepare GEN-level datasets
  cout<<"Starting adaptive-binning histogram setting..."<<endl;
  int xBin;
  for (int iCand=0; iCand<genEntries; ++iCand) {
    t_gen->GetEntry(iCand);
	  if (iCand%100000000==0){
	    cout << "entries" << iCand << endl;
    }
    cos_theta_k = genCosThetaK;
    cos_theta_l = genCosThetaL;
    phi_kst_mumu = genPhi;
    f_1s = 1 - cos_theta_k * cos_theta_k;
    M_6s = (1 - cos_theta_k * cos_theta_k) * cos_theta_l;
    M_6c = cos_theta_k * cos_theta_k * cos_theta_l;
    f_3 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Cos(2 * phi_kst_mumu);
    f_4 = 4 * cos_theta_k * cos_theta_l * TMath::Cos(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
    f_5 = 2 * cos_theta_k * TMath::Cos(phi_kst_mumu)* TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
    f_7 = 2 * cos_theta_k * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
    f_8 = 4 * cos_theta_k * cos_theta_l * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
    f_9 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Sin(2 * phi_kst_mumu);
    // find q2 bin of current candidate
	  if ( ( genDimuMass2 < binBorders[q2Bin+1] ) &&
        ( genDimuMass2 > binBorders[q2Bin]   ) ) {
			//if ( genSignHasFSR<0.5 ) {
      if ( fabs(genmupEta)<2.5 && fabs(genmumEta)<2.5 &&
          fabs(genkstTrkpEta)<2.5 && fabs(genkstTrkmEta)<2.5 &&
          genmupPt>2.5 && genmumPt>2.5 &&
          genkstTrkpPt>0.4 && genkstTrkmPt>0.4){
        //if (int(eventN_Dou)%2==parity){				
        gennumentries = gennumentries+1;
        fi[0].push_back(f_1s);
        fi[1].push_back(M_6s);
        fi[2].push_back(M_6c);
        fi[3].push_back(f_3);
        fi[4].push_back(f_4);
        fi[5].push_back(f_5);
        fi[6].push_back(f_7);
        fi[7].push_back(f_8);
        fi[8].push_back(f_9);
          //}
		    //}
		  }
    }
  }



  for (int i=0;i<9;i++){
    sort(fi[i].begin(),fi[i].end());
  }


  int countsPerBin=gennumentries/binnum;
  double binEdgef1s[binnum+1];
  double binEdgem6s[binnum+1];
  double binEdgem6c[binnum+1];
  double binEdgef3[binnum+1];
  double binEdgef4[binnum+1];
  double binEdgef5[binnum+1];
  double binEdgef7[binnum+1];
  double binEdgef8[binnum+1];
  double binEdgef9[binnum+1];

  binEdgef1s[0] = 0;
  cout << "f1s edge 0 is " << binEdgef1s[0] << endl;
  binEdgem6s[0] = -1;
  binEdgem6c[0] = -1;
  binEdgef3[0] = -1;
  binEdgef4[0] = -1;
  binEdgef5[0] = -1;
  binEdgef7[0] = -1;
  binEdgef8[0] = -1;
  binEdgef9[0] = -1;

  for (int edge=1;edge<binnum;edge++){
    binEdgef1s[edge] = fi[0].at(countsPerBin * edge);
    cout << "f1s edge " << edge << " is " <<binEdgef1s[edge] << endl;
    binEdgem6s[edge] = fi[1].at(countsPerBin * edge);
    binEdgem6c[edge] = fi[2].at(countsPerBin * edge);
    binEdgef3[edge] = fi[3].at(countsPerBin * edge);
    binEdgef4[edge] = fi[4].at(countsPerBin * edge);
    binEdgef5[edge] = fi[5].at(countsPerBin * edge);
    binEdgef7[edge] = fi[6].at(countsPerBin * edge);
    binEdgef8[edge] = fi[7].at(countsPerBin * edge);
    binEdgef9[edge] = fi[8].at(countsPerBin * edge);
  }

  binEdgef1s[binnum] = 1;
  cout << "f1s edge " << binnum << " is " << binEdgef1s[binnum] << endl;
  binEdgem6s[binnum] = 1;
  binEdgem6c[binnum] = 1;
  binEdgef3[binnum] = 1;
  binEdgef4[binnum] = 1;
  binEdgef5[binnum] = 1;
  binEdgef7[binnum] = 1;
  binEdgef8[binnum] = 1;
  binEdgef9[binnum] = 1;


  TH1D f1sstyle("f1sstyle","",binnum,binEdgef1s);
  TH1D m6sstyle("m6sstyle","",binnum,binEdgem6s);
  TH1D m6cstyle("m6cstyle","",binnum,binEdgem6c);
  TH1D f3style("f3style","",binnum,binEdgef3);
  TH1D f4style("f4style","",binnum,binEdgef4);
  TH1D f5style("f5style","",binnum,binEdgef5);
  TH1D f7style("f7style","",binnum,binEdgef7);
  TH1D f8style("f8style","",binnum,binEdgef8);
  TH1D f9style("f9style","",binnum,binEdgef9);


  for (int k=0;k<9;k++){
    f[0][k]=new TH1D(f1sstyle);
    f[1][k]=new TH1D(m6sstyle);
    f[2][k]=new TH1D(m6cstyle);
    f[3][k]=new TH1D(f3style);
    f[4][k]=new TH1D(f4style);
    f[5][k]=new TH1D(f5style);
    f[6][k]=new TH1D(f7style);
    f[7][k]=new TH1D(f8style);
    f[8][k]=new TH1D(f9style);
  }
  f[0][0]->SetName("haha");
  cout << "histo " << f[0][0]->GetName() << endl;


  for (int i=0;i<9;i++){
    for (int k=0;k<9;k++){
        f[i][k]->SetName(Form("%s%s",name[i].c_str(),kind[k].c_str()));
    }
  }
  // Prepare GEN-level datasets
  cout<<"Starting GEN datasets filling..."<<endl;
  //  int xBin;
  for (int iCand=0; iCand<genEntries; ++iCand) {
    t_gen->GetEntry(iCand);
    if (iCand%10000000==0){
      cout << "entries" << iCand << endl;
    }
    cos_theta_k = genCosThetaK;
    cos_theta_l = genCosThetaL;
    phi_kst_mumu = genPhi;
    f_1s = 1 - cos_theta_k * cos_theta_k;
    M_6s = (1 - cos_theta_k * cos_theta_k) * cos_theta_l;
    M_6c = cos_theta_k * cos_theta_k * cos_theta_l;
    f_3 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Cos(2 * phi_kst_mumu);
    f_4 = 4 * cos_theta_k * cos_theta_l * TMath::Cos(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
    f_5 = 2 * cos_theta_k * TMath::Cos(phi_kst_mumu)* TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));                                                                                                                       
    f_7 = 2 * cos_theta_k * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
    f_8 = 4 * cos_theta_k * cos_theta_l * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
    f_9 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Sin(2 * phi_kst_mumu);
    // find q2 bin of current candidate
    if ( ( genDimuMass2 < binBorders[q2Bin+1] ) &&
      ( genDimuMass2 > binBorders[q2Bin]) &&
      ( int(eventN_Dou)%2==parity)) {
      if ( genSignHasFSR<0.5 ) {
        f[0][0]->Fill(f_1s);
        f[1][0]->Fill(M_6s);
        f[2][0]->Fill(M_6c);
        f[3][0]->Fill(f_3);
        f[4][0]->Fill(f_4);
        f[5][0]->Fill(f_5);
        f[6][0]->Fill(f_7);
        f[7][0]->Fill(f_8);
        f[8][0]->Fill(f_9);      
      }

      if ( fabs(genmupEta)<2.5 && fabs(genmumEta)<2.5 &&
        fabs(genkstTrkpEta)<2.5 && fabs(genkstTrkmEta)<2.5 &&
        genmupPt>2.5 && genmumPt>2.5 &&
        genkstTrkpPt>0.4 && genkstTrkmPt>0.4){
        f[0][1]->Fill(f_1s);
        f[1][1]->Fill(M_6s);
        f[2][1]->Fill(M_6c);
        f[3][1]->Fill(f_3);
        f[4][1]->Fill(f_4);
        f[5][1]->Fill(f_5);
        f[6][1]->Fill(f_7);
        f[7][1]->Fill(f_8);
        f[8][1]->Fill(f_9); 
      }
    }
  }
 
  // Prepare denominator dataset
  cout<<"Starting denominator dataset filling..."<<endl;
  for (int iCand=0; iCand<denEntries; ++iCand) {
    t_den->GetEntry(iCand);
	  if (iCand%1000000==0){
	    cout << "entries" << iCand << endl;
	  }
	  cos_theta_k = genCosThetaK;
    cos_theta_l = genCosThetaL;
    phi_kst_mumu = genPhi;
    f_1s = 1 - cos_theta_k * cos_theta_k;
    M_6s = (1 - cos_theta_k * cos_theta_k) * cos_theta_l;
    M_6c = cos_theta_k * cos_theta_k * cos_theta_l;
    f_3 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Cos(2 * phi_kst_mumu);
    f_4 = 4 * cos_theta_k * cos_theta_l * TMath::Cos(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
    f_5 = 2 * cos_theta_k * TMath::Cos(phi_kst_mumu)* TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
    f_7 = 2 * cos_theta_k * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
    f_8 = 4 * cos_theta_k * cos_theta_l * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
    f_9 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Sin(2 * phi_kst_mumu);
    if ( ( genDimuMass2 < binBorders[q2Bin+1] ) &&
      ( genDimuMass2 > binBorders[q2Bin]   ) &&
      (eventN%2==parity)){
      f[0][2]->Fill(f_1s,PUweight);
      f[1][2]->Fill(M_6s,PUweight);
      f[2][2]->Fill(M_6c,PUweight);
      f[3][2]->Fill(f_3,PUweight);
      f[4][2]->Fill(f_4,PUweight);
      f[5][2]->Fill(f_5,PUweight);
      f[6][2]->Fill(f_7,PUweight);
      f[7][2]->Fill(f_8,PUweight);
      f[8][2]->Fill(f_9,PUweight);
	  }
  }

  // Prepare numerator dataset
  cout<<"Starting numerator dataset filling..."<<endl;
  for (int iCand=0; iCand<numEntries; ++iCand) {
    t_num->GetEntry(iCand);
    if (iCand%1000000==0){
	    cout << "entries" << iCand <<endl;
	  }
	  // anti-radiation cut
    if ( recoDimuMass < PDGJpsiMass ) { // below Jpsi
      if ( fabs( recoB0Mass - PDGB0Mass - recoDimuMass + PDGJpsiMass ) < 0.18 ) continue;
    } else if ( recoDimuMass > PDGPsiPrimeMass ) { // above PsiPrime
      if ( fabs( recoB0Mass - PDGB0Mass - recoDimuMass + PDGPsiPrimeMass ) < 0.08 ) continue;
    } else { // between the resonances
      if ( fabs( recoB0Mass - PDGB0Mass - recoDimuMass + PDGJpsiMass ) < 0.08 ) continue;
      if ( fabs( recoB0Mass - PDGB0Mass - recoDimuMass + PDGPsiPrimeMass ) < 0.09 ) continue;
    }
    cos_theta_k = recoCosThetaK;
    cos_theta_l = recoCosThetaL;
    phi_kst_mumu = recoPhi;
    f_1s = 1 - cos_theta_k * cos_theta_k;
    M_6s = (1 - cos_theta_k * cos_theta_k) * cos_theta_l;
    M_6c = cos_theta_k * cos_theta_k * cos_theta_l;
    f_3 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Cos(2 * phi_kst_mumu);
    f_4 = 4 * cos_theta_k * cos_theta_l * TMath::Cos(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
    f_5 = 2 * cos_theta_k * TMath::Cos(phi_kst_mumu)* TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
    f_7 = 2 * cos_theta_k * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
    f_8 = 4 * cos_theta_k * cos_theta_l * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
    f_9 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Sin(2 * phi_kst_mumu);

    if ( ( pow(recoDimuMass,2) < binBorders[q2Bin+1] ) &&
      ( pow(recoDimuMass,2) > binBorders[q2Bin]   ) &&
      (eventN%2==parity)){
		  if (genSignal != tagB0+1) { // correctly tagged events
        f[0][3]->Fill(f_1s,PUweight);
        f[1][3]->Fill(M_6s,PUweight);
        f[2][3]->Fill(M_6c,PUweight);
        f[3][3]->Fill(f_3,PUweight);
        f[4][3]->Fill(f_4,PUweight);
        f[5][3]->Fill(f_5,PUweight);
        f[6][3]->Fill(f_7,PUweight);
        f[7][3]->Fill(f_8,PUweight);
        f[8][3]->Fill(f_9,PUweight);
			}
			else{
        M_6s = -M_6s;
        M_6c = -M_6c;
        f_5 = -f_5;
        f_8 = -f_8;
        f_9 = -f_9;
        f[0][4]->Fill(f_1s,PUweight);
        f[1][4]->Fill(M_6s,PUweight);
        f[2][4]->Fill(M_6c,PUweight);
        f[3][4]->Fill(f_3,PUweight);
        f[4][4]->Fill(f_4,PUweight);
        f[5][4]->Fill(f_5,PUweight);
        f[6][4]->Fill(f_7,PUweight);
        f[7][4]->Fill(f_8,PUweight);
        f[8][4]->Fill(f_9,PUweight);
			}
	  }
  }

  double effC[9],effW[9];
  cout << "begin calculate efficiency" << endl;
  double genDen[9],genNum[9],den[9],ctRECO[9],wtRECO[9],cteff[9],ctefferror[9],wtefferror[9],wteff[9],ctweight[9],wtweight[9];
  for (int i=1;i<=binnum;i++){
    for (int k=0;k<9;k++){
      genDen[k] = f[k][0]->GetBinContent(i);

      genNum[k] = f[k][1]->GetBinContent(i);
      den[k] = f[k][2]->GetBinContent(i);

      ctRECO[k] = f[k][3]->GetBinContent(i);
      wtRECO[k] = f[k][4]->GetBinContent(i);
      cteff[k] = ctRECO[k]*genNum[k]/genDen[k]/den[k];
      wteff[k] = wtRECO[k]*genNum[k]/genDen[k]/den[k];
      ctefferror[k] = cteff[k]*(TMath::Sqrt((genDen[k]-genNum[k])/(genDen[k]*genNum[k]))+TMath::Sqrt((den[k]-ctRECO[k])/(den[k]*ctRECO[k])));
      wtefferror[k] = wteff[k]*(TMath::Sqrt((genDen[k]-genNum[k])/(genDen[k]*genNum[k]))+TMath::Sqrt((den[k]-wtRECO[k])/(den[k]*wtRECO[k])));
      if (cteff[k]!=0){
        ctweight[k]=1/cteff[k];
      }
      if (wteff[k]!=0){
        wtweight[k]=1/wteff[k];
      }


      f[k][5]->SetBinContent(i,cteff[k]);
      f[k][5]->SetBinError(i,ctefferror[k]);
      f[k][6]->SetBinContent(i,wteff[k]);
      f[k][6]->SetBinError(i,wtefferror[k]);
      f[k][7]->SetBinContent(i,ctweight[k]);

      f[k][8]->SetBinContent(i,wtweight[k]);

    }
  }
  for (int i=0;i<9;i++){
    for (int k=0;k<9;k++){
      f[i][k]->Write();
    }
  }



  cout << "end year " << year << " bin" << q2Bin << " parity"<< parity << endl;



}
	



