#include <TTree.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TH1.h>
#include <Math/Functor.h>
#include <TMath.h>
#include <RooAbsPdf.h>
#include <RooWorkspace.h>
#include <RooDataSet.h>
#include <RooArgSet.h>
#include <RooArgList.h>
#include <RooRealVar.h>
#include <RooDataHist.h>
#include <RooHistFunc.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TSystem.h>
#include <TH3.h>
#include <TCanvas.h>
#include <TChain.h>

#include <sstream>
#include <time.h>
#include <ctime>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <tuple>

using std::cout;
using std::endl;
using std::vector;
using std::stringstream;
using std::string;
using namespace RooFit;
using namespace std;

// #############
// # Variables #
// #############
double q2Min = 0.;
double q2Max = 0.;
int q2Bin = -1;
double mumuMass2 = -1.;

double wight ;
double M_6s ;
double M_6c ;
double f_1s;
double f_3 ;
double f_4 ;
double f_5;
double f_6s ;
double f_6c;
double f_7 ;
double f_8 ;
double f_9 ;

double Mi ;
double frac[9] = {0.118,0.120,0.122,0.125,0.127,0.129,0.128,0.130,0.131};
//for 2017 all:2841914, mistag:368423
// ##############
// # Parameters #
// ##############
double value;
double error;
double Fl = 0.;
double AFB = 0.;
double S3 = 0.;
double S4 = 0.;
double S5 = 0.;
double S6 = 0.;
double S7 = 0.;
double S8 = 0.;
double S9 = 0.;

double P1 = 0.;
double P2 = 0.;
double P3 = 0.;
double P4p = 0.;
double P5p = 0.;
double P6p = 0.;
double P8p = 0.;

string paralist[15] = {"FlS", "AFBS", "S3S", "S4S", "S5S", "S7S", "S8S", "S9S", "P1S", "P2S", "P3S", "P4pS", "P5pS", "P6pS", "P8pS"};
//int bins = 100;


// ####################
// # angular varables #
// ####################
// #################
// # Reco-Level MC #
// #################
double cos_theta_k = 0;
TBranch *b_cos_theta_k = 0;

double cos_theta_l = 0;
TBranch *b_cos_theta_l = 0;

double phi_kst_mumu = 0;
TBranch *b_phi_kst_mumu = 0;

double mumuMass = 0;
TBranch *b_mumuMass = 0;

double tagB0 = 0;
TBranch *b_tagB0 = 0;

double genSignal = 0;
TBranch *b_genSignal = 0;


double recoB0Mass = 0;
TBranch *b_recoB0Mass = 0;

Long64_t eventN = 0;
TBranch *b_eventN;

double genSignHasFSR;
TBranch *b_genSignHasFSR;

double weight;

double PDGB0Mass = 5.27958;
double PDGJpsiMass = 3.096916;
double PDGPsiPrimeMass = 3.686109;

stringstream myString;

// ################
// # gen-Level MC #
// ################
//float cos_theta_k = 0;
//TBranch *b_cos_theta_k = 0;

//float cos_theta_l = 0;
//TBranch *b_cos_theta_l = 0;

//float phi_kst_mumu = 0;
//TBranch *b_phi_kst_mumu = 0;

double genQ = 0;
TBranch *b_genQ = 0;

void quzhi() {
	char a = '0' + q2Bin;
	switch (a) {
		case '0' :
			q2Min = 1.0;
			q2Max = 2.0;
			q2Bin = 0;
			break;

		case '1' :
			q2Min = 2.0;
			q2Max = 4.3;
			q2Bin = 1;
			break;

		case '2' :
			q2Min = 4.3;
			q2Max = 6.0;
			q2Bin = 2;
			break;

		case '3':
			q2Min = 6.0;
			q2Max = 8.68;
			q2Bin = 3;
			break;

			//        case '4':
			//            q2Min = 8.68;
			//            q2Max = 10.09;
			//            q2Bin = 4;
			//            break;

		case '5':
			q2Min = 10.09;
			q2Max = 12.86;
			q2Bin = 5;
			break;

			//        case '6':
			//            q2Min = 12.86;
			//            q2Max = 14.18;
			//            q2Bin = 6;
			//            break;

		case '7':
			q2Min = 14.18;
			q2Max = 16.0;
			q2Bin = 7;
			break;

			//        case '8':
			//            q2Min = 16.0;
			//            q2Max = 19.0;
			//            q2Bin = 8;
			//            break;

		default:
			break;
	}
}

//function to calculate all the values and errors of all the parameters
void CalValue(string Paratype, TH1D *f1s, TH1D *m6s, TH1D *m6c, TH1D *f3, TH1D *f4, TH1D *f5, TH1D *f7, TH1D *f8, TH1D *f9)
{

	if  (Paratype == "FlS") 
	{
		value = 2.0 - 2.5 * f1s->GetMean();
		error = 2.5 * f1s->GetRMS() /TMath::Sqrt(f1s->GetEntries() - 1);
	}
	else if  (Paratype == "AFBS")
	{
		value = 3.0/ 4.0 * (3.0 * m6s->GetMean() - 2.0 * m6c->GetMean());
		error = 3.0 / 4.0 * TMath::Sqrt(4.0 * TMath::Power(m6c->GetRMS(), 2) / (m6c->GetEntries() - 1) +
				9.0 * TMath::Power(m6s->GetRMS(), 2) / (m6s->GetEntries() - 1));
	}
	else if  (Paratype == "S3S")
	{
		value = 25.0 / 8.0 * f3->GetMean();
		error = 25.0 / 8.0 * f3->GetRMS() /TMath::Sqrt(f3->GetEntries() - 1);
	}
	else if  (Paratype == "S4S")
	{
		value = 25.0 / 8.0 * f4->GetMean();
		error = 25.0 / 8.0 * f4->GetRMS() /TMath::Sqrt(f4->GetEntries() - 1);
	}
	else if (Paratype == "S5S")
	{
		value = 2.5 * f5->GetMean();
		error = 2.5 * f5->GetRMS() /TMath::Sqrt(f5->GetEntries() - 1);
	}
	else if (Paratype == "S7S")
	{
		value = 2.5 * f7->GetMean();
		error = 2.5 * f7->GetRMS() /TMath::Sqrt(f7->GetEntries() - 1);
	}
	else if  (Paratype == "S8S")
	{
		value = 25.0 / 8.0 * f8->GetMean();
		error = 25.0 / 8.0 * f8->GetRMS() /TMath::Sqrt(f8->GetEntries() - 1);
	}
	else if  (Paratype == "S9S")
	{
		value = 25.0 / 8.0 * f9->GetMean();
		error = 25.0 / 8.0 * f9->GetRMS() /TMath::Sqrt(f9->GetEntries() - 1);
	}
	else if  (Paratype == "P1S")
	{
		value = 25./4. * f3->GetMean() / ( 1. - (2.0 - 2.5 * f1s->GetMean()) );
		error = TMath::Sqrt(TMath::Power(2/(1.0 - (2.0 - 2.5 * f1s->GetMean())) * (25/8 * f3->GetRMS() /TMath::Sqrt(f3->GetEntries() - 1)), 2 ) + TMath::Power( (25/4 * f3->GetMean() / TMath::Power( (1.0 - (2.0 - 2.5 * f1s->GetMean())), 2) ) * (2.5 * f1s->GetRMS() /TMath::Sqrt(f1s->GetEntries() - 1)), 2 ) );
	}
	else if  (Paratype == "P2S")
	{
		value = 2.5*0.5 * m6s->GetMean() / ( 1. - (2.0 - 2.5 * f1s->GetMean()) );
		error = 2.5*TMath::Sqrt( TMath::Power( 1/(2 * (1- (2.0 - 2.5 * f1s->GetMean()))) ,2) * (4 * TMath::Power(m6c->GetRMS(), 2) / (m6c->GetEntries() - 1) + 9 * TMath::Power(m6s->GetRMS(), 2) / (m6s->GetEntries() - 1)) + TMath::Power( (3.0 * m6s->GetMean() - 2.0 * m6c->GetMean())/( 2 * TMath::Power( (1-(2.0 - 2.5 * f1s->GetMean())) ,2 )) ,2) * (TMath::Power( 2.5 * f1s->GetRMS() ,2 ) / (f1s->GetEntries() - 1)) );
	}
	else if (Paratype == "P3S")
	{
		value = - 25.0/8.0 * f9->GetMean() / (1.0 - (2.0 - 2.5 * f1s->GetMean()));
		double a = TMath::Power( 2.5 * f1s->GetRMS() ,2 ) / (f1s->GetEntries() - 1);
		double b = TMath::Power( 25/8 * f9->GetRMS() ,2 ) / (f9->GetEntries() - 1);
		error = TMath::Sqrt(TMath::Power( 1.0 / (1.0 - (2.0 - 2.5 * f1s->GetMean())) ,2 ) * b + TMath::Power( 25/8 * f9->GetMean() / TMath::Power((1.0 - (2.0 - 2.5 * f1s->GetMean())) ,2) ,2) * a );
	}
	else if (Paratype == "P4pS")
	{
		value = 2*25/8 * f4->GetMean() /TMath::Sqrt((2.0 - 2.5 * f1s->GetMean()) * ( 1.0 - (2.0 - 2.5 * f1s->GetMean())) );
		double a = TMath::Power( 2.5 * f1s->GetRMS() ,2 ) / (f1s->GetEntries() - 1);
		double b = TMath::Power( 25/8 * f4->GetRMS() ,2 ) / (f4->GetEntries() - 1);
		error = 2*TMath::Sqrt( 1/ ((2.0 - 2.5 * f1s->GetMean()) * ( 1- (2.0 - 2.5 * f1s->GetMean()))) * b + TMath::Power((1 - 2 * (2.0 - 2.5 * f1s->GetMean())) * 25/8 * f4->GetMean() ,2) * a / ( 4 * TMath::Power((1 - (2.0 - 2.5 * f1s->GetMean())) ,3)) );
	}
	else if (Paratype == "P5pS")
	{
		double a = 2.0 - 2.5 * f1s->GetMean();
		value = 2.5 * f5->GetMean() /TMath::Sqrt(a * ( 1 - a) );
		double b = TMath::Power( 2.5 * f1s->GetRMS() ,2 ) / (f1s->GetEntries() - 1);
		double c = TMath::Power( 2.5 * f5->GetRMS() ,2 ) / (f5->GetEntries() - 1);
		error = TMath::Sqrt( 1/ (a * ( 1- a)) * c + TMath::Power((1 - 2 * a) * 2.5 * f5->GetMean() ,2) * b /( 4 * TMath::Power((1 - a) ,3)) );
	}
	else if (Paratype == "P6pS")
	{
		double a = 2.0 - 2.5 * f1s->GetMean();
		value = - 2.5 * f7->GetMean() /TMath::Sqrt(a * ( 1 - a) );
		double b = TMath::Power( 2.5 * f1s->GetRMS() ,2 ) / (f1s->GetEntries() - 1);
		double c = TMath::Power( 2.5 * f7->GetRMS() ,2 ) / (f7->GetEntries() - 1);
		error = TMath::Sqrt( 1/ (a * ( 1- a)) * c + TMath::Power((1 - 2 * a) * 2.5 * f7->GetMean() ,2) * b /(4 * TMath::Power((1 - a) ,3)) ); 
	}
	else if (Paratype == "P8pS")
	{
		double a = 2.0 - 2.5 * f1s->GetMean();
		value = 25/8 * f8->GetMean() /TMath::Sqrt(a * ( 1 - a) );
		double b = TMath::Power( 2.5 * f1s->GetRMS() ,2 ) / (f1s->GetEntries() - 1);
		double c = TMath::Power( 25/8 * f8->GetRMS() ,2 ) / (f8->GetEntries() - 1);
		error = TMath::Sqrt( 1/ (a * ( 1- a)) * c + TMath::Power((1 - 2 * a) * 25/8 * f8->GetMean() ,2) * b / ( 4 * TMath::Power((1 - a) ,3)) );
	}
	cout << "\n[Moment::GenCalValue]\t q2bin = " << q2Bin << ":" << Paratype << " = " << value << "+/- " << error  << endl;
}

// the code used to constrct histo from efficiency 
tuple<TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*,TH1D*> Constructhisto(int q2Bin){
   	stringstream myString;
	myString.clear();
	myString.str("");
	myString << "/afs/cern.ch/user/x/xuqin/work/B0KstMuMu/moment/1Deff/momentfinal/effnewphi/1Deffb" << q2Bin << "y6_ev.root" ;
	//any year any parity is okay as we only need the binning results
	cout << "\n[Moment::ReCalValue]\tTry to open " << myString.str().c_str() << endl;
	TFile *efffile = new TFile(myString.str().c_str(), "READ");
	if (efffile!=0){
		cout<< "file is existed" << endl;
	}
	TH1D *f1seffW = (TH1D*)efffile->Get("f1sWTeff");
	TH1D *m6seffW = (TH1D*)efffile->Get("m6sWTeff");
	TH1D *m6ceffW = (TH1D*)efffile->Get("m6cWTeff");
	TH1D *f3effW = (TH1D*)efffile->Get("f3WTeff");
	TH1D *f4effW = (TH1D*)efffile->Get("f4WTeff");
	TH1D *f5effW = (TH1D*)efffile->Get("f5WTeff");
	TH1D *f7effW = (TH1D*)efffile->Get("f7WTeff");
	TH1D *f8effW = (TH1D*)efffile->Get("f8WTeff");
	TH1D *f9effW = (TH1D*)efffile->Get("f9WTeff");


	int bins = f1seffW->GetNbinsX();
	//TAxis *axis = m6seffW->GetXaxis();
	double xbinsf1s[bins+1], xbinsm6s[bins+1], xbinsm6c[bins+1], xbinsf3[bins+1], xbinsf4[bins+1], xbinsf5[bins+1], xbinsf7[bins+1], xbinsf8[bins+1], xbinsf9[bins+1];
	for (int i=0;i<=bins;i++){
		xbinsf1s[i] = f1seffW->GetBinLowEdge(i+1);
		xbinsm6s[i] = m6seffW->GetBinLowEdge(i+1);
		xbinsm6c[i] = m6ceffW->GetBinLowEdge(i+1);
		xbinsf3[i] = f3effW->GetBinLowEdge(i+1);
		xbinsf4[i] = f4effW->GetBinLowEdge(i+1);
		xbinsf5[i] = f5effW->GetBinLowEdge(i+1);
		xbinsf7[i] = f7effW->GetBinLowEdge(i+1);
		xbinsf8[i] = f8effW->GetBinLowEdge(i+1);
		xbinsf9[i] = f9effW->GetBinLowEdge(i+1);
		//    cout << xbins[i] << endl;
	}


	TH1D *f1s = new TH1D("f1s","",bins,xbinsf1s);
	TH1D *m6s = new TH1D("m6s","",bins,xbinsm6s);
	TH1D *m6c = new TH1D("m6c","",bins,xbinsm6c);
	TH1D *f3 = new TH1D("f3","",bins,xbinsf3);
	TH1D *f4 = new TH1D("f4","",bins,xbinsf4);
	TH1D *f5 = new TH1D("f5","",bins,xbinsf5);
	TH1D *f7 = new TH1D("f7","",bins,xbinsf7);
	TH1D *f8 = new TH1D("f8","",bins,xbinsf8);
	TH1D *f9 = new TH1D("f9","",bins,xbinsf9); 
    return make_tuple(f1s,m6s,m6c,f3,f4,f5,f7,f8,f9);
}

void GenCalValue(int q2Bin, double q2Min, double q2Max) 
{
	value = 0.0;
	error= 0.0;
	M_6s = 0.;
	M_6c = 0.;
	f_1s = 0.;
	f_3 = 0.;
	f_4 = 0.;
	f_5 = 0.;
	f_6s = 0.;
	f_6c = 0.;
	f_7 = 0.;
	f_8 = 0.;
	f_9 = 0.;

    auto histo = Constructhisto(q2Bin);
	TH1D *f1s = get<0>(histo);
	TH1D *m6s = get<1>(histo);
	TH1D *m6c = get<2>(histo);
	TH1D *f3 = get<3>(histo);
	TH1D *f4 = get<4>(histo);
	TH1D *f5 = get<5>(histo);
	TH1D *f7 = get<6>(histo);
	TH1D *f8 = get<7>(histo);
	TH1D *f9 = get<8>(histo);

    string dataSetstring = "data_genDen";
    string dataSetstringodd = dataSetstring + Form("_od_b%i",q2Bin);
    string dataSetstringeven = dataSetstring + Form("_ev_b%i",q2Bin);

    // Load variables and dataset
    string filename = Form("effDataset_b%i_%i.root",q2Bin,2016);
    //any year is okay as they are the same for gen_nofilter sample
    TFile* fin = TFile::Open( filename.c_str() );
    if ( !fin || !fin->IsOpen() ) {
        cout<<"File not found: "<<filename<<endl;
        return;
    }
    RooWorkspace* wspodd = (RooWorkspace*)fin->Get(Form("ws_b%ip%i",q2Bin,1));
    if ( !wspodd || wspodd->IsZombie() ) {
        cout<<"Odd Workspace not found in file: "<<filename<<endl;
        return;
    }
    RooWorkspace* wspeven = (RooWorkspace*)fin->Get(Form("ws_b%ip%i",q2Bin,0));
    if ( !wspeven || wspeven->IsZombie() ) {
        cout<<"Even Workspace not found in file: "<<filename<<endl;
        return;
    }

    RooRealVar* ctK = wspodd->var("ctK");
    RooRealVar* ctL = wspodd->var("ctL");
    RooRealVar* phi = wspodd->var("phi");
    if ( !ctK || !ctL || !phi || ctK->IsZombie() || ctL->IsZombie() || phi->IsZombie() ) {
        cout<<"Variables not found in file: "<<filename<<endl;
        return;
    }

    RooDataSet* totdata = (RooDataSet*)wspodd->data(dataSetstringodd.c_str());
    if ( !totdata || totdata->IsZombie() ) {
        cout<<"Dataset "<<dataSetstringodd<<" not found in file: "<<filename<<endl;
        return;
    }
    cout << "total entries in odd data is " << totdata->numEntries() << endl;
    RooDataSet* totdata1 = (RooDataSet*)wspeven->data(dataSetstringeven.c_str());
    if ( !totdata1 || totdata1->IsZombie() ) {
        cout<<"Dataset "<<dataSetstringeven<<" not found in file: "<<filename<<endl;
        return;
    }
    cout << "total entries in even data is " << totdata1->numEntries() << endl;
    //merge odd and even data
    totdata->append(*totdata1);
    int totentries = totdata->numEntries();
    cout << "total entries in total data is " << totentries << endl;


    for (int i=0;i<totentries;i++){
        cos_theta_k = totdata->get(i)->getRealValue("ctK");
        cos_theta_l = totdata->get(i)->getRealValue("ctL");
        phi_kst_mumu = totdata->get(i)->getRealValue("phi");
        f_1s = 1 - cos_theta_k * cos_theta_k;
		M_6s = (1 - cos_theta_k * cos_theta_k) * cos_theta_l;
		M_6c = cos_theta_k * cos_theta_k * cos_theta_l;
		f_3 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Cos(2 * phi_kst_mumu);
		f_4 = 4 * cos_theta_k * cos_theta_l * TMath::Cos(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
		f_5 = 2 * cos_theta_k * TMath::Cos(phi_kst_mumu)* TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
		f_7 = 2 * cos_theta_k * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
		f_8 = 4 * cos_theta_k * cos_theta_l * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
		f_9 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Sin(2 * phi_kst_mumu);

        f1s->Fill(f_1s);
        m6s->Fill(M_6s);
        m6c->Fill(M_6c);
        f3->Fill(f_3);
        f4->Fill(f_4);
        f5->Fill(f_5);
        f7->Fill(f_7);
        f8->Fill(f_8);
        f9->Fill(f_9);  
    }

    for (int i =0;i<15;i++){
		string Paratype = paralist[i];
		CalValue(Paratype,f1s,m6s,m6c,f3,f4,f5,f7,f8,f9);
	}
}

void ReCalValue(int q2Bin, double q2Min, double q2Max, string TagType, int parity) {

	cout << "TagType = " << TagType;
	cout << "q2Bin = " << q2Bin;
	cout << "q2Min = " << q2Min;
	cout << "q2Max = " << q2Max;
	cout << "parity = " << parity;
	wight = 0.0;
	M_6s = 0.;
	M_6c = 0.;
	f_1s = 0.;
	f_3 = 0.;
	f_4 = 0.;
	f_5 = 0.;
	f_6s = 0.;
	f_6c = 0.;
	f_7 = 0.;
	f_8 = 0.;
    f_9 = 0.;
    auto histo = Constructhisto(q2Bin);
	TH1D *f1s = get<0>(histo);
	TH1D *m6s = get<1>(histo);
	TH1D *m6c = get<2>(histo);
	TH1D *f3 = get<3>(histo);
	TH1D *f4 = get<4>(histo);
	TH1D *f5 = get<5>(histo);
	TH1D *f7 = get<6>(histo);
	TH1D *f8 = get<7>(histo);
	TH1D *f9 = get<8>(histo);

    for (int year=6;year<9;year++){
        // ############################
		// # Read Effcicency Function #
		// ############################
		cout << "\n[Moment::ReCalValue]\tTry to fill year" << year << "reco mc sample" << endl; 
        myString.clear();
		myString.str("");
		//myString << "/afs/cern.ch/user/x/xuqin/work/B0KstMuMu/moment/1Deff/momentfinal/eff/1Deffb" << q2Bin << "y"<< year <<".root" ;
		myString << "/afs/cern.ch/user/x/xuqin/work/B0KstMuMu/moment/1Deff/momentfinal/effnewphi/1Deffb" << q2Bin << "y" << year << "_" << (parity==0? "od":"ev") << ".root" ;
		//odd efficiency for even event; even efficiency for odd event
		cout << "\n[Moment::ReCalValue]\tTry to open" << myString.str().c_str() << endl;
		TFile *efffile = new TFile(myString.str().c_str(), "READ");
		if (efffile!=0){
			cout<< "file is existed" << endl;
		}
		TH1D *f1seffC = (TH1D*)efffile->Get("f1sCTeff");
		TH1D *m6seffC = (TH1D*)efffile->Get("m6sCTeff");
		TH1D *m6ceffC = (TH1D*)efffile->Get("m6cCTeff");
		TH1D *f3effC = (TH1D*)efffile->Get("f3CTeff");
		TH1D *f4effC = (TH1D*)efffile->Get("f4CTeff");
		TH1D *f5effC = (TH1D*)efffile->Get("f5CTeff");
		TH1D *f7effC = (TH1D*)efffile->Get("f7CTeff");
		TH1D *f8effC = (TH1D*)efffile->Get("f8CTeff");
		TH1D *f9effC = (TH1D*)efffile->Get("f9CTeff");

		TH1D *f1seffW = (TH1D*)efffile->Get("f1sWTeff");
		TH1D *m6seffW = (TH1D*)efffile->Get("m6sWTeff");
		TH1D *m6ceffW = (TH1D*)efffile->Get("m6cWTeff");
		TH1D *f3effW = (TH1D*)efffile->Get("f3WTeff");
		TH1D *f4effW = (TH1D*)efffile->Get("f4WTeff");
		TH1D *f5effW = (TH1D*)efffile->Get("f5WTeff");
		TH1D *f7effW = (TH1D*)efffile->Get("f7WTeff");
		TH1D *f8effW = (TH1D*)efffile->Get("f8WTeff");
		TH1D *f9effW = (TH1D*)efffile->Get("f9WTeff");

		RooRealVar f1se("f1se","f1se",0,1);
		RooRealVar m6se("m6se","m6se",-1,1);
		RooRealVar m6ce("m6ce","m6ce",-1,1);
		RooRealVar f3e("f3e","f3e",-1,1);
		RooRealVar f4e("f4e","f4e",-1,1);
		RooRealVar f5e("f5e","f5e",-1,1);
		RooRealVar f7e("f7e","f7e",-1,1);
		RooRealVar f8e("f8e","f8e",-1,1);
		RooRealVar f9e("f9e","f9e",-1,1);

		RooDataHist* f1sCData = new RooDataHist("f1sCData","f1sCData",f1se,Import(*f1seffC));
		RooDataHist* m6sCData = new RooDataHist("m6sCData","m6sCData",m6se,Import(*m6seffC));
		RooDataHist* m6cCData = new RooDataHist("m6cCData","m6cCData",m6ce,Import(*m6ceffC));
		RooDataHist* f3CData = new RooDataHist("f3CData","f3CData",f3e,Import(*f3effC));
		RooDataHist* f4CData = new RooDataHist("f4CData","f4CData",f4e,Import(*f4effC));
		RooDataHist* f5CData = new RooDataHist("f5CData","f5CData",f5e,Import(*f5effC));
		RooDataHist* f7CData = new RooDataHist("f7CData","f7CData",f7e,Import(*f7effC));
		RooDataHist* f8CData = new RooDataHist("f8CData","f8CData",f8e,Import(*f8effC));
		RooDataHist* f9CData = new RooDataHist("f9CData","f9CData",f9e,Import(*f9effC));

		RooAbsReal* effCf1s = new RooHistFunc("effCf1s","effCf1s",f1se,*f1sCData,1);
		RooAbsReal* effCm6s = new RooHistFunc("effCm6s","effCm6s",m6se,*m6sCData,1);
		RooAbsReal* effCm6c = new RooHistFunc("effCm6c","effCm6c",m6ce,*m6cCData,1);
		RooAbsReal* effCf3 = new RooHistFunc("effCf3","effCf3",f3e,*f3CData,1);
		RooAbsReal* effCf4 = new RooHistFunc("effCf4","effCf4",f4e,*f4CData,1);
		RooAbsReal* effCf5 = new RooHistFunc("effCf5","effCf5",f5e,*f5CData,1);
		RooAbsReal* effCf7 = new RooHistFunc("effCf7","effCf7",f7e,*f7CData,1);
		RooAbsReal* effCf8 = new RooHistFunc("effCf8","effCf8",f8e,*f8CData,1);
		RooAbsReal* effCf9 = new RooHistFunc("effCf9","effCf9",f9e,*f9CData,1);

		RooDataHist* f1sWData = new RooDataHist("f1sWData","f1sWData",f1se,Import(*f1seffW));
		RooDataHist* m6sWData = new RooDataHist("m6sWData","m6sWData",m6se,Import(*m6seffW));
		RooDataHist* m6cWData = new RooDataHist("m6cWData","m6cWData",m6ce,Import(*m6ceffW));
		RooDataHist* f3WData = new RooDataHist("f3WData","f3WData",f3e,Import(*f3effW));
		RooDataHist* f4WData = new RooDataHist("f4WData","f4WData",f4e,Import(*f4effW));
		RooDataHist* f5WData = new RooDataHist("f5WData","f5WData",f5e,Import(*f5effW));
		RooDataHist* f7WData = new RooDataHist("f7WData","f7WData",f7e,Import(*f7effW));
		RooDataHist* f8WData = new RooDataHist("f8WData","f8WData",f8e,Import(*f8effW));
		RooDataHist* f9WData = new RooDataHist("f9WData","f9WData",f9e,Import(*f9effW));

		RooAbsReal* effWf1s = new RooHistFunc("effWf1s","effWf1s",f1se,*f1sWData,1);
		RooAbsReal* effWm6s = new RooHistFunc("effWm6s","effWm6s",m6se,*m6sWData,1);
		RooAbsReal* effWm6c = new RooHistFunc("effWm6c","effWm6c",m6ce,*m6cWData,1);
		RooAbsReal* effWf3 = new RooHistFunc("effWf3","effWf3",f3e,*f3WData,1);
		RooAbsReal* effWf4 = new RooHistFunc("effWf4","effWf4",f4e,*f4WData,1);
		RooAbsReal* effWf5 = new RooHistFunc("effWf5","effWf5",f5e,*f5WData,1);
		RooAbsReal* effWf7 = new RooHistFunc("effWf7","effWf7",f7e,*f7WData,1);
		RooAbsReal* effWf8 = new RooHistFunc("effWf8","effWf8",f8e,*f8WData,1);
		RooAbsReal* effWf9 = new RooHistFunc("effWf9","effWf9",f9e,*f9WData,1);

        double wightf1s,wightm6s,wightm6c,wightf3,wightf4,wightf5,wightf7,wightf8,wightf9;
		double wight1f1s,wight1m6s,wight1m6c,wight1f3,wight1f4,wight1f5,wight1f7,wight1f8,wight1f9;
		double wight2f1s,wight2m6s,wight2m6c,wight2f3,wight2f4,wight2f5,wight2f7,wight2f8,wight2f9;
        
        string dataSetstringCT = "data_ctRECO";
        string dataSetstringWT = "data_wtRECO";
        dataSetstringCT = dataSetstringCT + Form((parity==0?"_ev_b%i":"_od_b%i"),q2Bin);
        dataSetstringWT = dataSetstringWT + Form((parity==0?"_ev_b%i":"_od_b%i"),q2Bin);


        // Load variables and dataset
        string filename = Form("effDataset_b%i_201%i.root",q2Bin,year);
        TFile* fin = TFile::Open( filename.c_str() );
        if ( !fin || !fin->IsOpen() ) {
            cout<<"File not found: "<<filename<<endl;
            return;
        }
        RooWorkspace* wsp = (RooWorkspace*)fin->Get(Form("ws_b%ip%i",q2Bin,parity));
        if ( !wsp || wsp->IsZombie() ) {
            cout<<"Workspace not found in file: "<<filename<<endl;
        return;
        }
        RooRealVar* ctK = wsp->var("ctK");
        RooRealVar* ctL = wsp->var("ctL");
        RooRealVar* phi = wsp->var("phi");
        if ( !ctK || !ctL || !phi || ctK->IsZombie() || ctL->IsZombie() || phi->IsZombie() ) {
            cout<<"Variables not found in file: "<<filename<<endl;
            return;
        }

        if  (TagType == "good") {
            RooDataSet* totdata = (RooDataSet*)wsp->data(dataSetstringCT.c_str());
            if ( !totdata || totdata->IsZombie() ) {
                cout<<"Dataset "<<dataSetstringCT<<" not found in file: "<<filename<<endl;
                return;
            }
            int totentries = totdata->numEntries();
            cout << "total entries in total data is " << totentries << endl;
            for (int i=0;i<totentries;i++){
                cos_theta_k = totdata->get(i)->getRealValue("ctK");
                cos_theta_l = totdata->get(i)->getRealValue("ctL");
                phi_kst_mumu = totdata->get(i)->getRealValue("phi");
                weight = totdata->weight();
                f_1s = 1 - cos_theta_k * cos_theta_k;
                M_6s = (1 - cos_theta_k * cos_theta_k) * cos_theta_l;
                M_6c = cos_theta_k * cos_theta_k * cos_theta_l;
                f_3 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Cos(2 * phi_kst_mumu);
                f_4 = 4 * cos_theta_k * cos_theta_l * TMath::Cos(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
                f_5 = 2 * cos_theta_k * TMath::Cos(phi_kst_mumu)* TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
                f_7 = 2 * cos_theta_k * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
                f_8 = 4 * cos_theta_k * cos_theta_l * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
                f_9 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Sin(2 * phi_kst_mumu);

                f1se.setVal(f_1s);
                double efff1s = effCf1s->getValV();
                m6se.setVal(M_6s);
                double effm6s = effCm6s->getValV();
                m6ce.setVal(M_6c);
                double effm6c = effCm6c->getValV();
                f3e.setVal(f_3);
                double efff3 = effCf3->getValV();
                f4e.setVal(f_4);
                double efff4 = effCf4->getValV();
                f5e.setVal(f_5);
                double efff5 = effCf5->getValV();
                f7e.setVal(f_7);
                double efff7 = effCf7->getValV();
                f8e.setVal(f_8);
                double efff8 = effCf8->getValV();
                f9e.setVal(f_9);
                double efff9 = effCf9->getValV();
                if (efff1s ==0)                                                                                                       
                    wightf1s = 0;
                else
                    wightf1s = 1/efff1s;
                if (effm6s ==0)
                    wightm6s = 0;
                else 
                    wightm6s = 1/effm6s;
                if (effm6c ==0)
                    wightm6c = 0;
                else 
                    wightm6c = 1/effm6c;
                if (efff3 ==0)
                    wightf3 = 0;
                else 
                    wightf3 = 1/efff3;
                if (efff4 ==0)
                    wightf4 = 0;
                else        
                    wightf4 = 1/efff4;
                if (efff5 ==0)
                    wightf5 = 0;
                else        
                    wightf5 = 1/efff5;
                if (efff7 ==0)
                    wightf7 = 0;
                else        
                    wightf7 = 1/efff7;
                if (efff8 ==0)
                    wightf8 = 0;
                else        
                    wightf8 = 1/efff8;
                if (efff9 ==0)
                    wightf9 = 0;
                else        
                    wightf9 = 1/efff9;

                f1s->Fill(f_1s, wightf1s*weight);
                m6s->Fill(M_6s, wightm6s*weight);
                m6c->Fill(M_6c, wightm6c*weight);
                f3->Fill(f_3, wightf3*weight);
                f4->Fill(f_4, wightf4*weight);
                f5->Fill(f_5, wightf5*weight);
                f7->Fill(f_7, wightf7*weight);
                f8->Fill(f_8, wightf8*weight);
                f9->Fill(f_9, wightf9*weight);
            }
        }

        else if  (TagType == "mis") {
            RooDataSet* totdata = (RooDataSet*)wsp->data(dataSetstringWT.c_str());
            if ( !totdata || totdata->IsZombie() ) {
                cout<<"Dataset "<<dataSetstringWT<<" not found in file: "<<filename<<endl;
                return;
            }
            int totentries = totdata->numEntries();
            cout << "total entries in total data is " << totentries << endl;
            for (int i=0;i<totentries;i++){
                cos_theta_k = totdata->get(i)->getRealValue("ctK");
                cos_theta_l = totdata->get(i)->getRealValue("ctL");
                phi_kst_mumu = totdata->get(i)->getRealValue("phi");
                weight = totdata->weight();
                f_1s = 1 - cos_theta_k * cos_theta_k;
                M_6s = -(1 - cos_theta_k * cos_theta_k) * cos_theta_l;
                M_6c = -cos_theta_k * cos_theta_k * cos_theta_l;
                f_3 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Cos(2 * phi_kst_mumu);
                f_4 = 4 * cos_theta_k * cos_theta_l * TMath::Cos(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
                f_5 = -2 * cos_theta_k * TMath::Cos(phi_kst_mumu)* TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
                f_7 = 2 * cos_theta_k * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
                f_8 = -4 * cos_theta_k * cos_theta_l * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
                f_9 = -(1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Sin(2 * phi_kst_mumu);
                f1se.setVal(f_1s);
                double efff1s = effWf1s->getValV();
                m6se.setVal(M_6s);
                double effm6s = effWm6s->getValV();
                m6ce.setVal(M_6c);
                double effm6c = effWm6c->getValV();
                f3e.setVal(f_3);
                double efff3 = effWf3->getValV();
                f4e.setVal(f_4);
                double efff4 = effWf4->getValV();
                f5e.setVal(f_5);
                double efff5 = effWf5->getValV();
                f7e.setVal(f_7);
                double efff7 = effWf7->getValV();
                f8e.setVal(f_8);
                double efff8 = effWf8->getValV();
                f9e.setVal(f_9);
                double efff9 = effWf9->getValV();
                if (efff1s ==0)                                                                                                       
                    wightf1s = 0;
                else
                    wightf1s = 1/efff1s;
                if (effm6s ==0)
                    wightm6s = 0;
                else 
                    wightm6s = 1/effm6s;
                if (effm6c ==0)
                    wightm6c = 0;
                else 
                    wightm6c = 1/effm6c;
                if (efff3 ==0)
                    wightf3 = 0;
                else 
                    wightf3 = 1/efff3;
                if (efff4 ==0)
                    wightf4 = 0;
                else        
                    wightf4 = 1/efff4;
                if (efff5 ==0)
                    wightf5 = 0;
                else        
                    wightf5 = 1/efff5;
                if (efff7 ==0)
                    wightf7 = 0;
                else        
                    wightf7 = 1/efff7;
                if (efff8 ==0)
                    wightf8 = 0;
                else        
                    wightf8 = 1/efff8;
                if (efff9 ==0)
                    wightf9 = 0;
                else        
                    wightf9 = 1/efff9;

                f1s->Fill(f_1s, wightf1s*weight);
                m6s->Fill(M_6s, wightm6s*weight);
                m6c->Fill(M_6c, wightm6c*weight);
                f3->Fill(f_3, wightf3*weight);
                f4->Fill(f_4, wightf4*weight);
                f5->Fill(f_5, wightf5*weight);
                f7->Fill(f_7, wightf7*weight);
                f8->Fill(f_8, wightf8*weight);
                f9->Fill(f_9, wightf9*weight);
            }
        }

        else if (TagType == "mix"){
            RooDataSet* totdata = (RooDataSet*)wsp->data(dataSetstringCT.c_str());
            if ( !totdata || totdata->IsZombie() ) {
                cout<<"Dataset "<<dataSetstringCT<<" not found in file: "<<filename<<endl;
                return;
            }
            cout << "goodtag events number " << totdata->numEntries() << endl;
            RooDataSet* totdata1 = (RooDataSet*)wsp->data(dataSetstringWT.c_str());
            if ( !totdata1 || totdata1->IsZombie() ) {
                cout<<"Dataset "<<dataSetstringWT<<" not found in file: "<<filename<<endl;
                return;
            }
            cout << "mistag events number " << totdata->numEntries() << endl;
            totdata->append(*totdata1);

            int totentries = totdata->numEntries();
            cout << "total entries in total data is " << totentries << endl;
            for (int i=0;i<totentries;i++){
                cos_theta_k = totdata->get(i)->getRealValue("ctK");
                cos_theta_l = totdata->get(i)->getRealValue("ctL");
                phi_kst_mumu = totdata->get(i)->getRealValue("phi");
                weight = totdata->weight();            
                f_1s = 1 - cos_theta_k * cos_theta_k;
                M_6s = (1 - cos_theta_k * cos_theta_k) * cos_theta_l;
                M_6c = cos_theta_k * cos_theta_k * cos_theta_l;
                f_3 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Cos(2 * phi_kst_mumu);
                f_4 = 4 * cos_theta_k * cos_theta_l * TMath::Cos(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
                f_5 = 2 * cos_theta_k * TMath::Cos(phi_kst_mumu)* TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
                f_7 = 2 * cos_theta_k * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
                f_8 = 4 * cos_theta_k * cos_theta_l * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
                f_9 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Sin(2 * phi_kst_mumu);	

                f1se.setVal(f_1s);
                double effc1f1s = effCf1s->getValV();
                double effw1f1s = effWf1s->getValV();
                m6se.setVal(M_6s);
                double effc1m6s = effCm6s->getValV();
                double effw1m6s = effWm6s->getValV();
                m6ce.setVal(M_6c);
                double effc1m6c = effCm6c->getValV();
                double effw1m6c = effWm6c->getValV();
                f3e.setVal(f_3);
                double effc1f3 = effCf3->getValV();
                double effw1f3 = effWf3->getValV();
                f4e.setVal(f_4);
                double effc1f4 = effCf4->getValV();
                double effw1f4 = effWf4->getValV();
                f5e.setVal(f_5);
                double effc1f5 = effCf5->getValV();
                double effw1f5 = effWf5->getValV();
                f7e.setVal(f_7);
                double effc1f7 = effCf7->getValV();
                double effw1f7 = effWf7->getValV();
                f8e.setVal(f_8);
                double effc1f8 = effCf8->getValV();
                double effw1f8 = effWf8->getValV();
                f9e.setVal(f_9);
                double effc1f9 = effCf9->getValV();
                double effw1f9 = effWf9->getValV();

                m6se.setVal(-M_6s);
                double effc2m6s = effCm6s->getValV();
                double effw2m6s = effWm6s->getValV();
                m6ce.setVal(-M_6c);
                double effc2m6c = effCm6c->getValV();
                double effw2m6c = effWm6c->getValV();
                f5e.setVal(-f_5);
                double effc2f5 = effCf5->getValV();
                double effw2f5 = effWf5->getValV();
                f8e.setVal(-f_8);
                double effc2f8 = effCf8->getValV();
                double effw2f8 = effWf8->getValV();
                f9e.setVal(-f_9);
                double effc2f9 = effCf9->getValV();
                double effw2f9 = effWf9->getValV();


                double efff1s = effc1f1s+effw1f1s;
                double efff3 = effc1f3+effw1f3;
                double efff4 = effc1f4+effw1f4;
                double efff7 = effc1f7+effw1f7;
                if (efff1s ==0)
                    wightf1s = 0;
                else
                    wightf1s = 1/efff1s;
                if (efff3 ==0)
                    wightf3 = 0;
                else 
                    wightf3 = 1/efff3;
                if (efff4 ==0)
                    wightf4 = 0;
                else    
                    wightf4 = 1/efff4;
                if (efff7 ==0)
                    wightf7 = 0;
                else    
                    wightf7 = 1/efff7;

                double denom6s = effc1m6s*effc2m6s-effw1m6s*effw2m6s;
                if (denom6s==0){
                    wight1m6s=0;
                    wight2m6s=0;
                }
                else{
                    wight1m6s = effc2m6s/denom6s;
                    wight2m6s = -effw1m6s/denom6s;
                }
                double denom6c = effc1m6c*effc2m6c-effw1m6c*effw2m6c;
                if (denom6c==0){
                    wight1m6c=0;
                    wight2m6c=0;
                }
                else{
                    wight1m6c = effc2m6c/denom6c;
                    wight2m6c = -effw1m6c/denom6c;
                }


                double denof5 = effc1f5*effc2f5-effw1f5*effw2f5;
                if (denof5==0){
                    wight1f5=0;
                    wight2f5=0;
                }
                else{
                    wight1f5 = effc2f5/denof5;
                    wight2f5 = -effw1f5/denof5;
                }

                double denof8 = effc1f8*effc2f8-effw1f8*effw2f8;
                if (denof8==0){
                    wight1f8=0;
                    wight2f8=0;
                }
                else{
                    wight1f8 = effc2f8/denof8;
                    wight2f8 = -effw1f8/denof8;
                }
                double denof9 = effc1f9*effc2f9-effw1f9*effw2f9;
                if (denof9==0){
                    wight1f9=0;
                    wight2f9=0;
                }
                else{
                    wight1f9 = effc2f9/denof9;
                    wight2f9 = -effw1f9/denof9;
                }
                f1s->Fill(f_1s, wightf1s*weight);
                f3->Fill(f_3, wightf3*weight);
                f4->Fill(f_4, wightf4*weight);
                f7->Fill(f_7, wightf7*weight);

                m6s->Fill(M_6s, wight1m6s*weight);
                m6s->Fill(-M_6s,wight2m6s*weight);
                m6c->Fill(M_6c, wight1m6c*weight);
                m6c->Fill(-M_6c,wight2m6c*weight);
                f5->Fill(f_5, wight1f5*weight);
                f5->Fill(-f_5,wight2f5*weight);
                f8->Fill(f_8, wight1f8*weight);
                f8->Fill(-f_8,wight2f8*weight);
                f9->Fill(f_9, wight1f9*weight);
                f9->Fill(-f_9,wight2f9*weight);
            }
        }
    }

    for (int i =0;i<15;i++){
		string Paratype = paralist[i];
		CalValue(Paratype,f1s,m6s,m6c,f3,f4,f5,f7,f8,f9);
	}
}

int main(int argc, char **argv) {

	string TagType = "";
	string SampleType = "";
	int Parity;

	// ########################
	// # Give format of input #
	// ########################
	if (argc<2){
		cout << "\n[Moment::main]\tParameters missing\t" << endl;
		cout << "\n[Moment::main]\tPlease input correct format like following\t" << endl;
		cout << "\n./Moment\tSampletype\tTagtype\tParity\t" << endl;
		cout << "\nTagtype:mis,good,mix" << endl;
		cout << "\nSampletype:gen,reco\t" << endl;
		cout << "\nParity:0,1\t" << endl;
		return EXIT_FAILURE;
	}


	// ####################
	// # select Data & MC #
	// ####################
	SampleType = argv[1];
	if ((SampleType != "gen") && (SampleType != "reco")) {
		cout << "\n[Moment::main]\tIncorrect Sample Type\t" << SampleType << endl;
		return EXIT_FAILURE;
	}
	cout << "\n[Moment::main]\tSample Type = \t" << SampleType << endl;
	if (SampleType == "reco")
	{
		// ##################
		// # Check Tag Type #
		// ##################
		TagType = argv[2];
		Parity = atoi(argv[3]);
		if ((TagType != "good") && (TagType != "mis") && (TagType !="mix")) {
			cout << "\n[Moment::main]\tIncorrect Events Tag Type\t" << TagType << endl;
			return EXIT_FAILURE;
		}
		if ((Parity!=0) && (Parity!=1)){
			cout << "\n[Moment::main]\tIncorrect Events Parity\t" << Parity << endl;
		}
		cout << "\n[Moment::main]\tEvents Tag Type = \t" << TagType << endl;
		cout << "\n[Moment::main]\tEvents Parity = \t" << Parity << endl;
	}

	// ###########################
	// # calculate final results #
	// ###########################
	for (q2Bin = 0; q2Bin < 8; q2Bin++) {
		if (q2Bin == 4||q2Bin == 6){
			continue;
		}
		quzhi();
		if (SampleType == "gen") 
			GenCalValue(q2Bin, q2Min, q2Max);
		else if (SampleType == "reco") 
			ReCalValue(q2Bin, q2Min, q2Max, TagType, Parity);
	}
	return 0;
}















