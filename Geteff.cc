#include <TFile.h>
#include <RooWorkspace.h>
#include <RooDataSet.h>
#include <RooArgSet.h>
#include <RooArgList.h>
#include <RooRealVar.h>
#include <TH1.h>
#include <TMath.h>

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
    // q2-bin format: [0-8] for one bin
    // flags to mark which q2 bins should be filled

    cout << "year 201" << year << " q2Bin" << q2Bin << " binnum" << binnum << "parity" << parity << endl;
    string datasetString[5] = {"data_genDen","data_genNum","data_den","data_ctRECO","data_wtRECO"};
    string paraString = Form((parity==0?"_ev_b%i":"_od_b%i"),q2Bin);
    for (int i=0;i<5;i++){
        datasetString[i] = datasetString[i] + paraString;
    }

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

    RooDataSet* totdata[5];
    Int_t Entries[5];
    for (int i=0;i<5;i++){
        totdata[i] = (RooDataSet*)wsp->data(datasetString[i].c_str());
        if ( !totdata[i] || totdata[i]->IsZombie() ) {
            cout<<"Dataset "<<datasetString[i]<<" not found in file: "<<filename<<endl;
            return;
        }
        Entries[i] = totdata[i]->numEntries();
    }

    string filestring = Form(parity==0 ? "1Deffb%iy%i_ev.root" : "1Deffb%iy%i_od.root",q2Bin,year);
    cout << "efficiency will be saved in " << filestring << endl;
    TFile *fout = new TFile(filestring.c_str(),"recreate");
    TH1D *f[9][9];
    string name[9] = {"f1s","m6s","m6c","f3","f4","f5","f7","f8","f9"};
    string kind[9] = {"genDen","genNum","den","CTreco","WTreco","CTeff","WTeff","CTweight","WTweight"};
    double cos_theta_k,cos_theta_l,phi_kst_mumu,weight;
    double f_1s,M_6s,M_6c,f_3,f_4,f_5,f_7,f_8,f_9;

    vector<vector<double>> fi;
    for (int i=0;i<9;i++){
        fi.push_back(vector<double>());
    }
    cout<<"Starting adaptive-binning histogram setting using total genNum sample..."<<endl;
    //need to use all the genNum sample
    RooDataSet *totdata_genNum = totdata[1];
    RooWorkspace* wspextra = (RooWorkspace*)fin->Get(Form("ws_b%ip%i",q2Bin,1-parity));
    if ( !wsp || wspextra->IsZombie() ) {
        cout<<"Workspace not found in file: "<<filename<<endl;
        return;
    }
    string extradatastring = "data_genNum";
    extradatastring = extradatastring+Form((parity==0?"_od_b%i":"_ev_b%i"),q2Bin);
    RooDataSet *extradata = (RooDataSet*)wspextra->data(extradatastring.c_str());
    if ( !extradata || extradata->IsZombie() ) {
        cout<<"Dataset "<< extradatastring <<" not found in file: "<<filename<<endl;
        return;
    }
    totdata_genNum->append(*extradata);
    int totentries = totdata_genNum->numEntries();
    for (int i = 0; i<totentries; i++){
        cos_theta_k = totdata_genNum->get(i)->getRealValue("ctK");
        cos_theta_l = totdata_genNum->get(i)->getRealValue("ctL");
        phi_kst_mumu = totdata_genNum->get(i)->getRealValue("phi");
        f_1s = 1 - cos_theta_k * cos_theta_k;
		M_6s = (1 - cos_theta_k * cos_theta_k) * cos_theta_l;
		M_6c = cos_theta_k * cos_theta_k * cos_theta_l;
		f_3 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Cos(2 * phi_kst_mumu);
		f_4 = 4 * cos_theta_k * cos_theta_l * TMath::Cos(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
		f_5 = 2 * cos_theta_k * TMath::Cos(phi_kst_mumu)* TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
		f_7 = 2 * cos_theta_k * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
		f_8 = 4 * cos_theta_k * cos_theta_l * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
		f_9 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Sin(2 * phi_kst_mumu);
        fi[0].push_back(f_1s);
        fi[1].push_back(M_6s);
        fi[2].push_back(M_6c);
        fi[3].push_back(f_3);
        fi[4].push_back(f_4);
        fi[5].push_back(f_5);
        fi[6].push_back(f_7);
        fi[7].push_back(f_8);
        fi[8].push_back(f_9);
    }

    for (int i=0;i<9;i++){
        sort(fi[i].begin(),fi[i].end());
    }


    int countsPerBin=totentries/binnum;
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

    for (int i=0;i<9;i++){
        for (int k=0;k<9;k++){
            f[i][k]->SetName(Form("%s%s",name[i].c_str(),kind[k].c_str()));
        }
    }

    // Prepare GEN-level datasets
    cout<<"Starting datasets filling..."<<endl;
    for (int k = 0;k < 5; k++){
        cout<<"Starting datasets " << kind[k] << " filling" << endl;
        for (int i = 0;i < Entries[k]; i++){
            cos_theta_k = totdata[k]->get(i)->getRealValue("ctK");
            cos_theta_l = totdata[k]->get(i)->getRealValue("ctL");
            phi_kst_mumu = totdata[k]->get(i)->getRealValue("phi");
            weight = totdata[k]->weight();
            f_1s = 1 - cos_theta_k * cos_theta_k;
            M_6s = (1 - cos_theta_k * cos_theta_k) * cos_theta_l;
            M_6c = cos_theta_k * cos_theta_k * cos_theta_l;
            f_3 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Cos(2 * phi_kst_mumu);
            f_4 = 4 * cos_theta_k * cos_theta_l * TMath::Cos(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
            f_5 = 2 * cos_theta_k * TMath::Cos(phi_kst_mumu)* TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
            f_7 = 2 * cos_theta_k * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
            f_8 = 4 * cos_theta_k * cos_theta_l * TMath::Sin(phi_kst_mumu) * TMath::Sqrt((1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l));
            f_9 = (1 - cos_theta_k * cos_theta_k) * (1 - cos_theta_l * cos_theta_l) * TMath::Sin(2 * phi_kst_mumu);
            if (k==4){
                M_6s = -M_6s;
                M_6c = -M_6c;
                f_5 = -f_5;
                f_8 = -f_8;
                f_9 = -f_9;
            }
            f[0][k]->Fill(f_1s,weight);
            f[1][k]->Fill(M_6s,weight);
            f[2][k]->Fill(M_6c,weight);
            f[3][k]->Fill(f_3,weight);
            f[4][k]->Fill(f_4,weight);
            f[5][k]->Fill(f_5,weight);
            f[6][k]->Fill(f_7,weight);
            f[7][k]->Fill(f_8,weight);
            f[8][k]->Fill(f_9,weight);
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




