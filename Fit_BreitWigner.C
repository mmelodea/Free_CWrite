//Funcao Breit Wigner
Double_t bw(Double_t m4l, Double_t mZ, Double_t Lg)
{
  //WorkBook: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookHowToFit
  //Phys. Rev. D86 (2012) pag. 479, propagador na Eq.(3)
  return pow(Lg*mZ,2)/(pow((m4l*m4l-mZ*mZ),2) + pow(m4l,4)*pow(Lg/mZ,2));
}

#include "TF1.h"
#include "TH1.h"
#include "TMath.h"
#include "TFile.h"
#include "TLatex.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>

void Fit_BreitWigner()
{ 

 float xi,xf,L,M;
 int CL,dc,sample;
 
   cout<<"Data(1), MC(2) ou Merge(3)? "; cin >> sample;
   
if(sample == 1 || sample == 2){
 cout<< "Obter massa(1), largura(2) ou ambos(3)? "; cin >> dc;
 if(dc == 1){ cout<< "Informe a largura a ser usada: "; cin >> L; }
 if(dc == 2){ cout<< "Informe a massa a ser usada: "; cin >> M;   }
 
 TH1D* Z1 = new TH1D("","Z1 Reconstruido",40,40,120);
 if(dc == 3){ 
   TH2D* Lh = new TH2D("","Likelihood",80000,40,120,3000,0,3);
   Lh->GetXaxis()->SetTitle("M_Z (GeV)");
   Lh->GetYaxis()->SetTitle("#Gamma (GeV)");
            }
 
 if(sample == 1) TFile* data = TFile::Open("$HOME/Dropbox/my_new_analysis/minhas_rootplas/DADOS/Data.root");
 if(sample == 2) TFile* data = TFile::Open("$HOME/Dropbox/my_new_analysis/minhas_rootplas/MC/Data.root");
 TH1D* Z_1=(TH1D*)data->Get("h_Z1"); //Separa o histograma para Fit
 Int_t   division = Z_1->GetNbinsX();
 Float_t massMIN = Z_1->GetBinLowEdge(1);
 Float_t massMAX = Z_1->GetBinLowEdge(division+1);
 Float_t BIN_SIZE = Z_1->GetBinWidth(1);
 Z_1->SetLineColor(kRed);

 std::vector<float> input;
 float c;
 FILE *f = fopen("Z_masses_all.txt","r");
 while(fscanf(f,"%f",&c) !=EOF){ input.push_back(c); Z1->Fill(c); };
 fclose(f);
/*//----------------------------------------
 double mZ=91.188, Lg=2.495, lh=0., sum=0.;
 for(int z=0;z<input.size();z++){
   sum = TMath::Log(bw(input[z],mZ,Lg));
   lh = lh + sum;
                                }
 double Lh_Z_SM = lh; //Likelihood Z SM hypotesis*/
//========== Likelihoods test ============
//Precisao das verificacoes
 const int T=1000;
 double d[T];
 for(int s=1;s<=T;s++){ d[s-1]=s/double(T); }
 double lh_max=-1.E5, mZ_max=0., Lg_max=0., mZ=0., Lg=0., lh=0., sum=0.;
 std::vector<double> lhv, mZv, Lgv;
 int st1, fh1, st2, fh2;
 
 if(dc == 1){ Lg=L; st1=40; fh1=120; }
 if(dc == 2){ mZ=M; st1=0; fh1=3; }
   for(int j=st1;j<fh1;j++){
         for(int l=0;l<T;l++){
          if(dc == 1){ mZ = j + d[l]; }
          if(dc == 2){ Lg = j + d[l]; }
          lh=0., sum=0.;
          for(int i=0;i<input.size();i++){
           sum = TMath::Log(bw(input[i],mZ,Lg));
           lh = lh + sum;
                                         }
//         cout<<"mZ: "<<mZ<<"  |Lg: "<<Lg<<"  |lh: "<<lh<<endl;
          if(dc == 1) mZv.push_back(mZ);
          if(dc == 2) Lgv.push_back(Lg);
	  lhv.push_back(lh);
          //lhv.push_back(2*(Lh_Z_SM-lh));
          if(lh > lh_max){ lh_max=lh; mZ_max=mZ; Lg_max=Lg; }
                             }
                           }
  
 const int N = lhv.size();
 double lhV[N],mZV[N],LgV[N];
 if(dc==1){ for(int o=0;o<mZv.size();o++){ lhV[o]=lh_max-lhv[o]; mZV[o]=mZv[o]; } TGraph* m = new TGraph(N,mZV,lhV); }
 if(dc==2){ for(int o=0;o<Lgv.size();o++){ lhV[o]=lh_max-lhv[o]; LgV[o]=Lgv[o]; } TGraph* lg = new TGraph(N,LgV,lhV); }
   
 if(dc == 3){
 int st_mZ=40, fh_mZ=120, st_Lg=0, fh_Lg=3;
 for(int i=st_Lg;i<fh_Lg;i++){
        for(int lL=0;lL<T;lL++){ 
          Lg = i + d[lL];
          for(int j=st_mZ;j<fh_mZ;j++){
                for(int lm=0;lm<T;lm++){
                  mZ = j + d[lm];
	          double lh=0., sum=0.;
                  for(int k=0;k<input.size();k++){
                    sum = TMath::Log(bw(input[i],mZ,Lg));
                    lh = lh + sum;
                                                 }
                  Lh->Fill(mZ,Lg,lh);
                  if(lh > lh_max){ lh_max = lh; mZ_max = mZ; Lg_max = Lg; }
                                       }
	       	                      }
	    	               }
                             }
            }
// cout<<"Lh_Z_SM: "<<Lh_Z_SM<<endl;
 cout<<"mZ_max: "<<mZ_max<<" |Lg_max: "<<Lg_max<<" GeV"<<" |lh_c: "<<lh_max<<endl;
//================================================================
 TString title = (sample == 1)? "Data":"MC";  
 TCanvas *c1 = new TCanvas("c1",title,10,10,1000,500);
 c1->Divide(2,1);
 c1->cd(1);
 Z_1->Draw(); Z1->Draw("same");
 c1->cd(2);
 TLine* l1 = new TLine(90.5,0.5,93,0.5); l1->SetLineStyle(7);
 TLine* l2 = new TLine(90.5,2,93,2); l2->SetLineStyle(7);
 if(dc == 1){ m->GetYaxis()->SetTitle("2#Delta ln#mathcal{L}"); m->Draw("AL"); l1->Draw(); l2->Draw(); }
 if(dc == 2){ lg->GetYaxis()->SetTitle("2#Delta ln#mathcal{L}"); lg->Draw("AL"); l1->Draw(); l2->Draw(); }
 if(dc == 3) Lh->Draw("Colz");

 TFile ach("all.root","recreate");
 Z_1->Write();
 Z1->Write();
 if(dc == 1) m->Write();
 if(dc == 2) lg->Write();
 if(dc == 3) Lh->Write();
 ach.Close();
                              }

 if(sample == 3){
 cout<< "Mesclar os Likelihoods?(0/1) "; cin >> dc;
 if(dc == 1) continue;
 if(dc == 0){
   TMultiGraph* mg = new TMultiGraph();
   std::vector<string> channel, lc;
   channel.push_back("4e"); channel.push_back("4mu"); channel.push_back("all");
   lc.push_back("kRed"); lc.push_back("kBlue"); lc.push_back("kGreen");
   for(int i=0;i<channel.size();i++){
   TString path = "$HOME/Dropbox/my_new_analysis/"+channel[i]+".root";
   Color_t cor = lc[i];
   TFile* data = TFile::Open(path);
   TGraph* Z =(TGraph*)data->Get("Graph");
   Z->SetLineColor(cor);
   mg->Add(Z,"l");
                                    }
 mg->Draw("al");}
               }

}
