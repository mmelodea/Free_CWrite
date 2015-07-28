/*---------------------------------------------------
 * By Miqueias Melo
 * For use in other function do:
 * 1. Create a file.h with your function;
 * 2. Include it in this code (#include <file.h>);
 * 3. Modify the code calls function for your function name;
*/

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <TMath.h>
#include <TGraph.h>
//#include <my_fun.h>

double func(double alpha, double beta, double x_min, double x_max, double x)
{
  double a= (1+alpha*x+beta*x*x);
  double b= ((x_max-x_min)+0.5*alpha*(x_max*x_max-x_min*x_min)+(beta/3.)*(pow(x_max,3)-pow(x_min,3)));
  return a/b;
}


bi_ml()
{
      time_t date;
      struct tm* timeinfo;
            
  int ref,nl=0;
  const int n;
  double alpha,beta,x_min,x_max;
  
  cout << "Informe: alpha, beta, x_min, x_max, n_amostras"<<endl; /*, n_testes" << endl;*/
  cin >> alpha;
  cin >> beta;
  cin >> x_min;
  cin >> x_max;
  cin >> n;
  //cin >> ref;
  
  TH1D* h_teste = new TH1D("","Referencia",50,-1,1);
  h_teste->GetXaxis()->SetTitle("x"); h_teste->GetYaxis()->SetTitle("f(#alpha,#beta;x)");
  TH1D* h = new TH1D("","Meu",50,-1,1);
  TH2D* h2D = new TH2D("","beta x alpha",40,0,1,40,0,1);
  h2D->GetXaxis()->SetTitle("#alpha"); h2D->GetYaxis()->SetTitle("#beta");
  TH1D* h_alpha = new TH1D("","",40,0,1); h_alpha->GetXaxis()->SetTitle("#alpha");
  TH1D* h_beta = new TH1D("","",40,0,1);  h_beta->GetXaxis()->SetTitle("#beta");
  TH2D* h2D_lim_68 = new TH2D("","",500,0,1,500,0,1); h2D_lim_68->SetMarkerColor(kPink);
  h2D_lim_68->GetXaxis()->SetTitle("#alpha"); h2D_lim_68->GetYaxis()->SetTitle("#beta");
  TH2D* h2D_lim_86 = new TH2D("","Confidence Regions",500,0,1,500,0,1); h2D_lim_86->SetMarkerColor(kBlue);
  h2D_lim_86->GetXaxis()->SetTitle("#alpha"); h2D_lim_86->GetYaxis()->SetTitle("#beta");
  
 double sum_alpha=0.,sum_beta=0.,alpham=0.,betam=0.;
 //do
 //{
  double f_max= (func(alpha,beta,x_min,x_max,x_min)>func(alpha,beta,x_min,x_max,x_max))? func(alpha,beta,x_min,x_max,x_min):func(alpha,beta,x_min,x_max,x_max);
  //cout<<"f_max: "<<f_max<<endl;
  TRandom* r1 = new TRandom(NULL);
  TRandom* r2 = new TRandom(NULL+1);
  int p=0,f=0; double xf[n];
  
  time(&date);
  timeinfo = localtime(&date);
  cout<<"\nInicio: "<<asctime(timeinfo);

  do
  {
    p++;
    double rf_max = f_max*(r1->Rndm());
    double x = x_min+(x_max-x_min)*(r2->Rndm());
    double func_value = func(alpha,beta,x_min,x_max,x);
    if(rf_max<func_value) {h_teste->Fill(x,func_value/n); h->Fill(x,1); xf[f]=x; f=f+1; /*cout<<"x,func: "<<x<<"  |"<<func_value<<endl;*/}
  }
  while(f<n);
  
  //Likelihood tests
  const int N= 500; 
  double lnL_max=-1.E9,alpha_max = 0.,beta_max = 0.,r= 1./N,lnL[N][N];
  for(int i=1;i<N-1;i++)
  {
    //cout<<i-1<<", ";
    double alphai= i*r;
    for(int j=1;j<N-1;j++)
    {
      double betai=j*r;
      double lh= 0., sum= 0.;
      for(int k=0;k<=f;k++)
      {
	sum = TMath::Log(func(alphai,betai,x_min,x_max,xf[k]));
        lh = lh + sum;
        //sum = func(alphai,betai,x_min,x_max,xf[k]);
	//lh = lh*sum;
      }
      //double ln = TMath::Log(lh);
      if(lh>lnL_max) {lnL_max=lh; alpha_max = alphai; beta_max = betai;}
      lnL[i-1][j-1] = lh;
      //cout<<"lh: "<<lh<<endl;
            
    }
    //cout<<"alpha: "<<alpha;
    //gSystem->Exec("clear");
  }
  h2D->Fill(alpha_max,beta_max,1);
    
  nl= nl+1;
  //gSystem->Exec("clear");
  //cout<<"nl: "<<nl<<endl;
  
  cout<<"alpha_max: "<<alpha_max<<" |beta_max: "<<beta_max<<" |lnL_max:"<<lnL_max<<endl;
  
  h_alpha->Fill(alpha_max,1);
  h_beta->Fill(beta_max,1);
  //cout<<"Ok, here final! No:"<<nl<<endl;
  
  // Confidence Regions plots
  //--------------------------------------------------------
  double MLEia68=alpha, MLEfa68=alpha, MLEib68=beta, MLEfb68=beta;
  double MLEia86=alpha, MLEfa86=alpha, MLEib86=beta, MLEfb86=beta;
  double xv[N*N],yv[N*N],zv[N*N]; int T=0;
  for(int i=1;i<N-1;i++)
  {
    double alphai= i*r;
    for(int j=1;j<N-1;j++)
    {
      double betai=j*r;
      if(fabs(lnL[i-1][j-1]-(lnL_max-1))<0.01)
      {
	if(alphai<MLEia68)MLEia68=alphai;	if(alphai>MLEfa68)MLEfa68=alphai;
	if(betai<MLEib68)MLEib68=betai;		if(betai>MLEfb68)MLEfb68=betai;
	h2D_lim_68->Fill(alphai,betai,1); xv[T]=alphai; yv[T]=betai; zv[T]=1.; T=T+1;
      }
      if(fabs(lnL[i-1][j-1]-(lnL_max-2))<0.01)
      {
	if(alphai<MLEia86)MLEia86=alphai;	if(alphai>MLEfa86)MLEfa86=alphai;
	if(betai<MLEib86)MLEib86=betai;		if(betai>MLEfb86)MLEfb86=betai;
	h2D_lim_86->Fill(alphai,betai,1);
      }
    }
  }
  double ia,fa,ib,fb;
  ia=alpha_max-0.02; fa=alpha_max+0.02; ib=beta_max-0.02; fb=beta_max+0.02;
  double iap,fap,ibp,fbp;
  iap=alpha-0.02; fap=alpha+0.02; ibp=beta-0.02; fbp=beta+0.02;
  
  TLine* l1; TLine* l2; TLine* l3; TLine* l4;
  l1 = new TLine(ia,beta_max,fa,beta_max); l1->SetLineColor(kBlack);
  l2 = new TLine(alpha_max,ib,alpha_max,fb); l2->SetLineColor(kBlack);
  l3 = new TLine(iap,beta,fap,beta); l3->SetLineColor(kRed);
  l4 = new TLine(alpha,ibp,alpha,fbp); l4->SetLineColor(kRed);
  
  TLatex t;	t.SetTextSize(0.05);
  TString alphat,betat,lai68,lai86,laf68,laf86,lbi68,lbi86,lbf68,lbf86,alphatref,betatref;
  std::stringstream sac;
  sac << alpha_max; alphat=sac.str(); sac.str("");
  sac << alpha_max-MLEia68; lai68=sac.str(); sac.str("");	sac << MLEfa68-alpha_max; laf68=sac.str(); sac.str("");
  sac << alpha_max-MLEia86; lai86=sac.str(); sac.str("");	sac << MLEfa86-alpha_max; laf86=sac.str(); sac.str("");
  sac << beta_max; betat=sac.str(); sac.str("");
  sac << beta_max-MLEib68; lbi68=sac.str(); sac.str("");	sac << MLEfb68-beta_max; lbf68=sac.str(); sac.str("");
  sac << beta_max-MLEib86; lbi86=sac.str(); sac.str("");	sac << MLEfb86-beta_max; lbf86=sac.str(); sac.str("");
  
  sac << alpha; alphatref=sac.str(); sac.str("");	sac << beta; betatref=sac.str(); sac.str("");
  //--------------------------------------------------------
 //}
 //while(nl<ref);  
  
  TCanvas* c = new TCanvas();
  c->Divide(2,1);
  c->cd(1); /*h_teste->Fit("gaus","LL","",-1,1);*/ h_teste->Draw();  //Plot the function values and fit
  t.SetTextColor(kBlack); t.DrawLatex(-0.9,0.005,"#alpha = "+alphatref); t.DrawLatex(-0.9,0.015,"#beta = "+betatref);
  //c->cd(2); h->Draw();						 //Plot how many random max function values is below to random funtion value
  //c->cd(3); h2D->Draw("Colz");	//(Cont4 put color on base)	 //Plot 2D distribution for beta and alpha values to each ML
  //c->cd(4); h_alpha->Draw();					 //Plot distribution for alpha to each ML run
  //c->cd(5); h_beta->Draw();					 // "       "            beta      "       "
  c->cd(2); h2D_lim_68->Draw("P"); h2D_lim_86->Draw("same");	 //Plot 2D distribution for beta and alpha values for each ML value
  l1->Draw(); l2->Draw(); l3->Draw(); l4->Draw();
  t.SetTextColor(kRed);	t.DrawLatex(0.02,0.92,"63.2% CL");
  t.DrawLatex(0.02,0.87,"#alpha = "+alphat+"_{- "+lai68+"}^{+"+laf68+"}");
  t.DrawLatex(0.02,0.8,"#beta = "+betat+"_{- "+lbi68+"}^{+"+lbf68+"}");
  t.SetTextColor(kBlue); t.DrawLatex(0.5,0.92,"86.5% CL");
  t.DrawLatex(0.5,0.87,"#alpha = "+alphat+"_{- "+lai86+"}^{+"+laf86+"}");
  t.DrawLatex(0.5,0.8,"#beta = "+betat+"_{- "+lbi86+"}^{+"+lbf86+"}");
  
  
  /*TFile* bi_ml = new TFile("bi_ml.root","recreate");
  h_teste->Write();
  h->Write();
  h_alpha->Write();
  h_beta->Write();
  h2D->Write();
  h2D_lim_68->Write();
  h2D_lim_86->Write();*/
  
  //--- Informs MLE ---
  cout<<"\n********************** MLE ***************************"<<endl;
  cout<<"******************************************************"<<endl;
  cout<<" alpha: "<<alpha_max<<" -"<<alpha_max-MLEia68<<" +"<<MLEfa68-alpha_max<<"  (63.2% CL)"<<endl;
  cout<<"  beta: "<<beta_max<<" -"<<beta_max-MLEib68<<" +"<<MLEfb68-beta_max<<"  (63.2% CL)"<<endl;
  cout<<" alpha: "<<alpha_max<<" -"<<alpha_max-MLEia86<<" +"<<MLEfa86-alpha_max<<"  (86.5% CL)"<<endl;
  cout<<"  beta: "<<beta_max<<" -"<<beta_max-MLEib86<<" +"<<MLEfb86-beta_max<<"  (86.5% CL)"<<endl;
  cout<<"******************************************************"<<endl;
  
  time(&date);
  timeinfo = localtime(&date);
  cout<<"Termino: "<<asctime(timeinfo)<<endl;
  
  TCanvas* c2 = new TCanvas();
  TGraph* gteste = new TGraph(T,xv,yv);
  gteste->Draw("AP");

  return 0;
}//General End
