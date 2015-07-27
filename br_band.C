
void br_band(){
  
  TGraph *cl_obs = new TGraph();
  cl_obs->SetLineColor(kBlack);
  cl_obs->SetMarkerStyle(21);
  cl_obs->SetLineWidth(2);
    
  TGraph *cl_exp = new TGraph();
  cl_exp->SetLineColor(kBlue);
  cl_exp->SetLineStyle(2);
  cl_exp->SetLineWidth(2);

//Exclusion zone 68%
  TGraph *cl_exp_1sm = new TGraph(); 		
  cl_exp_1sm->SetLineColor(kGreen);	
  cl_exp_1sm->SetFillColor(kGreen);	
  cl_exp_1sm->SetLineWidth(2);	
  
  TGraph *cl_exp_1sp = new TGraph();
  cl_exp_1sp->SetLineColor(kGreen);
  cl_exp_1sp->SetFillColor(kGreen);	
  cl_exp_1sp->SetLineWidth(2);
  
//Exclusion zone 97.5%  
  TGraph *cl_exp_2sm = new TGraph(); 		
  cl_exp_2sm->SetLineColor(kYellow);	
  cl_exp_2sm->SetFillColor(kYellow);	
  cl_exp_2sm->SetLineWidth(2);	
  
  TGraph *cl_exp_2sp = new TGraph();
  cl_exp_2sp->SetLineColor(kYellow);
  cl_exp_2sp->SetFillColor(kYellow);	
  cl_exp_2sp->SetLineWidth(2);

  Double_t limit;
/*  
  TFile *f170 = TFile::Open("170/higgsCombineTest.Asymptotic.mH170.root");
  TTree *t170 = (TTree*)f170->Get("limit"); t170->SetDirectory(0); 
  t170->SetBranchAddress("limit",&limit);
  t170->GetEntry(0); cl_exp_2sm->SetPoint(0,170,limit);
  t170->GetEntry(1); cl_exp_1sm->SetPoint(0,170,limit);
  t170->GetEntry(2); cl_exp->SetPoint(0,170,limit);
  t170->GetEntry(3); cl_exp_1sp->SetPoint(0,170,limit);
  t170->GetEntry(4); cl_exp_2sp->SetPoint(0,170,limit);
  t170->GetEntry(5); cl_obs->SetPoint(0,170,limit);
*/
  TFile *f180 = TFile::Open("180/higgsCombineTest.Asymptotic.mH180.root");
  TTree *t180 = (TTree*)f180->Get("limit"); t180->SetDirectory(0); 
  t180->SetBranchAddress("limit",&limit);
  t180->GetEntry(0); cl_exp_2sm->SetPoint(0,180,limit);
  t180->GetEntry(1); cl_exp_1sm->SetPoint(0,180,limit);
  t180->GetEntry(2); cl_exp->SetPoint(0,180,limit);
  t180->GetEntry(3); cl_exp_1sp->SetPoint(0,180,limit);
  t180->GetEntry(4); cl_exp_2sp->SetPoint(0,180,limit);
  t180->GetEntry(5); cl_obs->SetPoint(0,180,limit);
  
  TFile *f190 = TFile::Open("190/higgsCombineTest.Asymptotic.mH190.root");
  TTree *t190 = (TTree*)f190->Get("limit"); t190->SetDirectory(0); 
  t190->SetBranchAddress("limit",&limit);
  t190->GetEntry(0); cl_exp_2sm->SetPoint(1,190,limit);
  t190->GetEntry(1); cl_exp_1sm->SetPoint(1,190,limit);
  t190->GetEntry(2); cl_exp->SetPoint(1,190,limit);
  t190->GetEntry(3); cl_exp_1sp->SetPoint(1,190,limit);
  t190->GetEntry(4); cl_exp_2sp->SetPoint(1,190,limit);
  t190->GetEntry(5); cl_obs->SetPoint(1,190,limit);

  TFile *f200 = TFile::Open("200/higgsCombineTest.Asymptotic.mH200.root");
  TTree *t200 = (TTree*)f200->Get("limit"); t200->SetDirectory(0); 
  t200->SetBranchAddress("limit",&limit);
  t200->GetEntry(0); cl_exp_2sm->SetPoint(2,200,limit);
  t200->GetEntry(1); cl_exp_1sm->SetPoint(2,200,limit);
  t200->GetEntry(2); cl_exp->SetPoint(2,200,limit);
  t200->GetEntry(3); cl_exp_1sp->SetPoint(2,200,limit);
  t200->GetEntry(4); cl_exp_2sp->SetPoint(2,200,limit);
  t200->GetEntry(5); cl_obs->SetPoint(2,200,limit);

  TFile *f250 = TFile::Open("250/higgsCombineTest.Asymptotic.mH250.root");
  TTree *t250 = (TTree*)f250->Get("limit"); t250->SetDirectory(0); 
  t250->SetBranchAddress("limit",&limit);
  t250->GetEntry(0); cl_exp_2sm->SetPoint(3,250,limit);
  t250->GetEntry(1); cl_exp_1sm->SetPoint(3,250,limit);
  t250->GetEntry(2); cl_exp->SetPoint(3,250,limit);
  t250->GetEntry(3); cl_exp_1sp->SetPoint(3,250,limit);
  t250->GetEntry(4); cl_exp_2sp->SetPoint(3,250,limit);
  t250->GetEntry(5); cl_obs->SetPoint(3,250,limit);

  TFile *f300 = TFile::Open("300/higgsCombineTest.Asymptotic.mH300.root");
  TTree *t300 = (TTree*)f300->Get("limit"); t300->SetDirectory(0); 
  t300->SetBranchAddress("limit",&limit);
  t300->GetEntry(0); cl_exp_2sm->SetPoint(4,300,limit);
  t300->GetEntry(1); cl_exp_1sm->SetPoint(4,300,limit);
  t300->GetEntry(2); cl_exp->SetPoint(4,300,limit);
  t300->GetEntry(3); cl_exp_1sp->SetPoint(4,300,limit);
  t300->GetEntry(4); cl_exp_2sp->SetPoint(4,300,limit);
  t300->GetEntry(5); cl_obs->SetPoint(4,300,limit);

  TFile *f350 = TFile::Open("350/higgsCombineTest.Asymptotic.mH350.root");
  TTree *t350 = (TTree*)f350->Get("limit"); t350->SetDirectory(0); 
  t350->SetBranchAddress("limit",&limit);
  t350->GetEntry(0); cl_exp_2sm->SetPoint(5,350,limit);
  t350->GetEntry(1); cl_exp_1sm->SetPoint(5,350,limit);
  t350->GetEntry(2); cl_exp->SetPoint(5,350,limit);
  t350->GetEntry(3); cl_exp_1sp->SetPoint(5,350,limit);
  t350->GetEntry(4); cl_exp_2sp->SetPoint(5,350,limit);
  t350->GetEntry(5); cl_obs->SetPoint(5,350,limit);

  TFile *f400 = TFile::Open("400/higgsCombineTest.Asymptotic.mH400.root");
  TTree *t400 = (TTree*)f400->Get("limit"); t400->SetDirectory(0); 
  t400->SetBranchAddress("limit",&limit);
  t400->GetEntry(0); cl_exp_2sm->SetPoint(6,400,limit);
  t400->GetEntry(1); cl_exp_1sm->SetPoint(6,400,limit);
  t400->GetEntry(2); cl_exp->SetPoint(6,400,limit);
  t400->GetEntry(3); cl_exp_1sp->SetPoint(6,400,limit);
  t400->GetEntry(4); cl_exp_2sp->SetPoint(6,400,limit);
  t400->GetEntry(5); cl_obs->SetPoint(6,400,limit);

  TFile *f450 = TFile::Open("450/higgsCombineTest.Asymptotic.mH450.root");
  TTree *t450 = (TTree*)f450->Get("limit"); t450->SetDirectory(0); 
  t450->SetBranchAddress("limit",&limit);
  t450->GetEntry(0); cl_exp_2sm->SetPoint(7,450,limit);
  t450->GetEntry(1); cl_exp_1sm->SetPoint(7,450,limit);
  t450->GetEntry(2); cl_exp->SetPoint(7,450,limit);
  t450->GetEntry(3); cl_exp_1sp->SetPoint(7,450,limit);
  t450->GetEntry(4); cl_exp_2sp->SetPoint(7,450,limit);
  t450->GetEntry(5); cl_obs->SetPoint(7,450,limit);

  TFile *f500 = TFile::Open("500/higgsCombineTest.Asymptotic.mH500.root");
  TTree *t500 = (TTree*)f500->Get("limit"); t500->SetDirectory(0); 
  t500->SetBranchAddress("limit",&limit);
  t500->GetEntry(0); cl_exp_2sm->SetPoint(8,500,limit);
  t500->GetEntry(1); cl_exp_1sm->SetPoint(8,500,limit);
  t500->GetEntry(2); cl_exp->SetPoint(8,500,limit);
  t500->GetEntry(3); cl_exp_1sp->SetPoint(8,500,limit);
  t500->GetEntry(4); cl_exp_2sp->SetPoint(8,500,limit);
  t500->GetEntry(5); cl_obs->SetPoint(8,500,limit);

  TFile *f550 = TFile::Open("550/higgsCombineTest.Asymptotic.mH550.root");
  TTree *t550 = (TTree*)f550->Get("limit"); t550->SetDirectory(0); 
  t550->SetBranchAddress("limit",&limit);
  t550->GetEntry(0); cl_exp_2sm->SetPoint(9,550,limit);
  t550->GetEntry(1); cl_exp_1sm->SetPoint(9,550,limit);
  t550->GetEntry(2); cl_exp->SetPoint(9,550,limit);
  t550->GetEntry(3); cl_exp_1sp->SetPoint(9,550,limit);
  t550->GetEntry(4); cl_exp_2sp->SetPoint(9,550,limit);
  t550->GetEntry(5); cl_obs->SetPoint(9,550,limit);

  TFile *f600 = TFile::Open("600/higgsCombineTest.Asymptotic.mH600.root");
  TTree *t600 = (TTree*)f600->Get("limit"); t600->SetDirectory(0); 
  t600->SetBranchAddress("limit",&limit);
  t600->GetEntry(0); cl_exp_2sm->SetPoint(10,600,limit);
  t600->GetEntry(1); cl_exp_1sm->SetPoint(10,600,limit);
  t600->GetEntry(2); cl_exp->SetPoint(10,600,limit);
  t600->GetEntry(3); cl_exp_1sp->SetPoint(10,600,limit);
  t600->GetEntry(4); cl_exp_2sp->SetPoint(10,600,limit);
  t600->GetEntry(5); cl_obs->SetPoint(10,600,limit);

  
  TLegend *leg = new TLegend(0.3,0.7,0.7,0.9);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(cl_obs,"Limite Observado 95%C.L.","l");
  leg->AddEntry(cl_exp,"Limite Esperado 95%C.L.","l");
  leg->AddEntry(cl_exp_1sm,"Limite Esperado #pm1#sigma","l");
  leg->AddEntry(cl_exp_2sm,"Limite Esperado #pm2#sigma","l");
  
  TLine *l1 = new TLine(150,1,640,1); l1->SetLineColor(kBlack); l1->SetLineStyle(2); l1->SetLineWidth(2); 
  TCanvas *c = new TCanvas();
  c->SetLogy();
  
  TF1 *exp_2sp = new TF1("exp_2sp","[0]*([1]*TMath::Landau(x,[2],[3])+[4]*TMath::Gaus(x,[5],[6])+[7]*TMath::Gaus(x,[8],[9])+[10]*TMath::Gaus(x,[11],[12]))",138,900);
  exp_2sp->FixParameter(0,1);
  exp_2sp->FixParameter(1,3.42182);
  exp_2sp->FixParameter(2,196.347);
  exp_2sp->FixParameter(3,5.3804);
  exp_2sp->FixParameter(4,1.2054);
  exp_2sp->FixParameter(5,228.32);
  exp_2sp->FixParameter(6,-77.1402);
  exp_2sp->FixParameter(7,26.7502);
  exp_2sp->FixParameter(8,400.86);
  exp_2sp->FixParameter(9,14.524);
  exp_2sp->FixParameter(10,2.66632);
  exp_2sp->FixParameter(11,793.723);
  exp_2sp->FixParameter(12,190.615);
  TH1D *hexp_2sp = new TH1D("hexp_2sp","",4700,150,620);
  hexp_2sp->SetLineColor(kYellow); hexp_2sp->SetFillColor(kYellow);
  for(int i=0; i<4700; i++){
    float value = exp_2sp->Eval(150+0.1*i);
    if((150+0.1*i)>=180 && (150+0.1*i)<=600) hexp_2sp->Fill(150+0.1*i+0.05,value);
  }
  exp_2sp->FixParameter(0,0.71);
  TH1D *hexp_1sp = new TH1D("hexp_1sp","",4700,150,620);
  hexp_1sp->SetLineColor(kGray); hexp_1sp->SetFillColor(kGreen);
  for(int i=0; i<4700; i++){
    float value = exp_2sp->Eval(150+0.1*i);
    if((150+0.1*i)>=180 && (150+0.1*i)<=600) hexp_1sp->Fill(150+0.1*i+0.05,value);
  }

  
  TF1 *exp_1sm = new TF1("exp_1sm","[0]*([1]*TMath::Landau(x,[2],[3])+[4]*TMath::Gaus(x,[5],[6])+[7]*TMath::Gaus(x,[8],[9])+[10]*TMath::Gaus(x,[11],[12]))",138,900);
  exp_1sm->FixParameter(0,0.299178);
  exp_1sm->FixParameter(1,3.26033);
  exp_1sm->FixParameter(2,197.492);
  exp_1sm->FixParameter(3,7.08859);
  exp_1sm->FixParameter(4,1.06731);
  exp_1sm->FixParameter(5,227.156);
  exp_1sm->FixParameter(6,-72.2686);
  exp_1sm->FixParameter(7,2.89731);
  exp_1sm->FixParameter(8,400.498);
  exp_1sm->FixParameter(9,17.0729);
  exp_1sm->FixParameter(10,12.6117);
  exp_1sm->FixParameter(11,1309.89);
  exp_1sm->FixParameter(12,344.295);
  TH1D *hexp_2sm = new TH1D("hexp_2sm","",4700,150,620);
  hexp_2sm->SetLineColor(kGray); hexp_2sm->SetFillColor(kGray);
  for(int i=0; i<4700; i++){
    float value = exp_1sm->Eval(150+0.1*i);
    if((150+0.1*i)>=180 && (150+0.1*i)<=600) hexp_2sm->Fill(150+0.1*i+0.05,value);
  }
  
  exp_1sm->FixParameter(0,0.4);
  exp_1sm->FixParameter(7,6.5);
  TH1D *hexp_1sm = new TH1D("hexp_1sm","",4700,150,620);
  hexp_1sm->SetLineColor(kGray); hexp_1sm->SetFillColor(kYellow);
  for(int i=0; i<4700; i++){
    float value = exp_1sm->Eval(150+0.1*i);
    if((150+0.1*i)>=180 && (150+0.1*i)<=600) hexp_1sm->Fill(150+0.1*i+0.05,value);
  }
  
  
  hexp_2sp->GetYaxis()->SetRangeUser(0.001,100);
  hexp_2sp->GetYaxis()->SetTitle("#sigma_{excluded}/#sigma_{SM}");
  hexp_2sp->GetXaxis()->SetTitle("M_{H}(GeV)");
  
  hexp_2sp->Draw();
  hexp_1sp->Draw("same");
  hexp_1sm->Draw("same");
  hexp_2sm->Draw("same");
  cl_exp->Draw("LP,same");
  cl_obs->Draw("LP,same");
  l1->Draw();
  leg->Draw();
  
}
