
double func(double x)
{
  Float_t l=4.;
  return (l/2)*sin(x);
}

void buffon_a()
{
gROOT->Reset();
c1=new TCanvas("c1","c1",10,10,600,600);


Int_t n=0;
const Float_t N=5000;
Float_t a=3.14,l=4., d=5.;
Double_t x[N], y[N],hy[N];
Double_t m,soma=0,pimed;

hist=new TH1D ("N_{exp}=2000","Experimento de Buffon",50,3,3.3);
hist->SetFillColor(4);
hist->SetFillStyle(3003);
hist->GetXaxis()->SetTitle("Valor de #pi");
hist->GetXaxis()->CenterTitle();
hist->GetYaxis()->SetTitle("Eventos");
hist->GetYaxis()->CenterTitle();

TRandom *r=new TRandom(NULL); //gerador de números aleatórios no intervalo [0,1] (Rndm ativa o gerador)

 do
  {n=n+1;
m=0;
//cout<<"n="<<n<<"\n";
for(Int_t i=1; i<=N; i++)
{
y[i]=(d/2)*(r->Rndm(NULL));
x[i]=a*(r->Rndm(NULL));
hy[i]=0;
if (y[i]<= func(x[i]))
  {m=m+1;}
} 
//cout<<"m="<<m<<"\n";
Float_t pi=(2*l/d)*(N/m);
//cout<<"pi="<<pi<<"\n";
hist->Fill(pi);
soma=soma+pi;
//cout<<"soma="<<soma<<"\n";
hist->Draw();
c1->Update();  
}
while(n<500);
pimed=soma/n;

cout<< "###################################\n"<<"O VALOR MÉDIO DE PI É: "<<pimed<<"\n###################################"<<endl;

/*
Estas linhas salvam o histograma gerado em um arquivo ".root" que pode ser lido (gerando igualmente o que se produz no programa, incluindo legendas, títulos e eixos) posteriormente. Os comandos de leitura são:

TFile f("arquivo.root");
TH1(F,D,...) *nome_do_histograma= (TH1(F,D,...)*)f.Get("nome_do_histograma");
nome_do_histograma->Draw(); 

O nome do histograma pode ser visto pelo comando: f.GetListOfKeys()->Print();
*/

TFile f("histbuffon.root","new");
hist->Write();
f.Close();

}
