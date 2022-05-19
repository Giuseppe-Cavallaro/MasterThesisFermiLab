#define NFolders 1
#define NTrees 1

//Change this to change the used folder
#define NPath 0
#define SelectedTree 0



void TreePlot(){
std::string FolderNames[NFolders] = {"flashmatchAna"};

std::string xname = "Middlez";
std::string yname = "Middley";
std::string weight = "GammaStep";
std::string gamma = "GammaScint";

std::string FileName = "/media/giuseppe/01D5DC4DD5EAB920/lightmap/analisi_mu_14500.root";
//std::string FileName = "/media/giuseppe/01D5DC4DD5EAB920/lightmap/analisi_mu_11500.root";

std::string TreeNames[NTrees] = {"FlashMatchTree"};

std::string TitleDef = "#frac{SumPE}{EnergyDepositedTotal}:Abs(TrueX), met#grave{a} dei punti disponibili; " + FolderNames[NPath];
TFile *fp = new TFile(FileName.c_str(), "UPDATE");
  //Check if file is opened
  if (fp->IsZombie()) {
     std::cout << "Error opening file" << std::endl;
     exit(-1);
  }
  TTree *T;
  //std::cout<<"hi!"<<endl;

  fp->cd(FolderNames[NPath].c_str());
  gDirectory->GetObject(TreeNames[SelectedTree].c_str(), T);

  std::vector<float> * v1=0;
  std::vector<float> * v2=0;
  std::vector<float> * v3=0;
  std::vector<float> * v4=0;
  Double_t x=0,y=0,w=0,w2=0;
  
  T->SetBranchAddress(xname.c_str(),&v1);
  T->SetBranchAddress(yname.c_str(),&v2);
  T->SetBranchAddress(weight.c_str(),&v3);
  T->SetBranchAddress(gamma.c_str(),&v4);
  
   TH2* h2 = new TH2F("h2", "Y-Z LightMap", 35, 0, 1350, 35, -550, 550);
   h2->GetXaxis()->SetTitle("Z [cm]");
   h2->GetYaxis()->SetTitle("Y [cm]");
   
   
  Int_t Dimension = T->GetEntries();
  Int_t n = 0;
  
  
  
  for (Int_t i = 0; i<Dimension; i++) {
    T -> GetEntry(i);
    n = v1->size();
    for (Int_t z = 0; z<n-2; z++) {
    x=v1->at(z);
    y=v2->at(z);
    w=(v3->at(z))/0.05;
    w2=0;
    for(Int_t j = 0; j<60; j++){
    w2+=v4->at(j+60*z);
    }
    w2=w2/0.05;
    w=w/w2;
    h2->Fill(x,y,w);
    //h2->TH2::Interpolate(x,y,w);
    //myfile_2 <<x<<","<< y <<","<< w <<"\n"; //write to file
   
  }
  
 }
  //std::cout<<"hi!"<<std::endl;
  
  
 /* ofstream myfile;
  myfile.open ("i1-i2-binyz.txt");
  myfile << "coord1 coord2 bincontent\n";
  
  Double_t binvalue=0;
  for (Int_t i = 0; i<35; i++) {
  
  for (Int_t j = 0; j<35; j++) {
  binvalue=h2->GetBinContent(i,j);
  myfile <<i<<","<< j <<","<< binvalue <<"\n"; //write to file
  }
  
  }
  
  myfile.close();*/
  //ofstream myfile;
  //myfile.open ("y-z-bin.txt");
  //myfile << "coord1 coord2 bincontent\n";
  
  //TGraph2D *g2 = new TGraph2D(h2);
  
  
  //double z1=0;
  
  /*ofstream myfile_2;
  myfile_2.open ("x-y-z-w.txt");
  myfile_2 << "x y z w\n";
  
  TGraph2D *g2 = new TGraph2D(h2);
  g2->SetTitle("Graph title; X axis title; Y axis title; Z axis title");
  
  for (Int_t i = 0; i<g2->TGraph2D::GetN(); i++) {
        
        myfile_2 <<g2->TGraph2D::GetX()[i]<< ","<< g2->TGraph2D::GetY()[i] <<","<< g2->TGraph2D::GetZ()[i]<<"\n"; 
        
        
  }
  
  myfile_2.close();*/
  
  //auto c1 = new TCanvas("c1","c1",1000, 1000);
  //g2 -> Draw("SURF");
  //z=g2->TGraph2D::Interpolate(x,y);
  auto c1 = new TCanvas("c1","c1",1000, 1000);
  //gStyle->SetPalette(1);
  //dt->Draw("surf1");
  gStyle->SetNumberContours(10);
  gStyle->SetPalette(55);
  //h2 -> Draw("LEGO2Z");
  h2 -> Draw("COLZ");
  
  auto c2 = new TCanvas("c2","c2",1000, 1000);
  gStyle->SetNumberContours(10);
  gStyle->SetPalette(55);
  h2 -> Draw("LEGO2Z");

}  
