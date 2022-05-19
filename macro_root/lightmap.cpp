#define NFolders 1
#define NTrees 1

//Change this to change the used folder
#define NPath 0
#define SelectedTree 0



void TreePlot(){
std::string FolderNames[NFolders] = {"flashmatchAna"};

std::string xname = "E_QL";
std::string yname = "TotEdep";
std::string zname = "Light";
std::string weight = "Fvis";
std::string gamma = "TrueCCNC";

//std::string FileName = "/media/giuseppe/01D5DC4DD5EAB920/lightmap/analisi_mu_33.root";
std::string FileName = "/home/giuseppe/Scrivania/try/analisi_mu_200_1.root";
//std::string FileName = "/media/giuseppe/01D5DC4DD5EAB920/lightmap/analisi_mu_14500.root";
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

  /*std::vector<float> * v1=0;
  std::vector<float> * v2=0;
  std::vector<float> * v3=0;
  std::vector<float> * v4=0;
  std::vector<float> * v5=0;*/
  
  float v1=0;
  float v2=0;
  float v3=0;
  float v4=0;
  int v5=0;
  Double_t x=0,y=0,z1=0,w=0,w2=0;
  
  T->SetBranchAddress(xname.c_str(),&v1);
  T->SetBranchAddress(yname.c_str(),&v2);
  T->SetBranchAddress(zname.c_str(),&v3);
  T->SetBranchAddress(weight.c_str(),&v4);
  T->SetBranchAddress(gamma.c_str(),&v5);
  
  Int_t Dimension = T->GetEntries();
  Int_t n = 0;
  
  //ofstream myfile;
  //myfile.open ("/home/giuseppe/Scrivania/mappus/54_500.txt");
  //myfile << "x y z w\n";
  
  TH1F *p = new TH1F("p","Profile",150,-2,2);
  TH1F *p1 = (TH1F*)(p->Clone());
  //&&(v4<(200*10^6))
  //Fvis!=0&&Light<200*10^6
  //EQ
  /*for (Int_t i = 0; i<Dimension; i++) {
    T -> GetEntry(i);
    x=v1;
    y=v2 ;
    w=(x-y)/y;
    cout<<x<<"e"<<y<<endl;  
    p->Fill(w);
    if(v5==0){
    x=v1;
    y=v2 ;
    w=(x-y)/y;  
    p1->Fill(w);}*/
    //p->Fill(x,y,w);
    //h2->TH2::Interpolate(x,y,w);
    //myfile <<x<<" "<< y <<" "<< z1 <<" "<< w <<"\n"; //write to file
    for (Int_t i = 0; i<Dimension; i++) {
    T -> GetEntry(i);
    if(v4!=0&&v3<200000000){
    x=v1;
    y=v2 ;
    w=(x-y)/y;
    p->Fill(w);
    if(v5==0){
    x=v1;
    y=v2; 
    w=(x-y)/y;  
    p1->Fill(w);}}
    
  
 }
 gStyle->SetOptFit(111111);
 auto c1 = new TCanvas("c1","c1",1000, 1000);
  //gStyle->SetPalette(1);
  //dt->Draw("surf1");
  gStyle->SetNumberContours(10);
  gStyle->SetPalette(55);
  //h2 -> Draw("LEGO2Z");
  p -> Draw("COLZ");
  p1 -> Draw("SAMES");
  p -> Fit("gaus","","sames");
  p1 -> Fit("gaus","","sames");
  
 
 //myfile.close();
 }
