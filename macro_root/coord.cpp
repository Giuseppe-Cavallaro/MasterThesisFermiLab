#define NFolders 1
#define NTrees 1

//Change this to change the used folder
#define NPath 0
#define SelectedTree 0



void TreePlot(){
std::string FolderNames[NFolders] = {"flashmatchAna"};

std::string xname = "Middlex";
std::string yname = "Middley";
std::string zname = "Middlez";
std::string weight = "GammaStep";
std::string gamma = "GammaScint";
//std::string gamma2 = "TotGammaScint";

std::string FileName = "/media/giuseppe/01D5DC4DD5EAB920/lightmap/analisi_mu_30.root";
//std::string FileName = "/media/giuseppe/Seagate/lightmap/analisi_mu_59.root";
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

  std::vector<float> * v1=0;
  std::vector<float> * v2=0;
  std::vector<float> * v3=0;
  std::vector<float> * v4=0;
  std::vector<float> * v5=0;
  std::vector<float> * v6=0;
  Double_t x=0,y=0,z1=0,w=0,w2=0,w3=0;
  Int_t dim_dep=0;
  Int_t dim_x=0;
  
   T->SetBranchAddress(xname.c_str(),&v1);
  T->SetBranchAddress(yname.c_str(),&v2);
  T->SetBranchAddress(weight.c_str(),&v3);
  T->SetBranchAddress(gamma.c_str(),&v4);
  T->SetBranchAddress(zname.c_str(),&v5);
//  T->SetBranchAddress(gamma2.c_str(),&v6);
  
  Int_t Dimension = T->GetEntries();
  Int_t n = 0;
  
  ofstream myfile;
  myfile.open ("/home/giuseppe/Scrivania/map/map_30.txt");
  
  
  
  
  
  //}
  
  //myfile << "x y z w\n";
  //for (Int_t i = 0; i<Dimension; i++) {
  
  for (Int_t i = 0; i<Dimension; i++) {
    T -> GetEntry(i);
    n = v1->size();
    dim_dep=v4->size();
    dim_x=dim_dep/n;
    //w3=(v6->at(i));
    for (Int_t z = 0; z<n-2; z++) {
    x=v1->at(z); //x
    y=v2->at(z); //y
    z1=v5->at(z); //z
    w=(v3->at(z)); //gammastep
    w2=0;
    //w3=(v6->at(i));
    //w=(w/0.05)/w3;
    for(Int_t j = 0; j<dim_x; j++){
    w2+=v4->at(j+dim_x*z);}
    //w2=w2;
    w=(w/0.05)/w2;
    //w2=0;
    
    //h2->Fill(x,y,w);
    //h2->TH2::Interpolate(x,y,w);
    myfile <<x<<" "<< y <<" "<< z1 <<" "<< w <<"\n"; //write to file
   
  }
  
 }
 
 myfile.close();
 }
