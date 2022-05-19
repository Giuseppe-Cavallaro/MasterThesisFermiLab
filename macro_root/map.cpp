#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>

using namespace std;

void TreePlot(){

ifstream file("/home/giuseppe/Scrivania/map/dati30000.txt");

  if (!file)
  {
    cerr << "cannot read the file example.txt :" 
         << strerror(errno) << endl;
    return -1;
  }

  vector<double> column1; 
  vector<double> column2; 
  vector<double> column3; 
  vector<double> column4; 
  double n1, n2, n3, n4,x,y,z,w;

  while (file >> n1 >> n2 >> n3 >> n4) {
    column1.push_back(n1);
    column2.push_back(n2);
    column3.push_back(n3);
    column4.push_back(n4);
  }

 TH3* h2 = new TH3F("h2", "X-Y-Z LightMap", 35, -350, 350, 35, -550, 550, 35, 0, 1350);
   h2->GetXaxis()->SetTitle("X [cm]");
   h2->GetYaxis()->SetTitle("Y [cm]");
   h2->GetZaxis()->SetTitle("Z [cm]");
   
  int Dimension=column1.size();
  //ofstream myfile;
  //myfile.open ("bin_31.txt");
  
 // for (Int_t i = 0; i<Dimension; i++) {
    
   // y= column2.at(i);
   // z= column3.at(i);
   // w= column4.at(i);

   // h2->Fill(z,y,w);
    //myfile <<y<<" "<< z <<" "<< w <<"\n"; //write to file
    //h2->TH2::Interpolate(x,y,w);
    //myfile_2 <<x<<","<< y <<","<< w <<"\n"; //write to file
   
  //}
  
   for (Int_t i = 0; i<Dimension; i++) {
    
    x=column1.at(i);
    y=column2.at(i);
    z=column3.at(i);
    w=(column4.at(i)*0.0125);
    h2->Fill(x,y,z,w);
    //h2->TH2::Interpolate(x,y,w);
    //myfile_2 <<x<<","<< y <<","<< w <<"\n"; //write to file
   
  }
  
 TAxis *xaxis = h2->GetXaxis();
 TAxis *yaxis = h2->GetYaxis();
 TAxis *zaxis = h2->GetZaxis();
 
 ofstream myfile;
 myfile.open ("/home/giuseppe/Scrivania/mappa30000.txt");
 //myfile << "x y z bincontent\n";
 
  for (int kbin=1; kbin<=zaxis->GetNbins(); ++kbin){
  for (int jbin=1; jbin<=yaxis->GetNbins(); ++jbin){
    for (int ibin=1; ibin<=xaxis->GetNbins(); ++ibin){
      //h_new->Fill(xaxis->GetBinCenter(ibin),yaxis->GetBinCenter(jbin),h->GetBinContent(ibin,jbin));
      myfile <<xaxis->GetBinCenter(ibin)<<" "<< yaxis->GetBinCenter(jbin) <<" "<< zaxis->GetBinCenter(kbin) <<" "<< h2->GetBinContent(ibin,jbin,kbin) <<"\n"; //write to file
    }
  }
  }
 
  myfile.close();
  
//auto c1 = new TCanvas("c1","c1",1000, 1000);
  //g2 -> Draw("SURF");
  //z=g2->TGraph2D::Interpolate(x,y);
  auto c1 = new TCanvas("c1","c1",1000, 1000);
  //gStyle->SetPalette(1);
  //dt->Draw("surf1");
  gStyle->SetNumberContours(10);
  gStyle->SetPalette(55);
  //h2 -> Draw("LEGO2Z");
  h2 -> Draw("BOX");

}    














