/**
 *  Author: Philipp Schlunder
 *  Licence: GPL v3
 */

#include <iostream>
#include <fstream>
#include <TRandom3.h>
// ROOT includes
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph2D.h>
#include <TH1.h>
#include <TH3.h>
#include <TLegend.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TTree.h>
// toy mc includes
#include "TEvent.h"
#include "TDetector.h"


#include "TBranch.h"
#include "TRandom2.h"
#include "TF3.h"
#include "TError.h"
#include "Fit/BinData.h"
#include "Fit/Fitter.h"
#include "Math/WrappedMultiTF1.h"

using namespace std;

vector<double> CrossPlane(vector<double> MP,
                vector<double> KL,
                vector<double> start,
                vector<double> direction);

void dumpVector(string phrase, vector<double> Vector);

void TestEvents() {
  double azimuthal = 45, polar = 50;
  vector<double> posA, posB, sensor;
  // Event Position
  posA.push_back(1);
  posA.push_back(2);
  posA.push_back(3);
  // Testing Position
  posB.push_back(5);
  posB.push_back(4);
  posB.push_back(8);
  // Sensor position
  sensor.push_back(2);
  sensor.push_back(0);
  sensor.push_back(1);
  
  TEvent *e1 = new TEvent();
  TEvent *e2 = new TEvent(10, posA, azimuthal);
  TEvent *e3 = new TEvent(10, posA, azimuthal, polar);
  cout << "--- Testing TEvent ---" << endl;
  cout << "Constructor default: " << endl;
  e1->Dump();
  cout << "Constructor 2d: " << endl;
  e2->Dump();
  cout << "Constructor 3d: " << endl;
  e3->Dump();
  cout << "Distance between posA and posB: d = "
       << e3->GetDistance(posA, posB) << endl;
  cout << "Energy at (" << posB.at(0) << ", "
                        << posB.at(1) << ", "
                        << posB.at(2) << ") = "
                        << e3->GetEnergy(posB) << endl;
  cout << "Perpendicular foot from (" << sensor.at(0) << ", "
       << sensor.at(1) << ", "
       << sensor.at(2) << ") to the event track of the 3d event:" << endl;
  cout << "pF = (" << e3->GetPerpFoot(sensor).at(0) << ", "
       << e3->GetPerpFoot(sensor).at(1) << ", "
       << e3->GetPerpFoot(sensor).at(2) << ")." << endl;
}

void TestDetector() {
  TRandom3 *grandom = new TRandom3(0);
  //TDetector *d1 = new TDetector();
  //TDetector *d2 = new TDetector(9, 8);
  //TDetector *d3 = new TDetector(9, 8, 0.5);
  TDetector *detector = new TDetector("../detector_config.txt");
  /*
  cout << "Constructor default:" << endl;
  d1->Dump();
  cout << "Constructor with standard permeability." << endl;
  d2->Dump();
  cout << "Constructor with predefined permeability." << endl;
  d3->Dump();*/
  //TDetector *detector = new TDetector(9, 8);
  //cout << "Edge length: " << detector->GetEdgeLength() << endl;
  /*
  vector<double> test;
  test.push_back(2);
  test.push_back(1);
  test.push_back(1);
  TSensor *sens = new TSensor(test, 0.2);
  detector->Dump();
  detector->AddSensor(*sens);*/
  
  unsigned int temp, NoE;
  temp = detector->GetNumberOfGroups();
  NoE = detector->GetGroup(1).size();
  cout << "Number of groups: " << temp << endl;
  cout << "Number of Elements: " << NoE << endl;
  for (unsigned int i = 0; i < temp; i++) {
    cout << "Id = " << i+1 << endl;
    for (unsigned int j = 0; j < detector->GetGroup(i+1).size(); j++) {
      detector->GetGroup(i+1).at(j)->Dump();
    }
  }

  vector<double> position, sensor;
  // Event Position
  position.push_back(1);
  position.push_back(0);
  position.push_back(0);
  TEvent *event = new TEvent(100, position, 0, 0);
  detector->Detect(*event, grandom->Rndm());
  cout << "Intensity = "
            << detector->GetSensors().at(0)->GetIntensity() << endl;
}

void TestGroupSensors() {
  TDetector *d = new TDetector("../detector_config.txt");
  vector<TSensor*> sensors_;
  sensors_ = d->GetSensors();
  vector<vector<TSensor*> > groups_;
  cout << "Loaded all sensors." << endl;

  // loop over sensors to find the max id
  unsigned int maxId = 0;
  for(unsigned int i = 0; i < sensors_.size(); i++) {
    if (sensors_.at(i)->GetGroup() > maxId) {
      maxId = sensors_.at(i)->GetGroup();
    }
  }
  cout << "max id = " << maxId << endl;
  // set group_ size to max id
  groups_.resize(maxId);
  // insert dummy vector with dummy sensor into groups_
  for (unsigned int i = 0; i < groups_.size(); i++) {
    vector<TSensor*> newId;
    newId.push_back(new TSensor());
    groups_.at(i) = newId;
  }

  // actual sorting process 
  for (unsigned int i = 0; i < sensors_.size(); i++) {
    unsigned int id = sensors_.at(i)->GetGroup();
    groups_.at(id - 1).push_back(sensors_.at(i));
  }
  // erase dummy vectors and sensors
  for (unsigned int i = 0; i < groups_.size(); i++) {
    groups_.at(i).erase(groups_.at(i).begin());
  }
}

void ViewData() {
  TFile *file = new TFile("result.root", "READ");
  TTree *tree;
  file->GetObject("run", tree);

  Double_t qtotal;
  tree->SetBranchAddress("Q_total", &qtotal);

  for (Int_t i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    cout << "Q_total = " << qtotal << endl;
  }

  /*
  vector<double> arrivalTimes;
  tree->SetBranchAddress("ArrivalTime", &arrivalTimes);
  cout << "Arrival Times from Measurement.root" << endl;
  for (unsigned int i = 0; i < 10; i++) {
    cout << "arrival time number " << i << ": " << arrivalTimes.at(i) << endl;
  }
  */

  //TCanvas *c1 = new TCanvas("c1", "Test", 800, 600);
  //c1->cd();
}

void ThreeDFit() {
     
  const int n = 1000; 
  double x[n], y[n], z[n], v[n]; 
  double ev = 0.1;

  // generate the data
  TRandom2 r; 
  for (int i = 0; i < n; ++i) { 
     x[i] = r.Uniform(0,10);
     y[i] = r.Uniform(0,10);
     z[i] = r.Uniform(0,10); 
     v[i] = sin(x[i] ) + cos(y[i]) + z[i] + r.Gaus(0,ev);         
  }

  // create a 3d binned data structure
  ROOT::Fit::BinData data(n,3); 
  double xx[3];
  for(int i = 0; i < n; ++i) {
     xx[0] = x[i]; 
     xx[1] = y[i]; 
     xx[2] = z[i]; 
     // add the 3d-data coordinate, the predictor value (v[i])  and its errors
     data.Add(xx, v[i], ev); 
  }

  TF3 * f3 = new TF3("f3","[0] * sin(x) + [1] * cos(y) + [2] * z",0,10,0,10,0,10);
  f3->SetParameters(2,2,2);
  ROOT::Fit::Fitter fitter;
  // wrapped the TF1 in a IParamMultiFunction interface for teh Fitter class
  ROOT::Math::WrappedMultiTF1 wf(*f3,3);
  fitter.SetFunction(wf); 
  //
  bool ret = fitter.Fit(data); 
  if (ret) { 
     const ROOT::Fit::FitResult & res = fitter.Result(); 
     // print result (should be around 1) 
     res.Print(std::cout);
     // copy all fit result info (values, chi2, etc..) in TF3
     f3->SetFitResult(res);
     // test fit p-value (chi2 probability)
     double prob = res.Prob();
     if (prob < 1.E-2) 
        Error("exampleFit3D","Bad data fit - fit p-value is %f",prob);
     else
        std::cout << "Good fit : p-value  = " << prob << std::endl;

  }
  else 
     Error("exampleFit3D","3D fit failed");
}


int main()
{/*
  TRandom3 *grandom = new TRandom3(0);
  TDetector *d1 = new TDetector(100, 27);

  vector<double> position;
  position.push_back(2);
  position.push_back(3);
  position.push_back(4);
  TEvent *testEvent = new TEvent(10, position, 45, 46);
  double randomNumber = grandom->Rndm();
  d1->Detect(*testEvent, randomNumber);

  vector<TSensor*> sensors;
  sensors = d1->GetSensors();
  
  TCanvas *c1 = new TCanvas("c1", "c1", 700, 900);
  c1->cd();
  TFile f("result.root");
  TTree *tree = (TTree*)f.Get("run");
  
  tree->Draw("SensorPosX:SensorPosY:SensorPosZ:Intensity");
  c1->SaveAs("ModulTest.pdf");
  c1->Write();*/
  /*
  tree->Draw("Intensity");
  c1->SaveAs("ModulTest.pdf");
  c1->Write();*/
  vector<double> MP_, KL_, start_, direction_, result_, intersection_;
  MP_.push_back(10);
  MP_.push_back(5);
  MP_.push_back(5);
  KL_.push_back(0);
  KL_.push_back(10);
  KL_.push_back(10);
  start_.push_back(9);
  start_.push_back(7);
  start_.push_back(4);
  direction_.push_back(-1);
  direction_.push_back(0);
  direction_.push_back(0);

  result_ = CrossPlane(MP_, KL_, start_, direction_);
  for (int i = 0; i < 3; i++) { intersection_.push_back(result_.at(i+1)); }
  cout << "Abstand " << result_.at(0) << endl;
  dumpVector("Schnittpunkt", intersection_);

  return 0;
}

vector<double> CrossPlane(vector<double> MP,
                vector<double> KL,
                vector<double> start,
                vector<double> direction) {
  vector<double> result;
  result.resize(4);
  double gradient;
  bool intersection = false;
  cout << "in crossplane" << endl;
  for (int i = 0; i < 3; i++) {
    double length = KL.at(i);
    if (length == 0) {
      cout << "kantenlänge ist 0" << endl;
      double dLinePlane = MP.at(i) - start.at(i);
      gradient = dLinePlane / direction.at(i);
      cout << "Gradient ist " << gradient << endl;
      if (gradient != 0) {
        intersection = true;
        cout << "Schnittpunkt: " << endl;
        for (int j = 1; j < 4; j++) {
          result.at(j) = start.at(j-1) + direction.at(j-1) * gradient;
          cout << "Pos " << j << " ist " << result.at(j) << endl;
        }
        result.at(0) = sqrt(dLinePlane * dLinePlane);
      }
    }
  }
  for (int i = 0; i < 3; i++) {
    if ((result.at(i+1) > (MP.at(i) + KL.at(i) * 0.5))
      || (result.at(i+1) < (MP.at(i) - KL.at(i) * 0.5))) {
      intersection = false;
    }
  }
  if (!intersection) {
    for (int j = 1; j < 4; j++) {
      result.at(j) = -1;
    }
  }

  return result;
}

void dumpVector(string phrase, vector<double> Vector) {
  cout << phrase << " (" << Vector.at(0) << ", "
       << Vector.at(1) << ", "
       << Vector.at(2) << ")" << std::endl;
}
