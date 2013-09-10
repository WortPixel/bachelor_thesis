/**
 *  Author: Philipp Schlunder
 *  Licence: GPL v3
 */

#include <iostream>
#include <fstream>
#include <string>
// boost includes
#include <boost/math/constants/constants.hpp>
#include <boost/program_options.hpp>
#include <boost/random.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/normal_distribution.hpp>
// ROOT includes
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph2D.h>
#include <TH1.h>
#include <TH3.h>
#include <TLegend.h>
#include <TPolyMarker3D.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TTree.h>
// toy mc includes
#include "TDetector.h"
#include "TEvent.h"

namespace po = boost::program_options;

// structure for storing distribution parameters
struct distribution {
  std::string type;
  double mean, RMS;
};

/**
 * @brief DegToRad calculates the radiant of a given degree
 * @param degree the value to be transformened into a radiant
 * @return the radiant to a given degree
 */
double DegToRad(double degree);

/**
 * @brief GetNewEnergy calculates an energy value that belongs tu an exponential
 * distribution with the negative exponent -energyExponent_. The distribution
 * starts at energyMin and has an open ending.
 * @param energyMin sets the lower energy border of the energy spectrum
 * @param energyExponent defines the negative exponent of the energy spectrum
 * @param random is a pointer to the seeded random number generator
 * @return a randomly generated energy value in the given energy spectrum
 */
double GetNewEnergy(double energyMin, double energyMax, double energyExponent,
                    TRandom3 *random);

/**
 * @brief GetNewTrack creates random distributed values that are needed for a
 * track
 * @param detectorDimensions contains all three dimensions of a detector
 * @param rng awaits the address of a random number generator
 * @return a 5 dimensional vector, the first 3 dimensions are the cartesian
 * spacial dimensions, the 4th is the aziumthal angle and the 5th is the polar.
 */
std::vector<double> GetNewTrack(std::vector<double> detectorDimensions,
                                boost::mt19937& rng);

/**
 * @brief LoadConfig reads the parameters of a given config file and sets the
 * gloabl variables
 * @desc_ holds the assignment between variables and option names
 * @vm_ is a map that stores the read options
 */
void LoadConfig(po::options_description& desc_, po::variables_map& vm_);


/**
 * @brief PlotDistributions creates plots of the different distributed values:
 * azimuthalAngle, polarAngle, Energy.
 * @param samples awaits a pointer to all created events
 */
void PlotDistributions(std::vector<TEvent*>& samples);

/**
 * @brief PlotDetector creates a plot of the sensor positions of the given
 * detector into plots/Detector.pdf
 * @param detector is the detector whos sensors positions are plotted
 */
void PlotDetector(TDetector &detector);

/**
 * @brief PlotIntensity creates a 3 dimensional plot with colored spheres at the
 * detector positions. Intensities are indicated via colors.
 */
void PlotIntensity();

/**
 * @brief PlotResults creates an event view of the deposited intensities
 * @param event contains the number of the event that should be displayed
 */
void PlotResults(unsigned int event, 
                 std::vector<TEvent*>& samples,
                 TDetector &detector,
                 std::vector<std::vector<double> > results);

/**
 * @brief SphericalToCartesian converts spherical coordiantes to cartesian
 * @param radius conatins the radial distance of the position from the center
 * @param azimuthal contains the azimuthal deflection
 * @param polar contains the polar deflection
 * @return cartesian coordinates of the given spherical position
 */
std::vector<double> SphericalToCartesian(double radius, double azimuthal,
                                         double polar);

/**
 * @DumpVector is a debugging method that couts a vector with a string
 * @param phrase awaits a string that describes the vector
 * @param position is the vector that should be displayed
 */
void DumpVector(std::string phrase, std::vector<double> position) {
  std::cout << phrase << " (" << position.at(0) << ", "
            << position.at(1) << ", "
            << position.at(2) << ")" << std::endl;
}

int main(int argc, char const *argv[]) {
  std::string version = "Version 0.9 from 2013-06-29 last edited by Philipp";
  // set global style for root
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1);
  // --- start: initializing parameters and setting default values ---
  // output:
  std::string output_ = "result";
  // detector:
  double edgeLength_ = 10;
  int numberOfSensors_ = 27;
  
  // distributions:
  int seed_ = 0;
  double energyExponent_ = 2.7; // negative exponent of the power function
  double energyMin_ = 1e1; // Minimum of the generated energy spectrum in GeV
  double energyMax_ = 1e5; // Maximum of the generated energy spectrum in GeV
  distribution azimuthalAngle_, polarAngle_, xAxis_, yAxis_, zAxis_;
  azimuthalAngle_.type = "normal";
  azimuthalAngle_.mean = 0;  // in degree
  azimuthalAngle_.RMS = 5;  // in degree
  polarAngle_.type = "normal";
  polarAngle_.mean = 0;  // in degree
  polarAngle_.RMS = 5;  // in degree
  xAxis_.type = "normal";
  xAxis_.mean = edgeLength_/2;  // in degree
  xAxis_.RMS = 1; // in degree
  yAxis_.type = "normal";
  yAxis_.mean = edgeLength_/2;  // in degree
  yAxis_.RMS = 1; // in degree
  zAxis_.type = "normal";
  zAxis_.mean = 0;  // in degree
  zAxis_.RMS = 0; // in degree

  int amountOfSamples_ = 1e4;
  int measuredSamples_ = 0;

  TDetector *detector_; // addresses the used detector
  std::vector<TEvent*> samples_;  // can store all created events
  // --- end: initializing parameters and setting default values ---

  // --- start: assigning and loading values out of config ---
  po::options_description desc_("Available parameters for the config.ini:");
  desc_.add_options()
    ("settings.output", po::value<std::string>(&output_), "output filename")
    ("settings.seed", po::value<int>(&seed_),
     "seed for the random number generators")
    ("settings.edgeLength", po::value<double>(&edgeLength_),
     "edge length of the detector")
    ("settings.numberOfSensors", po::value<int>(&numberOfSensors_),
     "number of sensors the detector should hold")
    ("settings.amountOfSamples", po::value<int>(&amountOfSamples_),
     "amount of samples that should be created")
    ("energy.exponent", po::value<double>(&energyExponent_),
     "negative value of exponent is used as exponent for energy distribution")
    ("energy.min", po::value<double>(&energyMin_),
     "minimum of the energy spectrum")
    ("energy.max", po::value<double>(&energyMax_),
     "maximum of the energy spectrum")
    ("azimuthalAngle.type",
     po::value<std::string>(&azimuthalAngle_.type),
     "distribution type for the azimuthal angle")
    ("azimuthalAngle.mean", po::value<double>(&azimuthalAngle_.mean),
     "mean of the azimuthal angle distribution")
    ("azimuthalAngle.RMS", po::value<double>(&azimuthalAngle_.RMS),
     "root mean square of the azimuthal angle distribution")
    ("polarAngle.type", po::value<std::string>(&polarAngle_.type),
     "distribution type for the polar angle")
    ("polarAngle.mean", po::value<double>(&polarAngle_.mean),
     "mean of the polar angle distribution")
    ("polarAngle.RMS", po::value<double>(&polarAngle_.RMS),
     "root mean square of the polar angle distribution")
    ("xAxis.type", po::value<std::string>(&xAxis_.type),
     "distribution type of the x-axis")
    ("xAxis.mean", po::value<double>(&xAxis_.mean),
     "mean of the x-axis distribution")
    ("xAxis.RMS", po::value<double>(&xAxis_.RMS),
     "root mean square of the x-axis distribution")
    ("yAxis.type", po::value<std::string>(&yAxis_.type),
     "distribution type of the y-axis")
    ("yAxis.mean", po::value<double>(&yAxis_.mean),
     "mean of the y-axis distribution")
    ("yAxis.RMS", po::value<double>(&yAxis_.RMS),
     "root mean square of the y-axis distribution")
    ("zAxis.type", po::value<std::string>(&zAxis_.type),
     "distribution type of the z-axis")
    ("zAxis.mean", po::value<double>(&zAxis_.mean),
     "mean of the z-axis distribution")
    ("zAxis.RMS", po::value<double>(&zAxis_.RMS),
     "root mean square of the z-axis distribution");

  po::options_description generic_("Generic options");
  generic_.add_options()
    ("version,v", "print version string")
    ("help,h", "produce help message");
  po::variables_map vmGeneric_;
  po::store(po::parse_command_line(argc, argv, generic_), vmGeneric_);
  po::notify(vmGeneric_);
  if (vmGeneric_.count("help")) {
    std::cout << generic_ << std::endl;
    return 1;
  }
  if (vmGeneric_.count("version")) {
    std::cout << version << std::endl;
    return 1;
  }
  // create map that will store the loaded values
  po::variables_map vm_;
  // load the assigned parameters of the config file into the map vm_
  LoadConfig(desc_, vm_);
  // --- end: assigning and loading values out of config ---

  // --- start: creating and ploting detector ---
  // create individual detector from detector_config if edgeLength and
  // numberOfSensors is 0
  if (edgeLength_ == 0 && numberOfSensors_ == 0) {
    detector_ = new TDetector("../detector_config.txt");
    edgeLength_ = detector_->GetEdgeLength();
    numberOfSensors_ = detector_->GetNumberOfSensors();
    std::cout << "Created detector from detector_config.txt." << std::endl;
  }
  // or just create a detector with equally distributed sensors
  else {
    detector_ = new TDetector(edgeLength_, numberOfSensors_);
    std::cout << "Created detector with an edge length of " << edgeLength_
              << " meters and " << numberOfSensors_ << " sensors." << std::endl;
  }
  //detector_->Dump();  // textual output of sensor positions and detector data
  // create a plot of the detectors sensor positions
  PlotDetector(*detector_);
  // --- end: creating and ploting detector ---

  // --- start: preparing distributions ---
  /**
    * Generating random numbers using the mersenne_twister algorithm provided
    * by boost and implemented after "Mersenne Twister: A 623-dimensionally
    * equidistributed uniform pseudo-random number generator", Makoto Matsumoto
    * and Takuji Nishimura, ACM Transactions on Modeling & Computer Simulation:
    * Special Issue on Uniform Random Number Generation, Vol. 8, No. 1, January
    * 1998, pp. 3-30.
    */
  boost::mt19937 rng(seed_);
  TRandom3 *grandom = new TRandom3(seed_);
  /**
   * prepare normal distributions for angle distributions
   * function boost uses:
   * p(x) = 1 / (sqrt(2 * pi * RMS)) * exp(- (x - mean) **2 / (2 * RMS ** 2))
   * where a normal_distribution awaits mean and RMS as parameters
   */
  boost::normal_distribution<> azimuthalAngleDistribution(
    DegToRad(azimuthalAngle_.mean),
    DegToRad(azimuthalAngle_.RMS));
  boost::normal_distribution<> polarAngleDistribution(
    DegToRad(polarAngle_.mean),
    DegToRad(polarAngle_.RMS));
  boost::normal_distribution<> xAxisDistribution(xAxis_.mean, xAxis_.RMS);
  boost::normal_distribution<> yAxisDistribution(yAxis_.mean, yAxis_.RMS);
  boost::normal_distribution<> zAxisDistribution(zAxis_.mean, zAxis_.RMS);
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> >
    GetNewAzimuthalAngle(rng, azimuthalAngleDistribution);
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> >
    GetNewPolarAngle(rng, polarAngleDistribution);
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> >
    GetNewXPosition(rng, xAxisDistribution);
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> >
    GetNewYPosition(rng, yAxisDistribution);
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> >
    GetNewZPosition(rng, zAxisDistribution);
  std::cout << "Loaded and prepared distributions." << std::endl;
  // --- end: preparing distributions ---

  // --- start: initializing parameters for detection ---
  int eventId = 0;
  std::vector<double> position; // starting position of event
  // set default value for event starting position in case none is set in config
  for (unsigned int i = 0; i < 3; i++) { position.push_back(0); }
  double energy, azimuthalAngle, polarAngle;  // event attributes
  double randomNumber = grandom->Rndm();
  double activatedSensors, activatedGroups, totalIntensity;
  double firstSensorX, firstSensorY, firstSensorZ;
  double centerOfIntensityX, centerOfIntensityY, centerOfIntensityZ;
  double effectiveTrackLength;
  double reconstructedAzimuthalAngle, reconstructedPolarAngle;
  // stores activated sensors per id
  std::vector<unsigned int> activatedSensorsPerGroup;
  std::vector<double> arrivalTime;  // stores arrival times of all sensors
  std::vector<double> sensorPosX, sensorPosY, sensorPosZ;
  std::vector<double> intensity_;  // stores detected intensities for one event
  std::vector<std::vector<double> > results_; // store all detected intensities
  
  // add file ending to output filename
  output_ += ".root";
  TFile *rootFile = new TFile(output_.c_str(), "RECREATE");
  TTree *data = new TTree("run", "raw data of an event");
  
  // assign placeholders with branch attributes
  data->Branch("1stSensor_XPosition", &firstSensorX,
               "1stSensor_XPosition/D");
  data->Branch("1stSensor_YPosition", &firstSensorY,
               "1stSensor_YPosition/D");
  data->Branch("1stSensor_ZPosition", &firstSensorZ,
               "1stSensor_ZPosition/D");
  data->Branch("ArrivalTime", &arrivalTime);
  data->Branch("CoI_XPosition", &centerOfIntensityX, "CoI_XPosition/D");
  data->Branch("CoI_YPosition", &centerOfIntensityY, "CoI_YPosition/D");
  data->Branch("CoI_ZPosition", &centerOfIntensityZ, "CoI_ZPosition/D");
  data->Branch("Event_ID", &eventId, "Event_ID/Is");
  data->Branch("Event_Seed", &seed_, "Event_Seed/I");
  data->Branch("Effective_Track_Length", &effectiveTrackLength,
               "Effective_Track_Length/D");
  data->Branch("Intensity", &intensity_);
  data->Branch("MC_Energy0", &energy, "MC_Energy0/D");
  data->Branch("MC_AzimuthalAngle", &azimuthalAngle,
               "MC_AzimuthalAngle/D");
  data->Branch("MC_PolarAngle", &polarAngle, "MC_PolarAngle/D");
  data->Branch("MC_XPosition", &position.at(0), "MC_XPosition/D");
  data->Branch("MC_YPosition", &position.at(1), "MC_YPosition/D");
  data->Branch("MC_ZPosition", &position.at(2), "MC_ZPosition/D");
  data->Branch("N_Groups", &activatedGroups, "N_Groups/D");
  data->Branch("N_Sensors_per_Group", &activatedSensorsPerGroup);
  data->Branch("N_Sensors", &activatedSensors, "N_Sensors/D");
  data->Branch("Q_total", &totalIntensity, "TotalIntensity/D");
  data->Branch("R_AzimuthalAngle", &reconstructedAzimuthalAngle,
               "R_AzimuthalAngle/D");
  data->Branch("R_PolarAngle", &reconstructedPolarAngle, "R_PolarAngle/D");
  // needing sensor positions for each event for a plot
  data->Branch("SensorPosX", &sensorPosX);
  data->Branch("SensorPosY", &sensorPosY);
  data->Branch("SensorPosZ", &sensorPosZ);
  // --- end: initializing parameters for detection ---

  // --- start: measurement for all samples ---
  std::cout << "Starting measurement." << std::endl;
  for (int i = 0; i < amountOfSamples_; i++) {
    // reset event features
    bool measured = false;
    double minArrivalTime = 10000, maxArrivalTime = 0;
    TSensor *firstSensor, *lastSensor;

    activatedSensors = 0;
    activatedSensorsPerGroup.resize(0);
    arrivalTime.resize(0);
    intensity_.resize(0);
    sensorPosX.resize(0);
    sensorPosY.resize(0);
    sensorPosZ.resize(0);
    totalIntensity = 0;

    // increase event id
    eventId++;

    // --- start: generating distributed random numbers ---
    energy = GetNewEnergy(energyMin_, energyMax_, energyExponent_, grandom);

    azimuthalAngle = GetNewAzimuthalAngle();
    polarAngle = GetNewPolarAngle();
    position.at(0) = GetNewXPosition();
    // make sure that the position values are within the detector
    while ((position.at(0) > edgeLength_) || (position.at(0) < 0)) {
      position.at(0) = GetNewXPosition();
    }
    position.at(1) = GetNewYPosition();
    while ((position.at(1) > edgeLength_) || (position.at(1) < 0)) {
      position.at(1) = GetNewYPosition();
    }
    position.at(2) = GetNewZPosition();
    while ((position.at(2) > edgeLength_) || (position.at(2) < 0)) {
      position.at(2) = GetNewZPosition();
    }
    // --- end: generating distributed random numbers ---

    // create event for detection
    TEvent *event = new TEvent();
    // fill its values with just generated random numbers
    event->SetEnergy(energy);
    event->SetAzimuthalAngle(azimuthalAngle);
    event->SetPolarAngle(polarAngle);
    event->SetPosition(position);
    // and store it for further handling
    samples_.push_back(event);

    // get new random number for detection
    randomNumber = grandom->Rndm();

    // detect the given event -> sensors contain measured intensity
    detector_->Detect(*event, randomNumber);

    // prepare activatedSensorsPerGroup
    for (unsigned int j = 0; j < detector_->GetNumberOfGroups(); j++) {
      activatedSensorsPerGroup.push_back(0);
    }

    // --- start: generating features out of measured data ---
    // calculate the center of intensity for this event
    centerOfIntensityX = detector_->GetCenterOfIntensity().at(0);
    centerOfIntensityY = detector_->GetCenterOfIntensity().at(1);
    centerOfIntensityZ = detector_->GetCenterOfIntensity().at(2);

    // read data from sensors and get various features
    for(int j = 0; j < numberOfSensors_; j++) {
      // get intensity of the addressed sensor
      double intensity = detector_->GetSensors().at(j)->GetIntensity();
      // store arrival time of addressed sensor
      double tempArrivalTime = detector_->GetSensors().at(j)->GetArrivalTime();
      // increase the counter for activated sensors if an intensity is detected
      if (intensity > 0) {
        if (!measured) {
          measured = true;
        }
        activatedSensors++;
        // and compare it with the temporary minimal arrival time
        // a value equal to zero means, the detector did not detect
        if ((tempArrivalTime != 0) && (tempArrivalTime < minArrivalTime)) {
          // replace minArrivalTime with current
          minArrivalTime = tempArrivalTime;
          // store z position = height of current sensor
          firstSensorX = detector_->GetSensors().at(j)->GetPosition().at(0);
          firstSensorY = detector_->GetSensors().at(j)->GetPosition().at(1);
          firstSensorZ = detector_->GetSensors().at(j)->GetPosition().at(2);
          firstSensor = detector_->GetSensors().at(j);
        }
        if ((tempArrivalTime != 0) && (tempArrivalTime > maxArrivalTime)) {
          // replace maxArrivalTime with current
          maxArrivalTime = tempArrivalTime;
          lastSensor = detector_->GetSensors().at(j);
        }

        // increase number of activated sensor counter for this group
        activatedSensorsPerGroup.at(detector_->GetSensors().at(j)->GetGroup()-1)
                                    += 1;
      }
      // store activation time of the addressed sensor
      arrivalTime.push_back(tempArrivalTime);
      // collect the total deposited intensity
      totalIntensity += intensity;
      // add given intensity to intensity-collection for this event
      intensity_.push_back(intensity);
      sensorPosX.push_back(detector_->GetSensors().at(j)->GetPosition().at(0));
      sensorPosY.push_back(detector_->GetSensors().at(j)->GetPosition().at(1));
      sensorPosZ.push_back(detector_->GetSensors().at(j)->GetPosition().at(2));
    }

    activatedGroups = 0;
    // collect data about number of activated groups
    for (unsigned int j = 0; j < activatedSensorsPerGroup.size(); j++) {
      if (activatedSensorsPerGroup.at(j) > 0) { activatedGroups += 1; }
    }

    // --- start: generate features of reconstructed track ---
    TEvent reconstructedTrack = detector_->ReconstructTrack();
    reconstructedAzimuthalAngle = reconstructedTrack.GetAzimuthalAngle();
    reconstructedPolarAngle = reconstructedTrack.GetPolarAngle();
    // first and last sensor position projected on track
    std::vector<double> firstPosition, lastPosition;
    if (measured) {
      measuredSamples_ += 1;
      firstPosition = reconstructedTrack.GetPerpFoot(firstSensor->GetPosition());
      lastPosition = reconstructedTrack.GetPerpFoot(lastSensor->GetPosition());
      effectiveTrackLength = reconstructedTrack.GetDistance(firstPosition,
                                                            lastPosition);
    }
    // --- end: generate features of reconstructed track ---

    // store measured intensities for this event
    results_.push_back(intensity_);

    // --- end: generating features out of measured data ---

    // add features for this event to the root tree
    data->Fill();
  }

  // write the tree into the root file and close it
  data->Write();
  delete data;
  rootFile->Close();
  // --- end: measurement for all samples ---

  std::cout << "Number of created events: " << samples_.size() << std::endl;
  std::cout << "Number of measured events: " << measuredSamples_ << std::endl;
  // generate plots for the distributed values of the events
  PlotDistributions(samples_);
  PlotResults(0, samples_, *detector_, results_);
  PlotIntensity();

  return 0;
}

double DegToRad(double degree) {
  // load pi from boosts constants with double precision
  double pi = boost::math::constants::pi<double>();

  return (degree / 180) * pi;
}

double GetNewEnergy(double energyMin, double energyMax, double energyExponent,
                    TRandom3 *random) {
  /**
   * generate a random number y_i that follows an exponential distribution g(x)
   * and has the probability density g(y) whereas x is a uniformly distributed
   * random number with a probability density f(x) that is 1 for 0 <= x < 1 and
   * f(x) = 0 for x < 0 or x >= 1.
   * Using g(y) = | dx / dy | * f(x) one receivs
   * x = G(y) = \int \limits_{-\infty}^{y} g(t) dt which leads to
   * y = G^(-1)(x) where G^(-1) generates exponential distributed values.
   * inspired by Prof. Dr. Dr. Wolfgang Rhodes script about "statistical methods
   * for data analysis" page 41, from 12th october 2012.
   */
  double base, exponent;
  exponent = 1 - energyExponent;
  base = (pow(energyMax, exponent) - pow(energyMin, exponent)) * random->Rndm()
       + pow(energyMin, exponent);

  return pow(base, 1/exponent);
}

std::vector<double> GetNewTrack(std::vector<double> detectorDimensions,
                                boost::mt19937& rng) {
  boost::uniform_01<> uniformDistribution;
  boost::variate_generator<boost::mt19937&, boost::uniform_01<> >
  GetUniform(rng, uniformDistribution);
  double pi = boost::math::constants::pi<double>();
  boost::uniform_real<> azimuthalDistribution((-1)*pi, pi);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
  GetAzAngle(rng, azimuthalDistribution);
  boost::uniform_real<> polarDistribution(0, pi);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
  GetPoAngle(rng, polarDistribution);

  std::vector<double> start, ending, result;
  result.resize(5);
  double radius = detectorDimensions.at(0);
  for (unsigned int i = 0; i < 3; i++) {
    if (detectorDimensions.at(i) > radius) {
      radius = detectorDimensions.at(i);
    }

    double randomNumber = GetUniform();
    ending.push_back(randomNumber * detectorDimensions.at(i));
  }
  double a, b;
  a = GetAzAngle();
  b = GetPoAngle();
  start = SphericalToCartesian(radius, a, b);
  for (unsigned int i = 0; i < 3; i++) {
    result.at(i) = ending.at(i) - start.at(i);
  }
  double newAzimuthal, newPolar;
  newPolar = acos(result.at(2) / radius);
  newAzimuthal = acos(result.at(0) / (radius * sin(newPolar)));
  result.at(0) = ending.at(0);
  result.at(1) = ending.at(1);
  result.at(2) = ending.at(2);
  result.at(3) = newAzimuthal;
  result.at(4) = newPolar;

  return result;
}

void LoadConfig(po::options_description& desc_, po::variables_map& vm_) {
  //--- reading settings ---
  // open config file
  std::ifstream config_file("../config.ini");
  // clear parameter storage before reading
  vm_ = po::variables_map();
  // read parameters from file into map
  po::store(po::parse_config_file(config_file, desc_), vm_);
  config_file.close();
  po::notify(vm_);
}

void PlotDistributions(std::vector<TEvent*>& samples) {
  double pi = boost::math::constants::pi<double>();

  TFile *f1 = new TFile("../plots/AzimuthalDistribution.root", "recreate");
  TCanvas *c1 = new TCanvas("c1", "Testing plots", 800, 600);
  c1->cd();
  TH1D *hAzimuthal = new TH1D("hAzimuthal", "AzimuthalAngle", 20, -pi/2, pi/2);
  for (unsigned int i = 0; i < samples.size(); ++i) {
    hAzimuthal->Fill(samples.at(i)->GetAzimuthalAngle());
  }
  hAzimuthal->GetXaxis()->SetTitle("Azimuthal Angle [rad]");
  hAzimuthal->GetYaxis()->SetTitle("Amount of occurence");
  hAzimuthal->SetTitle("Azimuthal angle distribution of created events");
  hAzimuthal->Draw();
  c1->SaveAs("../plots/AzimuthalDistribution.pdf");
  c1->Write();
  f1->Write();
  f1->Close();

  TFile *f2 = new TFile("../plots/PolarDistribution.root", "recreate");
  TCanvas *c2 = new TCanvas("c2", "Testing plots", 800, 600);
  c2->cd();
  TH1D *hPolar = new TH1D("hPolar", "PolarAngle", 20, -pi, pi);
  for (unsigned int i = 0; i < samples.size(); ++i) {
    hPolar->Fill(samples.at(i)->GetPolarAngle());
  }
  hPolar->GetXaxis()->SetTitle("Polar Angle [rad]");
  hPolar->GetYaxis()->SetTitle("Amount of occurence");
  hPolar->SetTitle("Polar angle distribution of created events");
  hPolar->Draw();
  c2->SaveAs("../plots/PolarDistribution.pdf");
  c2->Write();
  f2->Write();
  f2->Close();

  TFile *f3 = new TFile("../plots/EnergyDistribution.root", "recreate");
  TCanvas *c3 = new TCanvas("c3", "Testing plots", 800, 600);
  c3->cd();
  c3->SetLogy();
  TH1D *hEnergy = new TH1D("hEnergy", "Energy", 40, 0.0, 100.0);
  for (unsigned int i = 0; i < samples.size(); ++i) {
    hEnergy->Fill(samples.at(i)->GetEnergy());
  }
  hEnergy->GetXaxis()->SetTitle("Energy [GeV]");
  hEnergy->GetYaxis()->SetTitle("log(Amount of occurence)");
  hEnergy->SetTitle("Energy distribution of created events");
  hEnergy->Draw();
  c3->SaveAs("../plots/EnergyDistribution.pdf");
  c3->Write();
  f3->Write();
  f3->Close();
  delete c1;
  delete c2;
  delete c3;
  delete f1;
  delete f2;
  delete f3;
}

void PlotDetector(TDetector &detector) {
  TFile *f1 = new TFile("../plots/Detector.root", "recreate");
  TCanvas *c1 = new TCanvas("viewDetector", "Detector", 800, 600);
  c1->cd();
  double sensorX, sensorY, sensorZ;
  TGraph2D *g1 = new TGraph2D();
  for (unsigned int k = 0; k < detector.GetSensors().size(); k++) {
    sensorX = detector.GetSensors().at(k)->GetPosition().at(0);
    sensorY = detector.GetSensors().at(k)->GetPosition().at(1);
    sensorZ = detector.GetSensors().at(k)->GetPosition().at(2);
    g1->SetPoint(k, sensorX, sensorY, sensorZ);
  }
  // need to set the range for all axis
  g1->SetTitle("Sensor placement");
  g1->GetXaxis()->SetLimits(0., detector.GetEdgeLength());
  g1->GetXaxis()->SetTitle("x");
  g1->GetYaxis()->SetLimits(0., detector.GetEdgeLength());
  g1->GetYaxis()->SetTitle("y");
  g1->GetZaxis()->SetLimits(0., detector.GetEdgeLength());
  g1->GetZaxis()->SetTitle("z");
  g1->Draw("P0");
  c1->SaveAs("../plots/Detector.pdf");
  c1->Write();
  f1->Write();
  f1->Close();
  delete c1;
}

void PlotIntensity() {
  TCanvas *c1 = new TCanvas("c1", "c1", 900, 900);
  c1->cd();
  TFile f("result.root");
  TTree *tree = (TTree*)f.Get("run");

  tree->Draw("SensorPosX:SensorPosY:SensorPosZ:Intensity", "", "colz");
  c1->SaveAs("../plots/Intensity.pdf");
  c1->SaveAs("../plots/Intensity.root");
  f.Close();
  delete c1;
}

void PlotResults(unsigned int event,
                 std::vector<TEvent*>& samples,
                 TDetector &detector,
                 std::vector<std::vector<double> > results) {
  double edgeLength = detector.GetEdgeLength();
  TFile *f1 = new TFile("../plots/Measurement.root", "recreate");
  TCanvas *c1 = new TCanvas("c1", "Startingpoint", 800, 600);
  c1->cd();
  int BINS = 40;
  TH3D *hSensors = new TH3D("hSensors", "Sensors", BINS, 0.0, edgeLength, BINS,
                          0.0, edgeLength, BINS, 0.0, edgeLength),
         *hStart = new TH3D("hStart", "Start", BINS, 0.0, edgeLength, BINS, 0.0,
                          edgeLength, BINS, 0.0, edgeLength);
  for (unsigned int i = 0; i < detector.GetSensors().size(); i++) {
    double sensorX = detector.GetSensors().at(i)->GetPosition().at(0);
    double sensorY = detector.GetSensors().at(i)->GetPosition().at(1);
    double sensorZ = detector.GetSensors().at(i)->GetPosition().at(2);
    double sensorIntensity = results.at(event).at(i);
    hSensors->Fill(sensorX, sensorY, sensorZ, sensorIntensity);
  }
  for (unsigned int i = 0; i < samples.size(); ++i) {
    hStart->Fill(samples.at(i)->GetPosition().at(0),
                 samples.at(i)->GetPosition().at(1),
                 samples.at(i)->GetPosition().at(2));
  }
  hStart->GetXaxis()->SetTitle("x [m]");
  hStart->GetYaxis()->SetTitle("y [m]");
  hStart->GetZaxis()->SetTitle("z [m]");
  hStart->SetTitle("Measurement");
  hSensors->Draw("samelego");
  c1->SaveAs("../plots/Measurement.pdf");
  c1->Write();
  f1->Write();
  f1->Close();
  delete c1;
}

std::vector<double> SphericalToCartesian(double radius, double azimuthal,
                                         double polar) {
  std::vector<double> result;
  result.push_back(radius * cos(azimuthal) * sin(polar));
  result.push_back(radius * sin(azimuthal) * sin(polar));
  result.push_back(radius * cos(polar));

  return result;
}
