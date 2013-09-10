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

/**
 * @brief DegToRad calculates the radiant of a given degree
 * @param degree the value to be transformened into a radiant
 * @return the radiant to a given degree
 */
double DegToRad(double degree);

/**
 * @brief GetCrossPlane returns the intersections points between a given track
 * and a given plane.
 * @param centralPoint defines the centraly positioned point of a plane
 * @param edgeLength awaits the dimensions of the plane, use 0 for its thickness
 * @param start starting position of the track
 * @param direction direction of the track
 * @return a 4 dimensional vector. 1st dimension is the distance (plane, start),
 * where as the other 3 dimensions are the 3 cartesian coordiantes of the
 * intersection point
 */
std::vector<double> GetCrossPlane(std::vector<double> centralPoint,
                                  std::vector<double> edgeLength,
                                  std::vector<double> start,
                                  std::vector<double> direction);

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
 * @brief LoadConfig reads the parameters of a given config file and sets the
 * gloabl variables
 * @desc_ holds the assignment between variables and option names
 * @vm_ is a map that stores the read options
 */
void LoadConfig(po::options_description& desc_, po::variables_map& vm_);

std::vector<double> SphericalToCartesian(double r, double azimuthal,
                                         double polar);

std::vector<double> GetNewTrack(std::vector<double> detectorDimensions,
                                double azimuthalSpectrum,
                                double polarSpectrum,
                                boost::mt19937& rng);

/**
 * @brief dumpVector creates a well readable output of the given vector
 * @param phrase awaits a string describing the use of this vector
 * @param position awaits a 3 dimensional vector
 */
void dumpVector(std::string phrase, std::vector<double> position) {
  std::cout << phrase << " (" << position.at(0) << ", "
            << position.at(1) << ", "
            << position.at(2) << ")" << std::endl;
}

int main() {
  // set global style for root
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1);
  // --- start: initializing parameters and setting default values ---
  // detector:
  double edgeLength_ = 1100;
  std::vector<double> detectorDimensions_;
  detectorDimensions_.push_back(1100);
  detectorDimensions_.push_back(1100);
  detectorDimensions_.push_back(700);
  int numberOfSensors_ = 27;
  std::string output_ = "result";
  
  // distributions:
  int seed_ = 0;
  double energyExponent_ = 2.7; // negative exponent of the power function
  double energyMin_ = 1e1; // Minimum of the generated energy spectrum in GeV
  double energyMax_ = 1e5; // Maximum of the generated energy spectrum in GeV
  // factor of pi for azimuthal spectrum of starting positions
  double azimuthalSpectrum_ = 2; 
  // factor of pi for polar spectrum of starting positions
  double polarSpectrum_ = 0;

  int amountOfSamples_ = 1e4;
  int measuredSamples_ = 0;

  TDetector *detector_; // addresses the used detector
  std::vector<TEvent*> samples_;  // can store all created events
  // --- end: initializing parameters and setting default values ---

  // --- start: assigning and loading values out of config ---
  po::options_description desc_("parameters");
  desc_.add_options()
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
    ("origin.azimuthal", po::value<double>(&azimuthalSpectrum_),
      "factor of pi for azimuthal angle")
    ("origin.polar", po::value<double>(&polarSpectrum_),
      "factor of pi for polar angle")
    ("settings.output", po::value<std::string>(&output_), "output filename");

  // create map that will store the loaded values
  po::variables_map vm_;
  // load the assigned parameters of the config file into the map vm_
  LoadConfig(desc_, vm_);
  if ((azimuthalSpectrum_ < 0) || (azimuthalSpectrum_ > 2)) {
    azimuthalSpectrum_ = 2;
    std::cout << "Origin.azimuthal has to be in range [0,2]." << std::endl;
    std::cout << "Set it back to 2." << std::endl;
  }
  if ((polarSpectrum_ <= 0) || (polarSpectrum_ > 1)) {
    polarSpectrum_ = 1;
    std::cout << "Origin.polar has to be in range [0,1]." << std::endl;
    std::cout << "Set it back to 0." << std::endl;
  }
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
  // --- end: creating and ploting detector ---

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
  data->Branch("CoI_XPosition", &centerOfIntensityX, "CoI_XPosition/D");
  data->Branch("CoI_YPosition", &centerOfIntensityY, "CoI_YPosition/D");
  data->Branch("CoI_ZPosition", &centerOfIntensityZ, "CoI_ZPosition/D");
  data->Branch("Event_ID", &eventId, "Event_ID/Is");
  data->Branch("Event_Seed", &seed_, "Event_Seed/I");
  data->Branch("Effective_Track_Length", &effectiveTrackLength,
               "Effective_Track_Length/D");
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
  // --- end: initializing parameters for detection ---

  // --- start: measurement for all samples ---
  std::cout << "Starting measurement." << std::endl;
  for (int i = 0; i < amountOfSamples_; i++) {
    // reset event features
    bool measured = false;
    double minArrivalTime = 100, maxArrivalTime = 0;
    TSensor *firstSensor, *lastSensor;

    activatedSensors = 0;
    activatedSensorsPerGroup.resize(0);
    totalIntensity = 0;

    // increase event id
    eventId++;

    // --- start: generating distributed random numbers ---
    energy = GetNewEnergy(energyMin_, energyMax_, energyExponent_, grandom);

    std::vector<double> newPosition;
    newPosition = GetNewTrack(detectorDimensions_, azimuthalSpectrum_,
                              polarSpectrum_, rng);
    // check if GetNewTrack could produce a track that hit the detector
    // if not, request a new one
    while (newPosition.at(0) < 0) {
      newPosition = GetNewTrack(detectorDimensions_, azimuthalSpectrum_,
                                polarSpectrum_, rng);
    }
    azimuthalAngle = newPosition.at(3);
    polarAngle = newPosition.at(4);
    for (unsigned int j = 0; j < 3; j++) { position.at(j) = newPosition.at(j); }
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
      // collect the total deposited intensity
      totalIntensity += intensity;
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
    // ToDo: Abfangen falls keine Messung
    if (measured) {
      measuredSamples_ += 1;
      firstPosition = reconstructedTrack.GetPerpFoot(firstSensor->GetPosition());
      lastPosition = reconstructedTrack.GetPerpFoot(lastSensor->GetPosition());
      effectiveTrackLength = reconstructedTrack.GetDistance(firstPosition,
                                                            lastPosition); 
    }
    // --- end: generate features of reconstructed track ---

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

  return 0;
}

double DegToRad(double degree) {
  // load pi from boosts constants with double precision
  double pi = boost::math::constants::pi<double>();

  return (degree / 180) * pi;
}

std::vector<double> GetCrossPlane(std::vector<double> centralPoint,
                                  std::vector<double> edgeLength,
                                  std::vector<double> start,
                                  std::vector<double> direction) {
  std::vector<double> result;
  result.resize(4);
  // gradient of the track at that the intersection takes place
  double gradient;
  // pre assuem that there is no intersection
  bool intersection = false;
  //std::cout << "in crossplane" << std::endl;
  for (int i = 0; i < 3; i++) {
    // store current edgeLength
    double length = edgeLength.at(i);
    // if this one is zero calculate the distance of the starting point to
    // this plane
    if (length == 0) {
      // distance in coordinate of the plane (centralPoint, start)
      double dLinePlane = centralPoint.at(i) - start.at(i);
      // calculate the gradient, but only if the denominator is not 0
      if (direction.at(i) != 0) {
        gradient = dLinePlane / direction.at(i);
        if (gradient != 0) {
          intersection = true;
          for (int j = 1; j < 4; j++) {
            // calculate intersection point
            result.at(j) = start.at(j-1) + direction.at(j-1) * gradient;
          }
          result.at(0) = sqrt(dLinePlane * dLinePlane);
        }
      }
    }
  }
  // check if the intersection point is within the borders of a plane
  for (int i = 0; i < 3; i++) {
    if ((result.at(i+1) > (centralPoint.at(i) + edgeLength.at(i) * 0.5))
      || (result.at(i+1) < (centralPoint.at(i) - edgeLength.at(i) * 0.5))) {
      intersection = false;
    }
  }
  // if there is no intersection set the result to -1's 
  if (!intersection) {
    for (int j = 0; j < 4; j++) {
      result.at(j) = -1;
    }
  }

  return result;
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
                                double azimuthalSpectrum,
                                double polarSpectrum,
                                boost::mt19937& rng) {
  // --- start: predefining distributions ---
  double pi = boost::math::constants::pi<double>();

  boost::uniform_01<> uniDist;
  boost::variate_generator<boost::mt19937&, boost::uniform_01<> >
  GetUni(rng, uniDist);

  boost::uniform_real<> azDist(0, azimuthalSpectrum * pi);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
  GetAzAngle(rng, azDist);

  boost::uniform_real<> poDist(0, polarSpectrum * pi);
  boost::variate_generator<boost::mt19937&, boost::uniform_real<> >
  GetPoAngle(rng, poDist);
  // --- end: predefining distributions ----

  double az, norm = 0, po, radius = 0, randomNumber;
  std::vector<double> direction, inDetector, intersection, normal, start, temp,
                      result;
  std::vector<std::vector<double> > intersections;
  direction.resize(3);
  inDetector.resize(3);
  intersection.resize(3);
  normal.resize(3);
  start.resize(3);
  temp.resize(3);
  result.resize(5);

  for (int i = 0; i < 3; i++) {
    randomNumber = GetUni();
    inDetector.at(i) = randomNumber * detectorDimensions.at(i);
    radius += detectorDimensions.at(i) * detectorDimensions.at(i);
    start.at(i) = detectorDimensions.at(i) / 2;
  }
  //dumpVector("inDetector", inDetector);
  //dumpVector("start", start);
  radius = 0.5 * sqrt(radius);
  az = GetAzAngle();
  po = GetPoAngle();
  //std::cout << "Radius " << radius << " az " << az << " po " << po << std::endl;
  temp = SphericalToCartesian(radius, az, po);
  //dumpVector("temp", temp);
  for (int i = 0; i < 3; i++) {
    start.at(i) += temp.at(i);
    direction.at(i) = inDetector.at(i) - start.at(i);
    norm += direction.at(i) * direction.at(i);
  }
  //dumpVector("new start", start);
  //dumpVector("direction", direction);
  norm = sqrt(norm);
  // --- start and inDetector are ready
  // calculate azimuthal and polar angle from cartesian direction vector
  result.at(3) = atan2(direction.at(1), direction.at(0));
  result.at(4) = acos(direction.at(2) / norm);
  // preset starting position to -1 incase there is no intersection
  for (int j = 0; j < 3; j++) { result.at(j) = -1; }

  // --- start: finding intersection point between track and detector
  // preparing vectors with the central points and dimensions of each plane
  std::vector<std::vector<double> > centralPoints, edgeLengths;
  for (int i = 0; i < 3; i++) {
    std::vector<double> centralPoint1, centralPoint2, edgeLength;
    for (int j = 0; j < 3; j++) {
      if (i == j) {
        edgeLength.push_back(0);
        centralPoint1.push_back(0);
        centralPoint2.push_back(detectorDimensions.at(j));
      }
      else {
        edgeLength.push_back(detectorDimensions.at(j));
        centralPoint1.push_back(detectorDimensions.at(j) / 2);
        centralPoint2.push_back(detectorDimensions.at(j) / 2);
      }
    }
    edgeLengths.push_back(edgeLength);
    edgeLengths.push_back(edgeLength);
    centralPoints.push_back(centralPoint1);
    centralPoints.push_back(centralPoint2);
  }

  // get intersection points and corresponding distances for each plane
  for (int i = 0; i < 6; i++) {
    intersection.resize(4);
    intersection = GetCrossPlane(centralPoints.at(i), edgeLengths.at(i),
                                 start, direction);
    // only store information if there is an interection within the borders of
    // the addressed plane
    if (intersection.at(0) >= 0) {
      intersections.push_back(intersection);
    }
  }

  // due to a bug there are some situations with only 1 intersection per track
  // to avoid the influence of them, only tracks with 2 intersections are taken
  if (intersections.size() == 2) {
    int nearestIntersection = 0;
    // deside wich intersection point is the clothest to the starting point
    if (intersections.at(1).at(0) < intersections.at(0).at(0)) {
      nearestIntersection = 1;
    }
    for (int i = 0; i < 3; i++) {
      // add intersection point to result
      result.at(i) = intersections.at(nearestIntersection).at(i + 1);
    }
  }

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

std::vector<double> SphericalToCartesian(double r, double azimuthal,
                                         double polar) {
  std::vector<double> result;
  result.push_back(r * cos(azimuthal) * sin(polar));
  result.push_back(r * sin(azimuthal) * sin(polar));
  result.push_back(r * cos(polar));

  return result;
}
