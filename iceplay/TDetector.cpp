/**
 *  Author: Philipp Schlunder
 *  Licence: GPL v3
 */

#include "TDetector.h"

TDetector::TDetector() : edgeLength_(0) {
  sensors_.resize(0);
  groups_.resize(0);
  // put all sensors in one group since all are created with id = 1
  groups_.push_back(sensors_);
}

TDetector::TDetector(double edgeLength, int numberOfSensors)
  : edgeLength_(edgeLength) {
  sensors_.resize(0);
  groups_.resize(0);
  PlaceSensors(numberOfSensors, 1.0);
  // save sensor position into detector_config.txt
  WriteToFile();
  // put all sensors in one group since all are created with id = 1
  groups_.push_back(sensors_);
}

TDetector::TDetector(double edgeLength, int numberOfSensors,
                     double permeability)
  : edgeLength_(edgeLength) { 
  PlaceSensors(numberOfSensors, permeability);
  // save sensor position into detector_config.txt
  WriteToFile();
  // put all sensors in one group since all are created with id = 1
  groups_.push_back(sensors_);
}

TDetector::TDetector(std::string filename) : edgeLength_(0) {
  sensors_.resize(0);
  double x, y, z, permeability;
  unsigned int group;
  // open file with given filename for reading
  std::ifstream file (filename.c_str());
  // read the first line with commentary
  std::string commentLine;
  std::getline(file, commentLine);
  // read all values and put them together into a position
  while (file.good()) {
    file >> x >> y >> z >> permeability >> group;
    std::vector<double> tempPosition;
    tempPosition.push_back(x);
    tempPosition.push_back(y);
    tempPosition.push_back(z);
    // create a new sensor out of the given values and add it to the detector
    AddSensor(*(new TSensor(tempPosition, permeability, group)));
  }
  // sort sensors into groups
  GroupSensors();
  // set edgeLength_ to the highest used coordinate
  SetEdgeLength();
}

void TDetector::AddSensor(TSensor &newSensor) {
  // check wether the position of the newSensor is already taken
  if (!PositionTaken(newSensor.GetPosition())) {
    // if not, add the new sensor
    sensors_.push_back(&newSensor);
  }
}

void TDetector::Detect(TEvent &event, double randomNumber) {
  double speedOfLight = 3e8;  // in m / s
  // measure for each sensor
  for (unsigned int i = 0; i < sensors_.size(); i++) {
    std::vector<double> sensorPosition = sensors_.at(i)->GetPosition();
    std::vector<double> perpFoot = event.GetPerpFoot(sensorPosition);
    // use energy the event has at its nearest position to the sensor
    double energy = event.GetEnergy(perpFoot);
    // trigger sensor based on energy dependent acceptance
    if (GetAcceptance(energy) > randomNumber) {
      double permeability = sensors_.at(i)->GetPermeability();
      double sensorDistance = event.GetDistance(perpFoot, sensorPosition);
      double trackLength = event.GetDistance(event.GetPosition() ,perpFoot);
      // arrival time is the time the event needs to reach the nearest position
      // to the sensor + time the photons need from this position to the sensor
      sensors_.at(i)->SetArrivalTime((sensorDistance + trackLength)
                                     / speedOfLight);
      // calculate deposited energy and store it in the adressed sensor
      sensors_.at(i)->SetIntensity(energy * permeability
                                   / (1 + pow(sensorDistance, 2)));
    }
    // if the sensor did not trigger, make sure intensity is set to 0
    else {
      sensors_.at(i)->SetIntensity(0);
    }
  }
}

void TDetector::DetectAlso(TEvent &event, double randomNumber) {
  double speedOfLight = 3e8;  // in m / s
  // measure for each sensor
  for (unsigned int i = 0; i < sensors_.size(); i++) {
    std::vector<double> sensorPosition = sensors_.at(i)->GetPosition();
    std::vector<double> perpFoot = event.GetPerpFoot(sensorPosition);
    // use energy the event has at its nearest position to the sensor
    double energy = event.GetEnergy(perpFoot);
    // trigger sensor based on energy dependent acceptance
    if (GetAcceptance(energy) > randomNumber) {
      double permeability = sensors_.at(i)->GetPermeability();
      double sensorDistance = event.GetDistance(perpFoot, sensorPosition);
      double trackLength = event.GetDistance(event.GetPosition() ,perpFoot);
      // arrival time is the time the event needs to reach the nearest position
      // to the sensor + time the photons need from this position to the sensor
      sensors_.at(i)->SetArrivalTime((sensorDistance + trackLength)
                                     / speedOfLight);
      // calculate deposited energy and add it to the stored value
      sensors_.at(i)->SetIntensity(sensors_.at(i)->GetIntensity()
                                   + energy * permeability
                                   / (1 + pow(sensorDistance, 2)));
    }
  }
}

void TDetector::Dump() {
  std::cout << "EdgeLength: " << edgeLength_ << std::endl;
  for (unsigned int i = 0; i < sensors_.size(); i++) {
    std::cout << "Sensor " << i + 1 << ":" << std::endl;
    sensors_.at(i)->Dump();
  }
}

double TDetector::GetAcceptance(double energy) {
  // (1 - exp( -E / 2)) ^ 3
  return pow((1 - exp(- energy / 2)), 3);
}

std::vector<double> TDetector::GetCenterOfIntensity() {
  // need to calculate the center for each coordinate
  std::vector<double> center;
  // needing 3 coordiantes
  center.resize(3);

  // calculate center of intensity for each coordiante
  for (unsigned int pos = 0; pos < 3; pos++) {
    // sum of (with intensity) weighted positions
    double numerator = 0;
    // sum of weights (intensities)
    double denominator = 0;

    // calculate numerator and denominator
    for (unsigned int event = 0; event < sensors_.size(); event++) {
      double weight = sensors_.at(event)->GetIntensity();
      double position = sensors_.at(event)->GetPosition().at(pos);

      numerator += weight * position;
      denominator += weight;
    }
    // center of intensity = sum of weighted positions / sum of weights
    center.at(pos) = numerator / denominator;
  }

  return center;
}

double TDetector::GetEdgeLength() { return edgeLength_; }

const std::vector<TSensor*>& TDetector::GetGroup(unsigned int id) {
  if (id > 0) {
    // need to decrease the id because a vector starts with 0 and not with 1
    return groups_.at(id - 1);
  }
  else {
    std::cout << "Id has to exceed 0. Returning group with id 1." << std::endl;
    std::cout << "Keep in mind: The actual id is used, not its vector position."
              << std::endl;
    return groups_.at(0);
  }
}

unsigned int TDetector::GetNumberOfGroups() { return groups_.size(); }

unsigned int TDetector::GetNumberOfSensors() { return sensors_.size(); }

const std::vector<TSensor*>& TDetector::GetSensors() const {
  return sensors_;
}

void TDetector::GroupSensors() {
  unsigned int maxId = 0;

  // loop over sensors to find the max id
  for(unsigned int i = 0; i < sensors_.size(); i++) {
    if (sensors_.at(i)->GetGroup() > maxId) {
      maxId = sensors_.at(i)->GetGroup();
    }
  }

  // set group_ size to max id
  groups_.resize(maxId);

  // insert dummy vector with dummy sensor into groups_ to initialize it
  // need to push_back vectors of sensors to groups_ thats why one needs a dummy
  for (unsigned int i = 0; i < groups_.size(); i++) {
    std::vector<TSensor*> dummy;
    dummy.push_back(new TSensor());
    groups_.at(i) = dummy;
  }

  // actual sorting process 
  for (unsigned int i = 0; i < sensors_.size(); i++) {
    unsigned int id = sensors_.at(i)->GetGroup();
    // sort sensor into group with the given id
    // keep in mind a vector starts with 0 but id's start with 1
    groups_.at(id - 1).push_back(sensors_.at(i));
  }

  // erase dummy vectors and sensors
  for (unsigned int i = 0; i < groups_.size(); i++) {
    groups_.at(i).erase(groups_.at(i).begin());
  }
}

void TDetector::PlaceSensors(int numberOfSensors, double permeability) {
  sensors_.resize(0);
  // test if edgeLength exceeds 0
  if (edgeLength_ >= 1) {
    int LoopLimit = 0;  // limit for sensor coordinates
    double cubeRoot = pow(numberOfSensors,1./3.);
    LoopLimit = int(pow(numberOfSensors,1./3.));
    // There are LoopLimit^3 positions for sensors available, if numberOfSensors
    // exceeds this amount it is necessary to use (LoopLimit+1)^3 positions.
    if ((cubeRoot - LoopLimit) > 0) LoopLimit++;
    std::vector<double> tempPosition; // holds generated sensor positions
    for (int i = 0; i < 3; i++) tempPosition.push_back(0);  // default values
    double SensorCounter = 0; // counter for the amount of placed sensors
    // calculate distance for equidistant sensor placement
    double distance = (edgeLength_) / (LoopLimit + 1);
    // --- generate sensor positions ---
    for (int z = 1; z <= LoopLimit; z++) {
      for (int y = 1; y <= LoopLimit; y++) {
        for (int x = 1; x <= LoopLimit; x++) {
          // until enough sensors are created
          if (SensorCounter >= numberOfSensors) { break; }
          // each coordinate is dependend on the amount of existing sensors in
          // the given direction and the distance between each sensor
          tempPosition.at(0) = distance * x;
          tempPosition.at(1) = distance * y;
          tempPosition.at(2) = distance * z;
          std::vector<double> tempVector;
          for (int i = 0; i < 3; i++) {
            tempVector.push_back(tempPosition.at(i));
          }
          // add the new sensor to the detector
          sensors_.push_back(new TSensor(tempVector, permeability, 1));
          SensorCounter++;
        }
      }
    }
  }
  else {
    std::cout << "The edge length of the detector has to exceed 0."
              << std::endl;
  }
}

bool TDetector::PositionTaken(std::vector<double> testPosition) {
  // initial guess is that the position is not taken
  bool isTaken = false;
  // compare the coordinates of testPosition against all existing sensors
  for (unsigned int i = 0; i < sensors_.size(); i++) {
    if (testPosition.at(0) == sensors_.at(i)->GetPosition().at(0) &&
        testPosition.at(1) == sensors_.at(i)->GetPosition().at(1) &&
        testPosition.at(2) == sensors_.at(i)->GetPosition().at(2)) {
      sensors_.at(i)->Dump();
      // if there is a match, end the loop and change marker "isTaken" to true
      isTaken = true;
      // This might also be triggered if your detector_config.txt ends with an
      // empty line. If that's the case, just ignore the message.
      std::cout << "TDetector: The chosen position is already taken."
                << std::endl;
      break;
    }
  }
  return isTaken;
}

const TEvent& TDetector::ReconstructTrack() {
  double mean[3] = {0, 0, 0};  // storage for x, y and z position of the means
  double offset[2], gradient[2];  // storage for calculated offsets & gradients
  std::vector<double> centerOfIntensity = GetCenterOfIntensity();
  std::vector<double> positions[3];  // positions of activaed sensors

  // sum up and store coordinates of the activated sensors
  for (unsigned int i = 0; i < sensors_.size(); i++) {
    // check if sensors has been activated
    if (sensors_.at(i)->GetIntensity() > 0) {
      // for each coordinate
      for (unsigned int j = 0; j < 3; j++) {
        // get position
        double temp = sensors_.at(i)->GetPosition().at(j);
        // add it to the sum
        mean[j] += temp;
        // store position
        positions[j].push_back(temp);
      }
    }
  }

  // get number of activated sensors
  double activatedSensors = positions[0].size();
  // convert sum of coordinates to mean
  for (unsigned int i = 0; i < 3; i++) {
    // one does not simply divide by zero
    if (activatedSensors > 0) {
      mean[i] /= activatedSensors;
    }
  }

  // calculate offset and gradient of a 2 dimension line
  // first in the (x, y)-plane, then in the (y, z)-plane
  for (unsigned int i = 0; i < 2; i++) {
    double numerator = 0, denominator = 0;

    /**
     * using the analytic solution for a 2 dimensional linear regression
     * gradient = sum_i^n ((x_i - mean(x)) * (y_i - mean(y)))
     *          / sum_i^n ((x_i - mean(x)) ^ 2)
     * offset = mean(y) - gradient * mean(x)
     */ 
    for (unsigned int j = 0; j < activatedSensors; j++) {
      // calculate the numerator and denominator separatly
      numerator += (positions[0 + i].at(j) - mean[0 + i])
                 * (positions[1 + i].at(j) - mean[1 + i]);
      denominator += pow((positions[0 + i].at(j) - mean[0 + i]), 2.);
    }
    // and execute the division
    gradient[i] = (numerator) / (denominator);
    offset[i] = mean[1 + i] - gradient[i] * mean[0 + i];
  }

  /**
   * preparing a position that might be part of the track
   * y(x) = offset[0] + x * gradient[0]
   * z(y) = offset[1] + y * gradient[1]
   * => z(x) = offset[1] + (offset[0] + x * gradient[0]) * gradient[1]
   * x = 0:
   * => y = offset[0] 
   * => z = offset[1] + offset[0] * gradient[1]
   */
  std::vector<double> trackPosition;
  trackPosition.push_back(0);
  trackPosition.push_back(offset[0]);
  trackPosition.push_back(offset[0] * gradient[1] + offset[1]);

  /**
   * Store track information as TEvent object, because it provides all variables
   * needed, besides the additional energy parameter. This provides the huge 
   * advantage, that one can receive projections of positions (eg. from sensors)
   * to the reconstructed track via using TEvent.GetPerpFoot().
   */
  double azimuthalAngle, polarAngle;
  azimuthalAngle = atan(gradient[0]);
  polarAngle = (M_PI/2) - atan(gradient[1]);

  TEvent *reconstructedTrack = new TEvent(0,
                                          trackPosition,
                                          azimuthalAngle,
                                          polarAngle);

  return *reconstructedTrack;
}

void TDetector::ResetSensors() {
  for (unsigned int i = 0; i < sensors_.size(); i++) {
    sensors_.at(i)->SetIntensity(0);
  }
}

void TDetector::SetEdgeLength() {
  double newEdgeLength = 0;
  for (unsigned int sensor = 0; sensor < sensors_.size(); sensor++) {
    for (int position = 0; position < 3; position++) {
      if (sensors_.at(sensor)->GetPosition().at(position) > newEdgeLength) {
        newEdgeLength = sensors_.at(sensor)->GetPosition().at(position);
      }
    }
  }
  edgeLength_ = newEdgeLength;
}

void TDetector::SetEdgeLength(double edgeLength) {
  if (edgeLength > 0) {
    edgeLength_ = edgeLength;
  }
  else {
    std::cout << "Please use a positive edge length." << std::endl;
    std::cout << "No new edge length set." << std::endl;
  }
}

void TDetector::WriteToFile() {
  std::ofstream file("../detector_config.txt");
  file << "#x\t y\t z\t permeability\t group" << std::endl;
  for (unsigned int i = 0; i < sensors_.size(); i++) {
    file << sensors_.at(i)->GetPosition().at(0) << "\t"
         << sensors_.at(i)->GetPosition().at(1) << "\t"
         << sensors_.at(i)->GetPosition().at(2) << "\t"
         << sensors_.at(i)->GetPermeability() << "\t"
         << sensors_.at(i)->GetGroup() << std::endl;
  }
  file.close();
}

TDetector::~TDetector() { }
