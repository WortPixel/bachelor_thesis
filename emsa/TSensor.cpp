/**
 *  Author: Philipp Schlunder
 *  Licence: GPL v3
 */

#include "TSensor.h"

TSensor::TSensor()
  : arrivalTime_(0),
    permeability_(1),
    intensity_(0),
    group_(1) {
  // make sure position_ is a 3d vector
  position_.resize(0);
  while(position_.size() < 3) position_.push_back(0);
}

TSensor::TSensor(std::vector<double> position)
  : position_(position),
    arrivalTime_(0),
    permeability_(1),
    intensity_(0),
    group_(1) { 
  // make sure position_ is a 3d vector
  while(position_.size() < 3) position_.push_back(0);
}

TSensor::TSensor(std::vector<double> position, double permeability,
                 unsigned int group)
  : position_(position),
    arrivalTime_(0),
    permeability_(permeability),
    intensity_(0),
    group_(group) {
  // make sure position_ is a 3d vector
  while(position_.size() < 3) position_.push_back(0);
  // make sure permeability doesn't exceed 1.
  if(permeability > 1) {
    permeability_ = 1;
    std::cout << "The permeability should not exceed 1." << std::endl;
    std::cout << "Resettet permeability to 1." << std::endl;
  }
}

void TSensor::Dump() {
  std::cout << "TSensor:" << std::endl;
  std::cout << "Position (x, y, z / m): ("
            << position_.at(0) << ", "
            << position_.at(1) << ", "
            << position_.at(2) << ")" << std::endl;
  std::cout << "Permeability: " << permeability_ << std::endl;
  std::cout << "Group: " << group_ << std::endl;
  std::cout << "ArrivalTime / s: " << arrivalTime_ << std::endl;
  std::cout << "Intensity / GeV / m^2: " << intensity_ << std::endl;
  std::cout << "---" << std::endl;
}

double TSensor::GetArrivalTime() { return arrivalTime_; }

unsigned int TSensor::GetGroup() { return group_; }

double TSensor::GetIntensity() { return intensity_; }

double TSensor::GetPermeability() { return permeability_; }

std::vector<double> TSensor::GetPosition() { return position_; }

void TSensor::SetArrivalTime(double arrivalTime) { 
  if (arrivalTime >= 0) {
    arrivalTime_ = arrivalTime; 
  }
  else {
    arrivalTime_ = 0;
    std::cout << "ArrivalTime should be positive." << std::endl;
    std::cout << "Resetting arrival time to 0." << std::endl;
  }
}

void TSensor::SetGroup(unsigned int group) { group_ = group; }

void TSensor::SetIntensity(double intensity) {
  // if intensity is negativ, the sensor is to far away to measure something
  // so set the sensor value to 0.
  if(intensity < 0) { intensity_ = 0; }
  else { intensity_ = intensity; }
}

void TSensor::SetPermeability(double permeability) { 
  // make sure permeability doesn't exceed 1.
  if(permeability > 1) {
    permeability_ = 1;
    std::cout << "The permeability should not exceed 1." << std::endl;
    std::cout << "Resettet permeability to 1." << std::endl;
  }
  else {
    permeability_ = permeability;
  }
}

TSensor::~TSensor() { }
