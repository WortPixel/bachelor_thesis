#include "TEvent.h"

TEvent::TEvent() : energy_(0), azimuthalAngle_(0), polarAngle_(0) {	
  // make sure position_ is a 3d vector
  position_.resize(0);
  while(position_.size() < 3) position_.push_back(0);
}

TEvent::TEvent(double energy, std::vector<double> position,
               double azimuthalAngle)
  : energy_(energy),
    position_(position),
    azimuthalAngle_(azimuthalAngle),
    polarAngle_(0) {
  // make sure position_ is a 3d vector
  while(position_.size() < 3) position_.push_back(0);
}

TEvent::TEvent(double energy, std::vector<double> position,
               double azimuthalAngle, double polarAngle)
  : energy_(energy),
    position_(position),
    azimuthalAngle_(azimuthalAngle),
    polarAngle_(polarAngle) {
  // make sure position_ is a 3d vector
  while(position_.size() < 3) position_.push_back(0);
}

void TEvent::Dump() {
  std::cout << "TEvent:" << std::endl;
  std::cout << "Starting energy / GeV: " << energy_ << std::endl;
  std::cout << "Azimuthal and polar angle /rad: " << azimuthalAngle_
            << ", " << polarAngle_ << std::endl;
  std::cout << "Starting position (x, y, z / m): ("
            << position_.at(0) << ", "
            << position_.at(1) << ", "
            << position_.at(2) << ")" << std::endl;
}

double TEvent::GetAzimuthalAngle() { return azimuthalAngle_; }

double TEvent::GetDistance(std::vector<double> firstPosition,
                           std::vector<double> secondPosition) {
  std::vector<double> difference;
  // create difference vector between firstPosition and secondPosition
  for (int i = 0; i < 3; i++) {
      difference.push_back(secondPosition.at(i) - firstPosition.at(i));
  }
  // return distance
  return sqrt(GetInnerProduct(difference, difference));
}

double TEvent::GetEnergy() { return energy_; }

double TEvent::GetEnergy(std::vector<double> referencePosition) {
  // travelled distance of the event
  double distance = GetDistance(position_, referencePosition);
  // --- constant energy loss ---
  double remainingEnergy, energyLoss;
  /**
   * approximated energy loss per metre from
   * "PROPOSAL for muon propagation" by Koehne et al.
   */
  energyLoss = 0.259 + 3.64e-4 * energy_;
  // using linear energy loss
  remainingEnergy = energy_ - energyLoss * distance;
  // ----------------------------

  return remainingEnergy;
}

double TEvent::GetPolarAngle() { return polarAngle_; }

std::vector<double> TEvent::GetPosition() { return position_; }

std::vector<double> TEvent::GetPerpFoot(std::vector<double> referencePosition) {
  /**
   * create dropped perpendicular foot (perpFoot) first
   * g: S = position_ + t * eventDirection
   * t = <referencePosition - position_, eventDirection>
   *     / <eventDirection, eventDirection>
   */
  std::vector<double> eventDirection, difference, perpFoot;
  // create directional vector out of angles
  eventDirection.push_back(cos(azimuthalAngle_) * sin(polarAngle_));
  eventDirection.push_back(sin(azimuthalAngle_) * sin(polarAngle_));
  eventDirection.push_back(cos(polarAngle_));
  // create difference vector between sensorPosition and this event
  for (int i = 0; i < 3; ++i) {
      difference.push_back(referencePosition.at(i) - position_.at(i));
  }
  // calculate gradient of straight line
  double t = GetInnerProduct(difference, eventDirection)
           / GetInnerProduct(eventDirection, eventDirection);
  // calculate dropped perpendicular foot using the straight line equation g
  for (int i = 0; i < 3; ++i) {
      perpFoot.push_back(position_.at(i) + t * eventDirection.at(i));
  }

  return perpFoot;
}

void TEvent::SetAzimuthalAngle(double azimuthalAngle) {
  azimuthalAngle_ = azimuthalAngle;
}

void TEvent::SetEnergy(double energy) { energy_ = energy; }

void TEvent::SetPolarAngle(double polarAngle) { polarAngle_ = polarAngle; }

void TEvent::SetPosition(std::vector<double> position) {
  for(unsigned int i = 0; i < 3; i++) { position_.at(i) = position.at(i); }
}

double TEvent::GetInnerProduct(const std::vector<double> &aFactor,
                               const std::vector<double> &bFactor) {
  double result = 0;
  // check weather both vectors have the same dimension
  if (aFactor.size() == bFactor.size()) {
    for (unsigned int i = 0; i < aFactor.size(); ++i) {
      // calculate the scalar product
      result += aFactor.at(i) * bFactor.at(i);
    }
  }
  else { std::cout << "TEvent::GetInnerProduct():"
                   << "Factors differ in dimension." << std::endl;
  }

  return result;
}

TEvent::~TEvent() { }
