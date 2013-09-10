#ifndef TDETECTOR_H
#define TDETECTOR_H
#include <fstream>
#include <iostream>
#include <math.h>
#include <vector>
#include "TEvent.h"
#include "TSensor.h"

/**
 * The detector class
 * A detector is a 3 dimensional cubic that contains sensors.
 * The coordinates used for this detector are cartesian. Their origin is located
 * in one corner of the cubical-shaped detector. Sensors are placed only with
 * positive coordinates. Directions are implemented via azimuthal and polar
 * angles. The azimuthal angle (-π to π) is the angle between the directional
 * vector and the x-axis counter-clockwise. The polar angle is used for
 * rotations in relation to the z-axis (0 to π). That means if the azimuthal and
 * polar angle are zero the direction is (0, 0, 1).
 */

class TDetector{
  public:
    /**
     * @brief TDetector creates a detector with no sensors.
     */
    TDetector();

    /**
     * @brief TDetector creates a detector with a given size and a given amount
     * of sensors that belong to it
     * @param edgeLength defines the length of the edges of the squared
     * detector
     * @param numOfSensors defines the amount of sensors to be created
     * This constructor creates a detector and a given amount of sensors.
     * Those sensors are equaly distributed in space and have the same 
     * permeability, that is set to a default 1 (that equals 100%) and belong to
     * sensor group 1. 
     * Have a look at PlaceSensors for more information about sensor placement.
     * The sensor structure is written in a file ("detector_config.txt").
     */
    TDetector(double edgeLength, int numberOfSensors);

    /**
     * @brief TDetector creates a detector with a given size and a given amount
     * of sensors that belong to it
     * @param edgeLength defines the length of the edges of the squared
     * detector
     * @param numOfSensors defines the amount of sensors to be created
     * @param permeability is the measurement probability of a sensor
     * This constructor creates a detector and a given amount of sensors.
     * Those sensors are equaly distributed in space and have the same
     * permeability and belong to sensor group 1.
     * Have a look at PlaceSensors for more information about sensor placement.
     * The sensor structure is written in a file ("detector_config.txt").
     */
    TDetector(double edgeLength, int numberOfSensors,
              double permeability);

    /**
     * @brief TDetector creates a detector with a given size and predefined
     * sensors
     * @param filename contains the filename of the configuration file that
     * contains information about sensors (position and permeability)
     * This constructors creates a detector with sensors, that have a
     * predefined position and accectance.
     * The information for those sensors is provided by a configuration file.
     * This file needs to be a plain txt-file that awaits four columns.
     * The first three are for the position (x, y, z), the fourth is for the
     * permeability and the fith is for the id of the sensor group. So one row
     * represents one sensor.
     * Keep in mind that the edgeLength is not calculated, and is set to zero.
     */
    TDetector(std::string filename);

    /**
     * @brief AddSensor adds a referenced, not already existing sensor
     * @param newSensor is a reference to a sensor that should be added
     * Keep in mind that the edgeLength is not calculated, and is set to zero.
     */
    void AddSensor(TSensor &newSensor);

    /**
     * @brief this method calculates and stores measured intensities for each
     * sensor
     * @param event is the event that is should be detected
     * @param randomNumber is a random number between 0 and 1 that is used to
     * decide weather the event has been detected or not
     */
    void Detect(TEvent &event, double randomNumber);

    /**
     * @breif this method calculates measured intensities for each sensor and
     * adds the measured value to the existing intensity values of the sensors
     * @param event is the event that should be added to the measurement
     * @param randomNumber is a random number between 0 and 1 that is used to
     * decide weather the event has been detected or not
     */
    void DetectAlso(TEvent &event, double randomNumber);

    /**
     * @brief Dump prints all the values of all sensors and the edge length
     */
    void Dump();

    /**
     * @brief GetCenterOfIntensity calculates the center of intensity
     * @return position of the center of intensity
     * It is like the center of mass but weighted with intensities instead.
     */
    std::vector<double> GetCenterOfIntensity();

    /**
     * @brief GetEdgeLength returns the edge length of the detector
     * @return edge length of the cubical-shaped detector
     */
    double GetEdgeLength();

    /**
     * @brief GetGroup returns the adress of a vector that stores the group "id"
     * @param id specifies the group that should be returned
     * @return adress of the vector with all sensors of the groud "id" 
     */
    const std::vector<TSensor*>& GetGroup(unsigned int id);

    /**
     * @brief GetNumberOfGroups returns the number of groups of sensors
     * @return number of groups of sensors
     */
    unsigned int GetNumberOfGroups();

    /**
     * @brief GetNumberOfSensors returns the number of created sensors
     * @return number of sensors
     */
    unsigned int GetNumberOfSensors();

    /**
     * @brief GetSensors returns the adress of a vector that stores all sensors
     * with their values (eg. measured intensity)
     * @return adress of the vector with all the sensor references
     */
    const std::vector<TSensor*>& GetSensors() const;

    /**
     * @brief ResetSensors sets the intensity values of all sensors to zero
     */
    void ResetSensors();

    /**
     * @brief ReconstructTrack uses a simple approach for a 3 dimensional linear
     * regression. It calculates the offset and gradient for two 2 dimension
     * planes usgin the analytic solution for a 2 dimensional linear regression.
     * Those planes are (x,y) and (y,z).
     * @return a TEvent object, that hols all the needed information about the track
     * the energy value is default (0).
     */
    const TEvent& ReconstructTrack();

    /**
     * @brief SetEdgeLength resets the edge length to a given value
     * @param edgeLength contains the new edgeLength value
     */
    void SetEdgeLength(double edgeLength);

    /**
     * @brief ~TDestroy destroys the detector object.
     */
    ~TDetector();
  protected:
    // edgeLength of the cubical-shaped detector in metre
    double edgeLength_;
    // groups of sensors that belong together (based on their group id)
    std::vector<std::vector<TSensor*> > groups_;
    // sensors belonging to the detector that are used to detect
    std::vector<TSensor*> sensors_;

    /** 
     * @brief GetAcceptance calculates the acceptance for a given energy
     * @param energy is the energy value for that an acceptance is calculated
     * @return acceptance probability for the given energy
     */
    double GetAcceptance(double energy);

    /**
     * @brief GrouSensors sorts sensors based on their group id into groups_
     */
    void GroupSensors();

    /**
     * @brief PlaceSensors creates sensors with equidistant positions, given
     * permeability values and sensor group id 1.
     * @param numberOfSensors hold the amount of sensors that should be created
     * @param permeability holds the detection propability of a sensor
     * To generate sensor positions this method takes the 3rd root of
     * numberOfSensors as an amount of sensors per dimension. If this number is
     * not natural the next natural number is chosen. The origin of the
     * positions is the origin of the detector (0, 0, 0). Those positions are
     * filled layer by layer with the desired numberOfSensors. The distance is
     * calculated in a way that all possible positions have the same distance.
     * Eg. numberOfSensors = 5, edgeLength = 3:
     * (1, 1, 1), (1, 1, 2), (1, 2, 1), (1, 2, 2), (2, 1, 1)
     */
    void PlaceSensors(int numberOfSensors, double permeability);

    /**
     * @brief PositionTaken tests if the given position is already taken
     * @param testPosition is the position that is tested
     * @return true if the position is already taken, otherwise false
     */
    bool PositionTaken(std::vector<double> testPosition);

    /**
     * @brief SetEdgeLength sets edgeLength to the highest x, y, or z position
     */
    void SetEdgeLength();

    /**
     * @brief WriteToFile generates a file with the positions and permeability of
     * all sensors
     */
    void WriteToFile();
};

#endif // TDETECTOR_H
