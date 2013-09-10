#ifndef TSENSOR_H
#define TSENSOR_H
#include <iostream>
#include <vector>

/**
 * An instance of TSensor is meant to simulate a sensor for energy detection.
 * This sensor has a position (in cartesian coordiantes)m a permeability, that
 * describes the percentage of energy that might be detected and is therefore 
 * used as some kind of material constant. Furthermore it contains an intensity
 * that is used to store the measured energy and an arrival time that stores the
 * time it took to measure the event since it entered the detector. Every sensor
 * belongs to a group, so that one can define a bunch of sensors that has the
 * same permeabilities or is related in another way.
 */

class TSensor{
  public:
    /**
     * @brief TSensor creates an object with a 3 dimensional null position,
     * permeability and group set to 1 and arrivalTime set to 0.
     */
    TSensor();

    /**
     * @brief TSensor creates an object with a given position, permeability,
     * group set to 1 and arrivalTime set to 0. Missing dimensions of the
     * position are filled up with zeros.
     * @param position contains coordinates of the sensors position
     */
    TSensor(std::vector<double> position);

    /**
     * @brief TSensor creates an object with a given position, permeability and 
     * group. The arrivalTime is set to 0. Missing dimensions of the position
     * are filled up with zeros.
     * @param position contains coordinates of the sensors position
     * @param permeability contains the measurement probability of the sensor
     * @param groupd contains the groupd id for the new sensor
     */
     TSensor(std::vector<double> position, double permeability,
             unsigned int group);

    /**
     * @brief Dump prints the values of the protected attributes
     */
    void Dump();

    /**
     * @brief GetArrivalTime returns the duration from the moment the event
     * entered the detector and the moment where it is measured
     * @return event duration: entering detector - beeing measured
     */
    double GetArrivalTime();

    /**
     * @brief GetGroup returns the stored it of the adressed sensor
     * @return id of the sensor
     */
    unsigned int GetGroup();

    /**
     * @brief GetIntensity returns the stored measured intensity of the adressed
     * sensor
     * @return stored measured intensity of the adressed sensor
     */
    double GetIntensity();

    /**
     * @brief GetPermeability returns the measurement probability of the
     * adressed sensor
     * @return measurement probability of the adressed sensor
     */
    double GetPermeability();

    /**
     * @brief GetPosition returns the coordinates of the adressed sensor
     * @return position of the adressed sensor
     */
    //TPosition GetPosition();
    std::vector<double> GetPosition();

    /**
     * @brief SetArrivalTime sets the duration the event needs to be measured
     * @param arrivalTime is the duration the event needs in seconds. If a
     * negative value is submitted it is set to 0.
     */
    void SetArrivalTime(double arrivalTime);

    /**
     * @brief SetGroup changes the group id of the adressed sensor
     * @param group contains the new id of the adressed sensor
     */
    void SetGroup(unsigned int group);

    /**
     * @brief SetIntensity changes the intensity of the adressed sensor
     * @param intensity conaints the new measured intensity for the
     * adressed sensor. The stored intensity is set to 0 if a negative intensity
     * is submitted.
     */
    void SetIntensity(double intensity);

    /**
     * @brief SetPermeability changes the permeability of the adressed sensor
     * @param permeability conaints the new measurement probability for the
     * adressed sensor. If it is negative, the value is set to 0.
     */
    void SetPermeability(double permeability);

    /**
     * @brief ~TSensor destroys the sensor object.
     */
    ~TSensor();
  protected:
    std::vector<double> position_;  // cartesian coordinates in metre
    double arrivalTime_;  // duration:event enters detector - gets measured in s
    double permeability_;   // probability of a sensor to function correctly
    double intensity_;  // measured intensity in GeV
    unsigned int group_;  // id of the group the sensor belongs to
};

#endif // TSENSOR_H
