/**
 *  Author: Philipp Schlunder
 *  Licence: GPL v3
 */

#ifndef TEVENT_H
#define TEVENT_H
#include <iostream>
#include <math.h>
#include <vector>

/**
 * An event simulates a particle that has a starting position (cartesian
 * coordinates), an energy, and a direction that is given by two angles.
 * Energy loss can be calculated to any given position. It is assumed that a
 * particle travels on a straight line and experiences a constant energy loss.
 * The calculated energy is the remaining energy the particle has left while it
 * passes the given position at its closest point.
 */

class TEvent{
public:
    /**
     * @brief TEvent creates a new event with a zero energy, zero angles and
     * a position that is placed at the origin of the coordinate (and detector)
     * system
     */
    TEvent();

    /**
     * @brief TEvent creates a new event with a starting energy, a given
     * position and direction (angle) for 2 dimensional calculations
     * @param energy contains the starting energy
     * @param position contains the coordiantes of the starting point
     * @param azimuthalAngle containts the aizumathal angle in relation to the
     * coordinates origin
     */
    TEvent(double energy, std::vector<double> position, double azimuthalAngle);

    /**
     * @brief TEvent creates a new event with a starting energy, a given
     * position and direction (angles) for 3 dimensional calculations 
     * @param energy contains the starting energy
     * @param position contains the coordiantes of the starting point
     * @param azimuthalAngle containts the aizumathal angle in relation to the
     * coordinates origin
     * @param polarAngle containts the polar angle in relation to the
     * coordinates origin
     */
    TEvent(double energy, std::vector<double> position, double azimuthalAngle,
           double polarAngle);

    /**
     * @brief Dump prints the values of the protected attributes
     */
    void Dump();

    /**
     * @brief GetAzimuthalAngle returns the azimuthal angle of the event
     * @return azimuthalAngle
     */
    double GetAzimuthalAngle();

    /**
     * @brief GetDistance returns the distance between firstPosition and
     * secondPosition
     * @param firstPosition is one of two needed position
     * @param secondPosition is the second position needed
     * @return the distance between firstPosition and secondPosition
     */
    double GetDistance(std::vector<double> firstPosition,
                       std::vector<double> secondPosition);

    /**
     * @brief GetEnergy returns the starting energie of the event
     * @return the initial energy of the event
     */
    double GetEnergy();

    /**
     * @brief GetEnergy calculates the remaining energy the event has at a given
     * position
     * @param referencePosition contains the coordinates for energy calculations
     * @return energy the event has left after propagating to the given position
     */
    double GetEnergy(std::vector<double> referencePosition);

    /**
     * @brief GetPerpFoot calculates the dropped perpendicular foot from a given
     * position to the linear tragectory of the event.
     * @param referencePosition is the position from which the perpendicular is
     * dropped
     * @return the position of the dropped perpendicular foot on the event line
     */
    std::vector<double> GetPerpFoot(std::vector<double> referencePosition);

    /**
     * @brief GetPolarAngle returns the polar angle of the event
     * @return polarAngle 
     */
    double GetPolarAngle();

    /**
     * @brief GetPosition returns the starting position of the event
     * @return cartesian coordinates of the starting position of the event
     */
    std::vector<double> GetPosition();

    /**
     * @brief SetAzimuthalAngle changes the azimuthal angle of the event
     * @param azimuthalAngle contains the new angle
     */
    void SetAzimuthalAngle(double azimuthalAngle);

    /**
     * @brief SetEnergy changes the starting energy of the event
     * @param energy contains the new energy value
     */
    void SetEnergy(double energy);

    /**
     * @brief SetPolarAngle changes the polar angle of the event
     * @param polarAngle contains the new angle
     */
    void SetPolarAngle(double polarAngle);

    /**
     * @brief SetPosition changes the starting position of the event
     * @param position contains the new starting position
     */
    void SetPosition(std::vector<double> position);

    /**
     * @brief ~TEvent destroys the event object
     */
    ~TEvent();
protected:
    double energy_; // energy in GeV
    std::vector<double> position_;  // cartesian coordinates in metre
    double azimuthalAngle_; // angle in xy-plane in radiant
    double polarAngle_; // angle in yz-plane in radiant

    /**
     * @brief GetInnerProduct calculates the inner product of the euclidean
     * space
     * @param aFactor is the first factor
     * @param bFactor is the second factor
     * @return the inner product of two vectors in the euclidean space
     */
    double GetInnerProduct(const std::vector<double> &aFactor,
                           const std::vector<double> &bFactor);
};

#endif // TEVENT_H
