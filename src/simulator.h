/*
 * simulator.h
 *
 *  Created on: September 5, 2018
 *      Author: Sean Flynn
 */

#ifndef SIMULATOR_H_
#define SIMULATOR_H_

#include <vector>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <fstream>

#include "macgrid.h"

/**
* This class defines a FLIP fluid simulator. It caches each frame of simulation as it simulates
* similar to other Houdini simulators. This means that there is a vector of particles and MacGrid
* object pointers that coorespond to each frame that has been simulated. It contains methods for
* the main FLIP algorithm as well as methods for changing visible simulation parameters and for
* adding particles to the simulation. Most of the detailed simulation logic is found in the MacGrid
* class while this class contains the high level algorithm.
*/
class Simulator {
public:

    /**
    * Constructs a Simulator object.
    *
    * @param grid The MacGrid where the simulation takes place.
    * @param gravity The gravity constant.
    * @param st_const The surface tension constant.
    * @param flip_ration The flip to pic ratio.
    */
    Simulator(MacGrid *grid, UT_Vector3 gravity, double st_const, double flip_ratio);

    /**
    * Deletes this Simulator and any associated data.
    */
    virtual ~Simulator();

    /**
    * Getter methods
    */
    const UT_Vector3& get_gravity(void) { return _gravity_; }
    double get_st_const(void) { return _st_const_; }
    double get_flip_ratio(void) { return _flip_ratio_; }
    const MacGrid* get_grid(void) { return _grids_[_grids_.size() - 1]; }
    const MacGrid* get_grid(size_t frame)
        { simulate_flip_to_frame(frame); return _grids_[frame]; }
    const std::vector<std::vector<UT_Vector3>>& get_particle_positions(void) { return _particle_positions_; }
    const std::vector<std::vector<UT_Vector3>>& get_particle_velocities(void) { return _particle_velocities_; }
    const std::vector<UT_Vector3>& get_particle_positions(size_t frame)
        { return _particle_positions_[frame]; }
    const std::vector<UT_Vector3>& get_particle_velocities(size_t frame)
        { return _particle_velocities_[frame]; }

    /**
    * Setter methods
    */
    void set_gravity(UT_Vector3& gravity)
        { _gravity_[0] = gravity[0]; _gravity_[1] = gravity[1]; _gravity_[2] = gravity[2]; }
    void set_st_const(double st_const) { _st_const_ = st_const; }
    void set_flip_ratio(double flip_ratio) { _flip_ratio_ = flip_ratio; }

    /**
    * Simulates and caches the FLIP fluid up until the specified frame or returns if the simulation
    * is already cached to that frame.
    *
    * @param frame The frame to simulate up until.
    */
    void simulate_flip_to_frame(size_t frame);

    /**
    * Adds the gravity constant to the provided particle velocities.
    *
    * @param[out] particle_velocities The particle velocities to modify by adding gravity.
    * @param t The timestep.
    */
    void apply_gravity_to_particles(std::vector<UT_Vector3>& particle_velocities, double t) const;

    /**
    * Adds a particle to the specified simulation frame with the provided position and velocity.
    *
    * @param px The x position of the particle in world units.
    * @param py The y position of the particle in world units.
    * @param pz The z position of the particle in world units. NOTE--this implementation is 2D,
    *   so the z position will be ignored.
    * @param ux The x component of velocity of the particle.
    * @param uy The y component of velocity of the particle.
    * @param uz The z component of velocity of the particle. NOTE--this implementation is 2D,
    *   so the z velocity will be ignored.
    * @param frame The simulation frame to add the particle to.
    */
    void add_particle(double px, double py, double pz,
                      double ux, double uy, double uz, size_t frame);

    /**
    * Creates a box particles at the specified simulation frame with the provided position,
    * dimensions, and initial velocity.
    *
    * @param min_x The minimum x bound of the box.
    * @param max_x The maximum x bound of the box.
    * @param min_y The minimum y bound of the box.
    * @param max_y The maximum y bound of the box.
    * @param rand_min The minimum amount to randomly perturb the particle locations.
    * @param rand_max The maximum amount to randomly perturb the particle locations.
    * @param init_vel The initial velocity of the particles in the box.
    * @param frame The simulation frame to add the particles to.
    */
    void fill_box_with_particles(size_t min_x, size_t max_x,
                                 size_t min_y, size_t max_y,
                                 double rand_min, double rand_max,
                                 UT_Vector3& init_vel, size_t frame);

private:

    //simulation constants
    UT_Vector3 _gravity_;
    double _st_const_;
    double _flip_ratio_;

    //vector of cached simulation grids
    std::vector<MacGrid*> _grids_;

    //vector of cached particles
    std::vector<std::vector<UT_Vector3>> _particle_positions_;
    std::vector<std::vector<UT_Vector3>> _particle_velocities_;
};

#endif
