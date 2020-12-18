/*
 * simulator.C
 *
 *  Created on: September 5, 2018
 *      Author: Sean Flynn
 */

#include "simulator.h"

using namespace std;

//constants
const int LARGE_NUMBER = 1000000000;
const int LARGE_NEG_NUMBER = LARGE_NUMBER * -1;
const double EPSILON = .00001;
const double FLUID_DENSITY = 1.0;

//Constructs a Simulator object.
Simulator::Simulator(MacGrid *grid, const UT_Vector3 gravity, const double st_const,
                     const double flip_ratio) : 
    _gravity_(gravity), _st_const_(st_const), _flip_ratio_(flip_ratio)
{
    _grids_.push_back(grid);
    _particle_positions_.emplace_back();
    _particle_velocities_.emplace_back();

    UT_Vector3 init_vel1(fRand(-2.0, 2.0), fRand(-2.0, 2.0), 0.0);
    fill_box_with_particles(12, 18, 22, 28, 0, 0, init_vel1, 0);
    UT_Vector3 init_vel2(fRand(-2.0, 2.0), fRand(-2.0, 2.0), 0.0);
    fill_box_with_particles(46, 52, 22, 28, 0, 0, init_vel2, 0);
}

//Deletes this Simulator and any associated data.
Simulator::~Simulator()
{
    for(int i = 0; i < _grids_.size(); i++)
        delete _grids_[i];
}

//Simulates and caches the FLIP fluid up until the specified frame or returns if the simulation
//is already cached to that frame.
void Simulator::simulate_flip_to_frame(const size_t frame)
{
    if(frame < _grids_.size())
        return;

    size_t start_frame = _grids_.size() - 1;
    size_t end_frame = frame;
    double cur_time = start_frame;
    double max_u, ts;
    MacGrid *grid;
    vector<UT_Vector3> particle_positions;
    vector<UT_Vector3> particle_velocities;
    
    for(size_t f = start_frame; f < end_frame; f++)
    {
        grid = new MacGrid(*(_grids_[f]));
        particle_positions = _particle_positions_[f];
        particle_velocities = _particle_velocities_[f];

        while(cur_time < f + 1)
        {
            max_u = grid->max_u();
            ts = grid->voxel_size() / max_u;
            ts = min(f + 1 - cur_time, ts);

            grid->advect_particles(particle_positions, ts);
            grid->enforce_particle_bounds(particle_positions, particle_velocities);
            grid->compute_sdf(particle_positions);
            grid->update_buffer(1, particle_positions);
            apply_gravity_to_particles(particle_velocities, ts);
            grid->pvel_to_grid(particle_positions, particle_velocities);
            grid->set_boundary_velocities();
            grid->extrapolate_velocity(4);
            grid->set_boundary_velocities();
            grid->pressure_projection(FLUID_DENSITY, _st_const_, ts);
            grid->extrapolate_velocity(4);
            grid->set_boundary_velocities();
            grid->grid_to_pvel(particle_positions, particle_velocities, _flip_ratio_);
            
            cur_time += ts;
        }
        _grids_.push_back(grid);
        _particle_positions_.push_back(particle_positions);
        _particle_velocities_.push_back(particle_velocities);
    }
}

//Adds the gravity constant to the provided particle velocities.
void Simulator::apply_gravity_to_particles(vector<UT_Vector3>& particle_velocities,
                                           const double t) const
{
    for(size_t p = 0; p < particle_velocities.size(); p++)
    {
        for(size_t dim = 0; dim < 3; dim++)
            particle_velocities[p][dim] += _gravity_[dim] * t;
    }
}

//Adds a particle to the specified simulation frame with the provided position and velocity.
void Simulator::add_particle(const double px, const double py, const double pz,
                             const double ux, const double uy, const double uz, const size_t frame)
{
    _particle_positions_[frame].push_back(UT_Vector3(px, py, pz));
    _particle_velocities_[frame].push_back(UT_Vector3(ux, uy, uz));
}

//Creates a box of particles at the specified simulation frame with the provided position,
//dimensions, and initial velocity.
void Simulator::fill_box_with_particles(const size_t min_x, const size_t max_x,
                                        const size_t min_y, const size_t max_y,
                                        const double rand_min, const double rand_max,
                                        const UT_Vector3& init_vel, const size_t frame)
{
    //TODO--make this dependent on the voxel size rather than hard coded to a voxel size of 1
    double offset =  .25;
    double x_pos, y_pos;

    for(size_t j = min_y; j < max_y; j++)
    {
        for(size_t i = min_x; i < max_x; i++)
        {
            x_pos = (i + offset + fRand(rand_min, rand_max));
            y_pos = (j + offset + fRand(rand_min, rand_max));
            add_particle(x_pos, y_pos, 0, init_vel[0], init_vel[1], init_vel[2], frame);

            x_pos = (i + .5 + offset + fRand(rand_min, rand_max));
            y_pos = (j + offset + fRand(rand_min, rand_max));
            add_particle(x_pos, y_pos, 0, init_vel[0], init_vel[1], init_vel[2], frame);

            x_pos = (i + offset + fRand(rand_min, rand_max));
            y_pos = (j + .5 + offset + fRand(rand_min, rand_max));
            add_particle(x_pos, y_pos, 0, init_vel[0], init_vel[1], init_vel[2], frame);

            x_pos = (i + .5 + offset + fRand(rand_min, rand_max));
            y_pos = (j + .5 + offset + fRand(rand_min, rand_max));
            add_particle(x_pos, y_pos, 0, init_vel[0], init_vel[1], init_vel[2], frame);
        }
    }
}
