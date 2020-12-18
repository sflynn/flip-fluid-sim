/*
 * macgrid.h
 *
 *  Created on: September 5, 2018
 *      Author: Sean Flynn
 */

#ifndef MACGRID_H_

#define MACGRID_H_

#include <vector>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include <eigen3/Eigen/Sparse>
#include <SOP/SOP_Node.h>

//The types of scalars that can be stored in a MacGrid, allows interpolation code reuse.
enum ScalarType {UX, UY, UZ, UX_TEMP, UY_TEMP, UZ_TEMP, U, SD, CURV, PRESSURE};

//The types of velocity fields in a MacGrid.
enum VelocityType {CUR, TEMP};

//The types of cells in a MacGrid. NOTE--the use of SOLID cells is not currently implemented.
enum CellType {FLUID, AIR, SOLID, UNUSED};

/**
* This class defines a uniform Marker and Cell (MAC) grid. It is used to store and compute data
* for fluid simulation. The axes match those of Houdini (positive x points to the right, positive
* y is up, and positive z is pointing out of the viewport).
*/
class MacGrid
{
public:

    /**
    * Constructs a MacGrid object.
    *
    * @param width The width in voxels of the grid.
    * @param height The height in voxels of the grid.
    * @param depth The depth in voxels of the grid. NOTE--this implementation is currently 2D,
        so this parameter is not used, but we provide it to allow easy extension to 3D.
    * @param voxel_size The width of an individual voxel in world units.
    * @param particle_radius The radius of FLIP particles in world units.
    */
    MacGrid(const size_t width, const size_t height, const size_t depth,
            const double voxel_size, const double particle_radius);

    /**
    * Constructs a MacGrid that is a copy of the input MacGrid.
    * 
    * @param grid The MacGrid to copy.
    */
    MacGrid(const MacGrid &grid);

    /**
    * Deletes this MacGrid and any associated data.
    */
    virtual ~MacGrid();

    /**
    * Getter methods
    */
    size_t width(void) const { return _width_; }
    size_t height(void) const { return _height_; }
    size_t depth(void) const { return _depth_; }
    size_t total_cells(void) const { return _total_cells_; }
    double voxel_size(void) const { return _voxel_size_; }
    double particle_radius(void) const { return _particle_radius_; }
    double sdf(const size_t index) const { return _sdf_[index]; }
    double sdf(const size_t i, const size_t j) const { return _sdf_[index_1D(i, j)]; }
    double u(const size_t dim, const size_t index) const { return _u_[dim][index]; }
    double u(const size_t dim, const size_t i, const size_t j) const
        { return _u_[dim][index_1D(i, j)]; }
    double old_u(const size_t dim, const size_t index) const { return _old_u_[dim][index]; }
    double old_u(const size_t dim, const size_t i, const size_t j) const
        { return _old_u_[dim][index_1D(i, j)]; }
    double temp_u(const size_t dim, const size_t index) const { return _temp_u_[dim][index]; }
    double temp_u(const size_t dim, const size_t i, const size_t j) const
        { return _temp_u_[dim][index_1D(i, j)]; }
    double w(const size_t dim, const size_t index) const { return _w_[dim][index]; }
    double w(const size_t dim, const size_t i, const size_t j) const
        { return _w_[dim][index_1D(i, j)]; }
    double pressure(const size_t index) const { return _pressure_[index]; }
    double pressure(const size_t i, const size_t j) const { return _pressure_[index_1D(i, j)]; }
    double ux(const size_t index) const { return _u_[0][index]; }
    double ux(const size_t i, const size_t j) const { return _u_[0][index_1D(i, j)]; }
    double uy(const size_t index) const { return _u_[1][index]; }
    double uy(const size_t i, const size_t j) const { return _u_[1][index_1D(i, j)]; }
    double uz(const size_t index) const { return _u_[2][index]; }
    double uz(const size_t i, const size_t j) const { return _u_[2][index_1D(i, j)]; }
    int layer(const size_t index) const { return _layer_[index]; }
    int layer(const size_t i, const size_t j) const { return _layer_[index_1D(i, j)]; }
    CellType type(const size_t index) const { return _type_[index]; }
    CellType type(const size_t i, const size_t j) const { return _type_[index_1D(i, j)]; }

    /**
    * Setter methods
    */
    void set_width(const size_t width);
    void set_height(const size_t height);
    void set_depth(const size_t depth);
    void set_u(const size_t dim, const size_t index, const double val) { _u_[dim][index] = val; }
    void set_u(const size_t dim, const size_t i, const size_t j, const double val)
        { _u_[dim][index_1D(i, j)] = val; }
    void set_old_u(const size_t dim, const size_t index, const double val)
        { _old_u_[dim][index] = val; } 
    void set_old_u(const size_t dim, const size_t i, const size_t j, const double val)
        { _old_u_[dim][index_1D(i, j)] = val; }
    void set_temp_u(const size_t dim, const size_t index, const double val)
        { _temp_u_[dim][index] = val; } 
    void set_temp_u(const size_t dim, const size_t i, const size_t j,const double val)
        { _temp_u_[dim][index_1D(i, j)] = val; }
    void set_w(const size_t dim, const size_t index, const double val) { _w_[dim][index] = val; }
    void set_w(const size_t dim, const size_t i, const size_t j, const double val)
        { _w_[dim][index_1D(i, j)] = val; }
    void set_pressure(const size_t index, const double val) { _pressure_[index] = val; }
    void set_pressure(const size_t i, const size_t j, const double val)
        { _pressure_[index_1D(i, j)] = val; }

    /**
    * Converts 2D indices into a 1D index. A MacGrid stores its data as 1D arrays, so a conversion
    * is often convenient.
    * 
    * @param i The index into the ith voxel in the x dimension.
    * @param j The index into the jth voxel in the y dimension.
    * @return The 1D index representing i, j.
    */
    int index_1D(const size_t i, const size_t j) const
    {
        if(i >= _width_ || j >= _height_)
        {
            return -1;
        }

        return j * _width_ + i;
    }

    /**
    * Linearly interpolates a velocity vector at the specified position on the MacGrid.
    *
    * @param x The x location of the desired velocity in world units.
    * @param y The y location of the desired velocity in world units. 
    * @param z The z location of the desired velocity in world units. NOTE--this parameter is
    *   not used in this 2D implementation.
    * @param u_type The type of velocity to interpolate,
    * @param[out] result The resulting interpolated velocity will be stored here.
    */
    void get_velocity(const double x, const double y, const double z,
                      const VelocityType u_type, UT_Vector3& result) const;

    /**
    * Trace a particle at a specified position using the MacGrid's current velocity field to the
    * resulting position.
    *
    * @param pos The start position of the particle.
    * @param t The timestep.
    * @param[out] result The resulting traced end position of the particle will be stored here.
    */
    void trace_particle(const UT_Vector3& pos, const double t, UT_Vector3& result) const;

    /**
    * Advect particles using the MacGrid's current velocity field to new positions.
    *
    * @param[in,out] particles The vector of particles to advect. These particles will be modified
    *   with the resulting positions.
    * @param t The timestep.
    */
    void advect_particles(std::vector<UT_Vector3>& particles, const double t) const;

    /**
    * Applies the specified gravity force uniformly to the fluid cells in the MacGrid.
    *
    * @param gravity The gravity force to apply in houdini units
    * @param t The timestep
    */
    void apply_gravity_forces(const double gravity, const double t);

    /**
    * Computes the pressure field and subtracts it from the velocity field to ensure that the
    * velocity field is divergence free.
    *
    * @param fluid_density The density of the fluid.
    * @param st_const The surface tension constant.
    * @param t The timestep
    */
    void pressure_projection(const double fluid_density, const double st_const, const double t);

    /**
    * Checks that the provided particles are within the bounds of this MacGrid and that the
    * particle velocities are not pointing into the boundaries
    *
    * @param[in,out] p_positions The positions of the particles to check. These positions will be
    *   modified if they are out of bounds.
    * @param[in,out] p_velocities The velocities of the particles to check. These velocities will
    *   modified if they point into the boundaries.
    */
    void enforce_particle_bounds(std::vector<UT_Vector3>& p_positions,
                                 std::vector<UT_Vector3>& p_velocities) const;

    /**
    * Sets the velocities at the MacGrid boundaries to zero so that fluid will not flow into the
    * boundaries.
    */
    void set_boundary_velocities(void);

    /**
    * Extrapolates the velocity in fluid cells into the surrounding cells.
    *
    * @param kcfl The size in number of voxels to build out the extrapolated velocity buffer.
    */
    void extrapolate_velocity(const size_t kcfl);

    /**
    * Updates the MacGrid velocity field with the provided particle velocities as part of the FLIP
    * algorithm.
    *
    * @param p_positions The positions of the particles.
    * @param p_velocities The velocities of the particles to tranfer to the grid.
    */
    void pvel_to_grid(const std::vector<UT_Vector3>& p_positions,
                      const std::vector<UT_Vector3>& p_velocities);

    /**
    * Updates the provided particle velocities with the velocity field stored in this MacGrid as
    * as part of FLIP algorithm.
    *
    * @param p_positions The positions of the particles to update.
    * @param[out] p_velocities The particle velocities to be updated.
    * @param flip_ratio The ratio of flip to pic. This affects how inviscid the resulting fluid is.
    */
    void grid_to_pvel(std::vector<UT_Vector3>& p_positions,
                      std::vector<UT_Vector3>& p_velocities,
                      double flip_ratio) const;

    /**
    * Updates the cell types in this MacGrid to reflect the provided particles including building a
    * buffer of air cells around the fluid cells.
    *
    * @param kcfl The size in number of voxels to build out the buffer of air cells around the
    *   fluid cells.
    * @param particles The fluid particles to update the grid with.
    */
    void update_buffer(const int kcfl, const std::vector<UT_Vector3>& particles);

    /**
    * Computes the divergence of the velocity at the given index.
    *
    * @param i The index into the ith voxel in the x dimension.
    * @param j The index into the jth voxel in the y dimension.
    * @return The divergence of the velocity field at index i, j.
    */
    double compute_divergence(const size_t i, const size_t j) const;

    /**
    * Computes the mean curvature of the signed distance field at the given index.
    *
    * @param i The index into the ith voxel in the x dimension.
    * @param j The index into the jth voxel in the y dimension.
    * @return The mean curvature of the signed distance field at index i, j
    */
    double compute_curvature(const double x, const double y) const;

    /**
    * Interpolates the signed distance to the fluid surface at the provided position.
    * 
    * @param x The x location of the desired velocity in world units.
    * @param y The y location of the desired velocity in world units.
    * @return The interpolated value of signed distance value at position x, y
    */
    double get_sd(double x, double y) const;

    /**
    * Interpolates the pressure at the provided position.
    *
    * @param x The x location of the desired velocity in world units.
    * @param y The y location of the desired velocity in world units.
    * @return The interpolated value of pressure at position x, y
    */
    double get_pressure(double x, double y) const;

    /**
    * Computes the signed distance field (sdf) for this MacGrid with the provided particles. The
    * sdf represents the distance to the fluid surface. A positive value indicates a distance
    * outside of the fluid, and a negative value is inside the fluid surface.
    *
    * @param The FLIP particle locations that represent the fluid.
    */
    void compute_sdf(const std::vector<UT_Vector3>& particles);

    /**
    * Returns the maximum value of velocity in the MacGrid.
    *
    * @return The maximum value of velocity in the grid.
    */
    double max_u(void) const;

    /**
    * Returns the minimum value of signed distance in the MacGrid.
    *
    * @return The minimum value of signed distance in the grid.
    */
    double min_sdf(void) const;

private:

    //voxel count variables
    size_t _width_;
    size_t _height_;
    size_t _depth_;
    size_t _total_cells_;

    //variables related to the voxel size
    double _voxel_size_;
    double _half_voxel_size_;
    double _particle_radius_;
    double _w_to_v_;

    //the number of times to sweep the signed distance field
    int _sdf_sweep_count_;

    //fields NOTE--the 3rd dimension in the vectors is unused in this 2D implementation
    double* _sdf_;
    double* _pressure_;
    double* _u_[3];
    double* _old_u_[3];
    double* _temp_u_[3];
    double* _w_[3];
    int* _layer_;
    CellType* _type_;

    //pressure solve variables
    Eigen::SparseMatrix<double>* _A_;
    Eigen::VectorXd* _p_;
    Eigen::VectorXd* _b_;

    /**
    * Using the current configuration of the grid, the pressure field is computed and stored in
    * the grid.
    *
    * @param d The (fluid density * voxel size) / timestep.
    * @param st_const The surface tension constant.
    * @param t The timestep
    */
    void _compute_pressures(const double d, const double st_const);

    /**
    * Stores the pressures in this MacGrid that have been previously computed in the _p_ vector.
    */
    void _store_pressures(void);

    /**
    * Subtracts the gradient of the pressure field from the velocity field.
    *
    * @param d The timestep / (fluid_density * _voxel_size_)
    * @param st_const The surface tension constant.
    */
    void _apply_pressure_gradient(const double d, const double st_const);

    /**
    * Gives each fluid cell a label cooresponding to its row in the pressure solve matrix
    *
    * @return The id of the labeled fluid cell.
    */
    size_t _label_fluid_cells(void);

    /**
    * Solves the previously initialized pressure solve matrix with the resulting pressures in _p_.
    */
    void _solve_pressure_matrix(void);

    /**
    * Builds the pressure solve matrix _A_ with the current configuration of the MacGrid.
    *
    * @param d (fluid_density * _voxel_size_) / timestep
    * @param st_const The surface tension constant
    */
    void _build_pressure_matrix(const double d, const double st_const);

    /**
    * Swaps the temporary velocity field. The values of _u_ will be set to the values currently in
    * _temp_u_.
    */
    void _swap_temp_u(void);

    /**
    * Divides the values in _u_ by the weights in _w_. This is used in the FLIP pvel_to_grid
    * method.
    */
    void _divide_vel_by_weights(void);

    /**
    * Gets the voxel indices of the 4 neighbors to the provided index.
    *
    * @param index The index to get the neighbors for.
    * @param[out] indices The array where the neighbor indices will be stored.
    */
    void _get_neighbor_indices(const size_t index, int* indices) const;

    /**
    * Gets the voxel indices of the 4 neighbors to the provided location.
    *
    * @param i The index into the ith voxel in the x dimension.
    * @param j The index into the jth voxel in the y dimension.
    * @param[out] indices The array where the neighbor voxel indices will be stored.
    */
    void _get_neighbor_indices(const size_t i, const size_t j, int* indices) const;

    /**
    * Gets the index of the neighbor given a dimension and direction (negative or positive).
    *
    * @param i The index into the ith voxel in the x dimension.
    * @param j The index into the jth voxel in the y dimension.
    * @param dim The dimension to get the neighbor in (0=x, 1=y)
    * @param neg Flag indicating whether to get the neighber in the negative (true) direction or
    *   positive (false).
    * @return The neighbor index.
    */
    int _get_neighbor_index(const size_t i, const size_t j, const size_t dim,
                            const bool neg) const;

    /**
    * Set the layer field for all cells in the grid.
    *
    * @param layer The value to set for the layer field for all cells.
    */
    void _set_layer_all(const int layer);

    /**
    * Interpolates a scalar value of a specified type at the desired position.
    *
    * @param x The x location to interpolate in world units.
    * @param y The y location to interpolate in world units. 
    * @param z The z location to interpolate in world units. NOTE--this parameter is
    *   not used in this 2D implementation.
    * @param type The type to interpolate.
    * @return The interpolated value.
    */
    double _get_interp_val(const double x, const double y, const double z,
                           const ScalarType type) const;

    /**
    * Gets the scalars of a specified type on the grid at the specified voxel indices.
    *
    * @param i The index into the ith voxel in the x dimension.
    * @param j The index into the jth voxel in the y dimension.
    * @param k The index into the jth voxel in the z dimension. NOTE--this parameter is not used
    *   in this 2D implementation.
    * @param type The type of scalar to get.
    * @param[out] The array to store the resulting scalars in.
    */
    void _get_interp_scalars(const size_t i, const size_t j, const size_t k,
                             const ScalarType type, double* result) const;

    /**
    * Get the interpolation weights given a position in world units.
    *
    * @param x The x location to interpolate in world units.
    * @param y The y location to interpolate in world units. 
    * @param z The z location to interpolate in world units. NOTE--this parameter is
    *   not used in this 2D implementation.
    * @param[out] result The array to store the resulting weights in.
    */
    void _get_interp_weights(const double x, const double y, const double z, double* result) const;

    /**
    * Sets all the grid field values to default values.
    */
    void _reset_grid(void);

    /**
    * Sets the velocity and weights fields to 0.0.
    */
    void _reset_grid_velocities(void);

    /**
    * Propagates the signed distance field given an initialized signed distance field. This is
    * part of the signed distance field computation.
    *
    * @param iterations The number of voxels to propagate the sdf
    * @param neg Flag indicating whether to propagate the sdf inside the surface (true) or 
    *   outside (false)
    */
    void _sweep_sdf(const size_t iterations, const bool neg);

    /**
    * Checks if the specified voxel has a neighbor in the desired layer.
    *
    * @param i The index into the ith voxel in the x dimension.
    * @param j The index into the jth voxel in the y dimension.
    * @param mlayer The layer to check for.
    * @return Whether voxel i, j has a neighbor in layer mylayer (true) or not (false)
    */
    bool _has_neighbor_in_layer(const size_t i, const size_t j, const int mlayer) const;

    /**
    * Returns the minimum signed distance value in a neighbor to voxel i, j in layer mlayer.
    *
    * @param i The index into the ith voxel in the x dimension.
    * @param j The index into the jth voxel in the y dimension.
    * @param mlayer The layer to check for.
    * @return The minimum signed distance value in a neighbor to voxel i, j in layer mlayer,
    */
    double _min_neighbor_sdf_in_layer(const size_t i, const size_t j, const int mlayer) const;

    /**
    * Returns the maximum signed distance value in a neighbor to voxel i, j in layer mlayer.
    *
    * @param i The index into the ith voxel in the x dimension.
    * @param j The index into the jth voxel in the y dimension.
    * @param mlayer The layer to check for.
    * @return The maximum signed distance value in a neighbor to voxel i, j in layer mlayer,
    */
    double _max_neighbor_sdf_in_layer(const size_t i, const size_t j, const int mlayer) const;
};


//Generates a random double between fMin and fMax
double fRand(const double fMin, const double fMax);

#endif /* MACGRID_H_ */
