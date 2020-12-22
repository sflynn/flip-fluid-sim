/*
 * macgrid.C
 *
 *  Created on: September 5, 2018
 *      Author: Sean Flynn
 */

#include "macgrid.h"

#include <GU/GU_Detail.h>
#include <GEO/GEO_PrimPart.h>

using namespace std;

//constants
const int LARGE_NUMBER = 1000000000;
const int LARGE_NEG_NUMBER = LARGE_NUMBER * -1;
const double EPSILON = .00001;
const double THETA_MIN = .03;

//Generates a random double between fMin and fMax
double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

auto array_of_vecs = [] (size_t count) {
        return std::array<std::vector<double>, 3>{
                std::vector<double>(count), std::vector<double>(count),
                std::vector<double>(count) };
    };


//Constructs a MacGrid object.
MacGrid::MacGrid(size_t width, size_t height, size_t depth,
                 double voxel_size, double particle_radius) :
        _width_(width), _height_(height), _depth_(depth), _total_cells_(width * height * depth),
        _voxel_size_(voxel_size), _half_voxel_size_(voxel_size / 2.0),
        _particle_radius_(particle_radius), _w_to_v_(1.0 / voxel_size),
        _sdf_sweep_count_((int)ceil(max(width, height) * .7)),
        _sdf_(_total_cells_), _pressure_(_total_cells_),
        _u_(array_of_vecs(_total_cells_)),
        _old_u_(array_of_vecs(_total_cells_)),
        _temp_u_(array_of_vecs(_total_cells_)),
        _w_(array_of_vecs(_total_cells_)),
        _layer_(_total_cells_), _type_(_total_cells_)
{
    _reset_grid();
}

void MacGrid::set_width(size_t width)
{
    _width_ = width;
    _total_cells_ = _width_ * _height_ * _depth_;
}

void MacGrid::set_height(size_t height)
{
    _height_ = height;
    _total_cells_ = _width_ * _height_ * _depth_;
}

void MacGrid::set_depth(size_t depth)
{
    _depth_ = depth;
    _total_cells_ = _width_ * _height_ * _depth_;
}

//Linearly interpolates a velocity vector at the specified position on the MacGrid.
void MacGrid::get_velocity(double x, double y, double z,
                           VelocityType u_type, UT_Vector3& result) const
{
    double hdx = max(x - _half_voxel_size_, 0.0);
    double hdy = max(y - _half_voxel_size_, 0.0);
    double hdz = max(z - _half_voxel_size_, 0.0);

    if(x <= 0 || x >= _width_)
        result[0] = 0.0;
    else
    {
        if(u_type == CUR)
            result[0] = _get_interp_val(x, hdy, hdz, UX);
        else
            result[0] = _get_interp_val(x, hdy, hdz, UX_TEMP);
    }

    if(y <= 0 || y >= _height_)
        result[1] = 0.0;
    else
    {
        if(u_type == CUR)
            result[1] = _get_interp_val(hdx, y, hdz, UY);
        else
            result[1] = _get_interp_val(hdx, y, hdz, UY_TEMP);
    }

    result[2] = 0.0;
}

//Trace a particle at a specified position using the MacGrid's current velocity field to the
//resulting position.
void MacGrid::trace_particle(const UT_Vector3& pos, double t, UT_Vector3& result) const
{
    //RK2
    UT_Vector3 vel;

    get_velocity(pos[0], pos[1], pos[2], CUR, vel);

    double ht = .5 * t;
    get_velocity(pos[0] + ht * vel[0], pos[1] + ht * vel[1], pos[2] + ht * vel[2], CUR, vel);

    result = pos + vel * t;
}

//Advect particles using the MacGrid's current velocity field to new positions.
void MacGrid::advect_particles(vector<UT_Vector3>& particles, double t) const
{
    //TODO: make this parallel
    for(size_t i = 0; i < particles.size(); i++)
    {
        UT_Vector3 new_pos;
        trace_particle(particles[i], t, new_pos);
        particles[i][0] = new_pos[0];
        particles[i][1] = new_pos[1];
        particles[i][2] = new_pos[2];
    }
}

//Applies the specified gravity force uniformly to the fluid cells in the MacGrid.
void MacGrid::apply_gravity_forces(double gravity, double t)
{
    int index;
    double vg = gravity * t;
    for(size_t j = 0; j < _height_; j++)
    {
        for(size_t i = 0; i < _width_; i++)
        {
            index = index_1D(i, j);
            if(type(index) == FLUID || type(index) == AIR)
                set_u(1, index, uy(index) + vg);
        }
    }
}

//Gives each fluid cell a label cooresponding to its row in the pressure solve matrix
size_t MacGrid::_label_fluid_cells(void)
{
    size_t cell_id = 0;
    for(size_t i = 0; i < _total_cells_; i++)
    {
        if(_type_[i] == FLUID)
            _layer_[i] = cell_id++;
        else
            _layer_[i] = -1;
    }

    return cell_id;
}

//Solves the previously initialized pressure solve matrix with the resulting pressures in _p_.
void MacGrid::_solve_pressure_matrix(void)
{
    //cout << "Pressure solve matrix-------------------------" << endl;

    //cout << "left hand side: " << endl;
    //cout << *(this->_A_) << endl;

    //cout << "right hand side: " << endl;
    //cout << *(this->_b_) << endl;

    //cout << "\t\t\tMatrix solve..." << endl;

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > cg;
    cg.compute(*this->_A_);
    *this->_p_ = cg.solve(*this->_b_);

    //cout << "result: \n" << endl;
    //cout << *(this->_p_) << endl;

    //cout << "----------------------------------------------" << endl;

    //cout << "\t\t\t\t# iterations:    " << cg.iterations() << endl;
    //cout << "\t\t\t\testimated error: " << cg.error() << endl;
}

//Builds the pressure solve matrix _A_ with the current configuration of the MacGrid.
void MacGrid::_build_pressure_matrix(double d, double st_const)
{
    double ghost, g_p, rhs, theta, k, divergence, x, y;
    size_t non_solid_neighbors;
    int index;
    int n_indices[4];

    for(size_t j = 0; j < _height_; j++)
    {
        for(size_t i = 0; i < _width_; i++)
        {
            index = index_1D(i, j);
            if(_type_[index] != FLUID)
                continue;

            non_solid_neighbors = 0;
            ghost = 0.0;
            rhs = 0.0;

            _get_neighbor_indices(index, n_indices);
            for(size_t n = 0; n < 4; n++)
            {
                if(n_indices[n] == -1 || _type_[n_indices[n]] == SOLID)
                    continue;

                ++non_solid_neighbors;
                if(_type_[n_indices[n]] == FLUID)
                {
                    _A_->coeffRef(_layer_[index], _layer_[n_indices[n]]) = 1;
                }
                else
                {
                    theta = max(THETA_MIN, (_sdf_[index] /
                            (_sdf_[index] - _sdf_[n_indices[n]]))) * _voxel_size_;

                    if(n == 0)//left
                    {
                        x = i + _half_voxel_size_ - theta;
                        y = j + _half_voxel_size_;
                    }
                    else if(n == 1)//down
                    {
                        x = i + _half_voxel_size_;
                        y = j + _half_voxel_size_ - theta;
                    }
                    else if(n == 2)//right
                    {
                        x = i + _half_voxel_size_ + theta;
                        y = j + _half_voxel_size_;
                    }
                    else//up
                    {
                        x = i + _half_voxel_size_;
                        y = j + _half_voxel_size_ + theta;
                    }

                    k = compute_curvature(x, y);
                    ghost += _sdf_[n_indices[n]];
                    g_p = (st_const * k) / theta;
                    rhs += g_p;
                    _pressure_[n_indices[n]] = g_p;
                }
            }

            double val = (non_solid_neighbors * -1.0) + (ghost / _sdf_[index]);
            _A_->coeffRef(_layer_[index], _layer_[index]) = val;
            divergence = compute_divergence(i, j);
            _b_->operator()(_layer_[index]) = (d * divergence) - rhs;
        }
    }
}

//Computes the pressure field and subtracts it from the velocity field to ensure that the
//velocity field is divergence free.
void MacGrid::pressure_projection(double fluid_density, double st_const, double t)
{
    double d = (fluid_density * _voxel_size_) / t;
    _compute_pressures(d, st_const);
    d = t / (fluid_density * _voxel_size_);
    _apply_pressure_gradient(d, st_const);
}

//Using the current configuration of the grid, the pressure field is computed and stored in the
//grid.
void MacGrid::_compute_pressures(double d, double st_const)
{
    size_t fluid_cell_count = _label_fluid_cells();

    //allocate memory for matrix solve
    _A_ = new Eigen::SparseMatrix<double>(fluid_cell_count, fluid_cell_count);
    _p_ = new Eigen::VectorXd(fluid_cell_count);
    _b_ = new Eigen::VectorXd(fluid_cell_count);
    _A_->reserve(Eigen::VectorXi::Constant(fluid_cell_count, 20));

    _build_pressure_matrix(d, st_const);
    _solve_pressure_matrix();
    _store_pressures();

    delete _A_;
    delete _p_;
    delete _b_;
}

//Subtracts the gradient of the pressure field from the velocity field.
void MacGrid::_store_pressures(void)
{
    for(size_t index = 0; index < _total_cells_; index++)
        if(_layer_[index] != -1)
            _pressure_[index] = (double)(*_p_)[_layer_[index]];
}

//Subtracts the gradient of the pressure field from the velocity field.
void MacGrid::_apply_pressure_gradient(double d, double st_const)
{
    int index, n_index;
    double p_for, p_back, theta, k, grad;
    for(size_t dim = 0; dim < 2; dim++)
    {
        for(size_t j = 0; j < _height_; j++)
        {
            for(size_t i = 0; i < _width_; i++)
            {
                index = index_1D(i, j);
                n_index = _get_neighbor_index(i, j, dim, true);
                if(n_index != -1)
                {
                    if(_type_[index] != FLUID && _type_[n_index] != FLUID)
                        continue;

                    if(_type_[index] == FLUID)
                    {
                        p_for = _pressure_[index];
                        if(_type_[n_index] == FLUID)
                        {
                            p_back = _pressure_[n_index];
                        }
                        else if(_type_[n_index] == AIR)
                        {
                            theta = max(THETA_MIN, (_sdf_[index] /
                                    (_sdf_[index] - _sdf_[n_index]))) * _voxel_size_;
                            if(dim == 0)
                            {
                                k = compute_curvature(i + _half_voxel_size_ - theta,
                                                      j + _half_voxel_size_);
                            }
                            else
                            {
                                k = compute_curvature(i + _half_voxel_size_,
                                                      j + _half_voxel_size_ - theta);
                            }
                            
                            p_back = ((st_const * k) / theta) +
                                     ((_sdf_[n_index] / _sdf_[index]) * p_for);
                        }

                        grad = d * (p_for - p_back);
                        _u_[dim][index] -= grad;
                    }
                    else if(_type_[index] == AIR)
                    {
                        if(_type_[n_index] == FLUID)
                        {
                            p_back = _pressure_[n_index];
                            theta = max(THETA_MIN, (_sdf_[n_index] /
                                    (_sdf_[n_index] - _sdf_[index]))) * _voxel_size_;
                            if(dim == 0)
                            {
                                k = compute_curvature(i - _half_voxel_size_ + theta,
                                                      j + _half_voxel_size_);
                            }
                            else
                            {
                                k = compute_curvature(i + _half_voxel_size_,
                                                      j - _half_voxel_size_ + theta);
                            }

                            p_for = ((st_const * k) / theta) +
                                    ((_sdf_[index] / _sdf_[n_index]) * p_back);
                        }

                        grad = d * (p_for - p_back);
                        _u_[dim][index] -= d * (p_for - p_back);
                    }
                }
            }
        }
    }
}

//Checks that the provided particles are within the bounds of this MacGrid and that the
//particle velocities are not pointing into the boundaries
void MacGrid::enforce_particle_bounds(vector<UT_Vector3>& p_positions,
                                      vector<UT_Vector3>& p_velocities) const
{
    double mult = -.01;
    size_t max_dims[3] = {_width_, _height_, _depth_};

    for(size_t i = 0; i < p_positions.size(); i++)
    {
        for(size_t dim = 0; dim < 2; dim++)
        {
            if(p_positions[i][dim] < EPSILON)
            {
                p_positions[i][dim] = EPSILON;
                if(p_velocities[i][dim] < 0.0)
                    p_velocities[i][dim] *= mult;
            }

            if(p_positions[i][dim] >= max_dims[dim])
            {
                p_positions[i][dim] = max_dims[dim] - EPSILON;
                if(p_velocities[i][dim] > 0.0)
                    p_velocities[i][dim] *= mult;
            }
        }
    }
}

//Sets the velocities at the MacGrid boundaries to zero so that fluid will not flow into the
//boundaries.
void MacGrid::set_boundary_velocities(void)
{
    for(size_t j = 0; j < _height_; j++)
        set_u(0, 0, j, 0.0);

    for(size_t i = 0; i < _width_; i++)
        set_u(1, i, 0, 0.0);
}

//Extrapolates the velocity in fluid cells into the surrounding cells.
void MacGrid::extrapolate_velocity(size_t kcfl)
{
    double avg_u[3];
    size_t u_counts[3];
    int n_indices[4];

    for(size_t counter = 0; counter < kcfl; ++counter)
    {
        for(size_t dim = 0; dim < 2; dim++)
        {
            for(size_t index = 0; index < _total_cells_; index++)
            {
                _get_neighbor_indices(index, n_indices);

                _temp_u_[dim][index] = _u_[dim][index];
                avg_u[dim] = 0.0;
                u_counts[dim] = 0;
                for(size_t n = 0; n < 4; n++)
                {
                    if(n_indices[n] == -1)
                        continue;

                    if(_u_[dim][n_indices[n]] != 0.0)
                    {
                        avg_u[dim] += _u_[dim][n_indices[n]];
                        u_counts[dim]++;
                    }
                }

                if(_temp_u_[dim][index] == 0.0 && u_counts[dim] > 0)
                {
                    _temp_u_[dim][index] = avg_u[dim] / (double)u_counts[dim];
                }
            }
        }
        _swap_temp_u();
    }
}

//Updates the MacGrid velocity field with the provided particle velocities as part of the FLIP
//algorithm.
void MacGrid::pvel_to_grid(const vector<UT_Vector3>& p_positions,
                           const vector<UT_Vector3>& p_velocities)
{
    _reset_grid_velocities();

    UT_Vector3 px, py;
    int pi[2], pj[2];
    size_t x_index, y_index;
    double wx[4], wy[4];

    int indices[4][2] =
    {
            {0, 0},
            {1, 0},
            {0, 1},
            {1, 1}
    };

    for(size_t index = 0; index < p_positions.size(); index++)
    {
        px[0] = p_positions[index][0];
        px[1] = p_positions[index][1] - .5;
        pi[0] = floor(px[0]);
        pi[1] = floor(px[1]);
        _get_interp_weights(px[0], px[1], 0.0, wx);

        py[0] = p_positions[index][0] - .5;
        py[1] = p_positions[index][1];
        pj[0] = floor(py[0]);
        pj[1] = floor(py[1]);
        _get_interp_weights(py[0], py[1], 0.0, wy);

        for(size_t i = 0; i < 4; i++)
        {
            x_index = index_1D(pi[0] + 1 * indices[i][0], pi[1] + 1 * indices[i][1]);
            y_index = index_1D(pj[0] + 1 * indices[i][0], pj[1] + 1 * indices[i][1]);

            if(x_index != -1 && _type_[x_index] != SOLID)
            {
                set_u(0, x_index, u(0, x_index) + p_velocities[index][0] * wx[i]);
                set_w(0, x_index, w(0, x_index) + wx[i]);
            }
            if(y_index != -1 && _type_[y_index] != SOLID)
            {
                set_u(1, y_index, u(1, y_index) + p_velocities[index][1] * wy[i]);
                set_w(1, y_index, w(1, y_index) + wy[i]);
            }   
        }
    }
    _divide_vel_by_weights();
}

//Updates the provided particle velocities with the velocity field stored in this MacGrid as
//as part of FLIP algorithm.
void MacGrid::grid_to_pvel(vector<UT_Vector3>& p_positions, vector<UT_Vector3>& p_velocities,
                           double flip_ratio)
{
    //compute the difference between u and old u and store it in temp u
    for(size_t dim = 0; dim < 2; dim++)
        for(size_t i = 0; i < _total_cells_; i++)
            _temp_u_[dim][i] = _u_[dim][i] - _old_u_[dim][i];

    double pic_ratio = 1.0 - flip_ratio;

    UT_Vector3 u_cur, u_dif;
    for(size_t index = 0; index < p_positions.size(); index++)
    {
        get_velocity(p_positions[index][0],
                     p_positions[index][1],
                     p_positions[index][2],
                     CUR, u_cur);
        get_velocity(p_positions[index][0],
                     p_positions[index][1],
                     p_positions[index][2],
                     TEMP, u_dif);

        p_velocities[index][0] += u_dif[0];
        p_velocities[index][1] += u_dif[1];

        p_velocities[index][0] *= flip_ratio;
        p_velocities[index][1] *= flip_ratio;

        p_velocities[index][0] += (u_cur[0] * pic_ratio);
        p_velocities[index][1] += (u_cur[1] * pic_ratio);
    }
}

//Divides the values in _u_ by the weights in _w_. This is used in the FLIP pvel_to_grid method.
void MacGrid::_divide_vel_by_weights(void)
{
    for(size_t dim = 0; dim < 2; dim++)
    {
        for(size_t i = 0; i < _total_cells_; i++)
        {
            if(_u_[dim][i] != 0.0 && _w_[dim][i] != 0.0)
            {
                set_u(dim, i, u(dim, i) / w(dim, i));
            }
            _old_u_[dim][i] = _u_[dim][i];
        }
    }
}

//Swaps the temporary velocity field. The values of _u_ will be set to the values currently in
//_temp_u_.
void MacGrid::_swap_temp_u(void)
{
    for(size_t dim = 0; dim < 2; dim++)
        for(size_t i = 0; i < _total_cells_; i++)
            _u_[dim][i] = _temp_u_[dim][i];
}

//Gets the voxel indices of the 4 neighbors to the provided index.
void MacGrid::_get_neighbor_indices(size_t index, int* indices) const
{
    size_t j = floor((double)index / (double)_width_);
    size_t i = index - j * _width_;
    _get_neighbor_indices(i, j, indices);
}

//Gets the voxel indices of the 4 neighbors to the provided location.
void MacGrid::_get_neighbor_indices(size_t i, size_t j, int* indices) const
{
    if(i > 0)
        indices[0] = index_1D(i - 1, j);
    else
        indices[0] = -1;

    if(j > 0)
        indices[1] = index_1D(i, j - 1);
    else
        indices[1] = -1;

    indices[2] = index_1D(i + 1, j);
    indices[3] = index_1D(i, j + 1);
}

//Gets the index of the neighbor given a dimension and direction (negative or positive).
int MacGrid::_get_neighbor_index(size_t i, size_t j,
                                 size_t dim, bool neg) const
{
    if(dim == 0)
    {
        if(neg == true)
            return index_1D(i - 1, j);
        else
            return index_1D(i + 1, j);
    }
    else if(dim == 1)
    {
        if(neg == true)
            return index_1D(i, j - 1);
        else
            return index_1D(i, j + 1);
    }
    return -1;
}

//Updates the cell types in this MacGrid to reflect the provided particles including building a
//buffer of air cells around the fluid cells.
void MacGrid::update_buffer(int kcfl, const vector<UT_Vector3>& particles)
{
    _set_layer_all(-1);

    for(size_t i = 0; i < _total_cells_; i++)
    {
        if(_sdf_[i] < 0)
        {
            _type_[i] = FLUID;
            _layer_[i] = 0;
        }
    }

    int max_layer = max(2, kcfl);
    int n_indices[4];
    for(size_t index = 1; index < max_layer; index++)
    {
        for(size_t i = 0; i < _total_cells_; i++)
        {
            if(_layer_[i] == index - 1 && (_type_[i] == FLUID || _type_[i] == AIR))
            {
                _get_neighbor_indices(i, n_indices);
                
                for(size_t n = 0; n < 4; n++)
                {
                    if(n_indices[n] != -1 && _layer_[n_indices[n]] == -1 && 
                       _type_[n_indices[n]] != SOLID)
                    {
                        _type_[n_indices[n]] = AIR;
                        _layer_[n_indices[n]] = index;
                    }
                }
            }
        }
    }

    for(size_t i = 0; i < _total_cells_; i++)
    {
        if(_type_[i] != SOLID && _layer_[i] == -1)
        {
            _type_[i] = UNUSED;
            for(size_t dim = 0; dim < 2; dim++)
                set_u(dim, i, 0.0);
        }
    }
}

//Set the layer field for all cells in the grid.
void MacGrid::_set_layer_all(int layer)
{
    for(size_t i = 0; i < _total_cells_; i++)
        _layer_[i] = layer;
}

//Computes the divergence of the velocity at the given index.
double MacGrid::compute_divergence(size_t i, size_t j) const
{
    int index = index_1D(i, j);
    int n;
    double u[2], ux, uy;
    u[0] = 0.0;
    u[1] = 0.0;
    ux = 0.0;
    uy = 0.0;

    if(index != -1 && _type_[index] != SOLID)
    {
        n = index_1D(i - 1, j);
        if(n != -1 && _type_[n] != SOLID)
            u[0] = _u_[0][index];

        n = index_1D(i, j - 1);
        if(n != -1 && _type_[n] != SOLID)
            u[1] = _u_[1][index];
    }

    n = index_1D(i + 1, j);
    if(n != -1 && _type_[n] != SOLID)
        ux = _u_[0][n];

    n = index_1D(i, j + 1);
    if(n != -1 && _type_[n] != SOLID)
        uy = _u_[1][n];

    return ux - u[0] + uy - u[1];
}

//Computes the mean curvature of the signed distance field at the given index.
double MacGrid::compute_curvature(double x, double y) const
{
    double sd, sdl, sdr, sdd, sdu;

    sd = get_sd(x, y);
    sdl = get_sd(x - _voxel_size_, y);
    sdd = get_sd(x, y - _voxel_size_);
    sdr = get_sd(x + _voxel_size_, y);
    sdu = get_sd(x, y + _voxel_size_);
    return -4.0 * sd + sdl + sdd + sdr + sdu;
}

//Interpolates the signed distance to the fluid surface at the provided position.
double MacGrid::get_sd(double x, double y) const
{
    x = max(x - _half_voxel_size_, 0.0);
    x = min(x, _width_ - EPSILON);
    y = max(y - _half_voxel_size_, 0.0);
    y = min(y, _height_ - EPSILON);

    return _get_interp_val(x, y, 0.0, SD);
}

//Interpolates the pressure at the provided position.
double MacGrid::get_pressure(double x, double y) const
{
    x = max(x - _half_voxel_size_, 0.0);
    x = min(x, _width_ - EPSILON);
    y = max(y - _half_voxel_size_, 0.0);
    y = min(y, _height_ - EPSILON);

    return _get_interp_val(x, y, 0.0, PRESSURE);
}

//Interpolates a scalar value of a specified type at the desired position.
double MacGrid::_get_interp_val(double x, double y, double z,
                                ScalarType type) const
{
    double weights[4];
    _get_interp_weights(x, y, z, weights);

    int i = floor(x);
    int j = floor(y);
    int k = floor(z);

    double vals[4];
    _get_interp_scalars(i, j, k, type, vals);

    return weights[0] * vals[0] +
           weights[1] * vals[1] +
           weights[2] * vals[2] +
           weights[3] * vals[3];
}

//Gets the scalars of a specified type on the grid at the specified voxel indices.
void MacGrid::_get_interp_scalars(size_t i, size_t j, size_t k,
                                  ScalarType type, double* result) const
{
    if(type == SD || type == PRESSURE)
    {
        if(type == SD)
            result[0] = sdf(i, j);
        else
            result[0] = pressure(i, j);
        
        if(i + 1 >= _width_)
            result[1]  = result[0];
        else
        {
            if(type == SD)
                result[1] = sdf(i + 1, j);
            else
                result[1] = pressure(i + 1, j);
        }

        if(j + 1 >= _height_)
            result[2] = result[0];
        else
        {
            if(type == SD)
                result[2] = sdf(i, j + 1);
            else
                result[2] = pressure(i, j + 1);
        }

        if(i + 1 >= _width_ || j + 1 >= _height_)
        {
            if(i + 1 >= _width_ && j + 1 >= _height_)
                result[3] = result[0];
            else if(i + 1 >= _width_)
                result[3] = result[2];
            else if(j + 1 >= _height_)
                result[3] = result[1];
        }
        else
        {
            if(type == SD)
                result[3] = sdf(i + 1, j + 1);
            else
                result[3] = pressure(i + 1, j + 1);
        }
    }
    else if(type == UX || type == UX_TEMP)
    {
        if(type == UX)
            result[0] = u(0, i, j);
        else
            result[0] = temp_u(0, i, j);

        if(i + 1 >= _width_)
            result[1] = 0.0;
        else
        {
            if(type == UX)
                result[1] = u(0, i + 1, j);
            else
                result[1] = temp_u(0, i + 1, j);
        }

        if(j + 1 >= _height_)
            result[2] = result[0];
        else
        {
            if(type == UX)
                result[2] = u(0, i, j + 1);
            else
                result[2] = temp_u(0, i, j + 1);
        }

        if(i + 1 >= _width_ || j + 1 >= _height_)
        {
            if(i + 1 >= _width_)
                result[3] = 0.0;
            else if(j + 1 >= _height_)
            {
                result[3] = result[1];
            }
        }
        else
        {
            if(type == UX)
                result[3] = u(0, i + 1, j + 1);
            else
                result[3] = temp_u(0, i + 1, j + 1);
        }
    }
    else if(type == UY || type == UY_TEMP)
    {
        if(type == UY)
            result[0] = u(1, i, j);
        else
            result[0] = temp_u(1, i, j);

        if(i + 1 >= _width_)
            result[1] = result[0];
        else
        {
            if(type == UY)
                result[1] = u(1, i + 1, j);
            else
                result[1] = temp_u(1, i + 1, j);
        }

        if(j + 1 >= _height_)
            result[2] = 0.0;
        else
        {
            if(type == UY)
                result[2] = u(1, i, j + 1);
            else
                result[2] = temp_u(1, i, j + 1);
        }

        if(i + 1 >= _width_ || j + 1 >= _height_)
        {
            if(j + 1 >= _height_)
                result[3] = 0.0;
            else if(i + 1 >= _width_)
            {
                result[3] = result[1];
            }
        }
        else
        {
            if(type == UY)
                result[3] = u(1, i + 1, j + 1);
            else
                result[3] = temp_u(1, i + 1, j + 1);
        }
    }
}

//Get the interpolation weights given a position in Houdini units.
void MacGrid::_get_interp_weights(double x, double y, double z,
                                  double* result) const
{
    int i = floor(x);
    int j = floor(y);

    result[0] = (i + 1 - x) * (j + 1 - y);
    result[1] = (x - i) * (j + 1 - y);
    result[2] = (i + 1 - x) * (y - j);
    result[3] = (x - i) * (y - j);
}

//Returns the maximum value of velocity in the MacGrid.
double MacGrid::max_u(void) const
{
    double max = LARGE_NEG_NUMBER;
    for(size_t dim = 0; dim < 2; dim++)
        for(size_t i = 0; i < _total_cells_; i++)
            if(abs(_u_[dim][i]) > max)
                max = abs(_u_[dim][i]);

    return max;
}

//Returns the minimum value of signed distance in the MacGrid.
double MacGrid::min_sdf(void) const
{
    double min = LARGE_NUMBER;
    for(size_t i = 0; i < _total_cells_; i++)
        if(_sdf_[i] < min)
            min = _sdf_[i];

    return min;
}

//Sets the velocity and weights fields to 0.0.
void MacGrid::_reset_grid_velocities(void)
{
    for(size_t dim = 0; dim < 2; dim++)
    {
        for(size_t i = 0; i < _total_cells_; i++)
        {
            _u_[dim][i] = 0.0;
            _temp_u_[dim][i] = 0.0;
            _old_u_[dim][i] = 0.0;
            _w_[dim][i] = 0.0;
        }
    }
}

//Sets all the grid field values to default values.
void MacGrid::_reset_grid(void)
{
    _reset_grid_velocities();
    for(size_t i = 0; i < _total_cells_; i++)
    {
        _sdf_[i] = LARGE_NUMBER;
        _pressure_[i] = 0.0;
        _layer_[i] = -1;
        _type_[i] = UNUSED;
    }
}

//Computes the signed distance field (sdf) for this MacGrid with the provided particles. The
//sdf represents the distance to the fluid surface. A positive value indicates a distance
//outside of the fluid, and a negative value is inside the fluid surface.
void MacGrid::compute_sdf(const vector<UT_Vector3>& particles)
{
    UT_Vector3 pos;

    double d, dr, du, dru;
    size_t index;

    _reset_grid();

    int i, j;
    for (size_t p = 0; p < particles.size(); p++)
    {
        pos = particles[p];

        i = max(floor(pos[0] - .5), 0.0);
        j = max(floor(pos[1] - .5), 0.0);

        //set the cell that contains this particle's distance
        d = sqrt(pow((i) + .5 - pos[0], 2) +
                 pow((j) + .5 - pos[1], 2)) - _particle_radius_;
        index = index_1D(i, j);
        if(d < _sdf_[index])
        {
            _sdf_[index] = d;
        }

        //set the right neighbor's distance
        if(i < _width_ - 1)
        {
            dr = sqrt(pow((i + 1.5) - pos[0], 2) +
                      pow((j) + .5 - pos[1], 2)) - _particle_radius_;
            index = index_1D(i + 1, j);
            if(dr < _sdf_[index])
                _sdf_[index] = dr;
        }

        //set the up neighbor's distance
        if(j < _height_ - 1)
        {
            du = sqrt(pow((i) + .5 - pos[0], 2) +
                      pow((j + 1.5) - pos[1], 2)) - _particle_radius_;
            index = index_1D(i, j + 1);
            if(du < _sdf_[index])
                _sdf_[index] = du;

            //set the up-right neighbor's distance
            if(i < _width_ - 1)
            {   
                dru = sqrt(pow((i + 1.5) - pos[0], 2) +
                           pow((j + 1.5) - pos[1], 2)) - _particle_radius_;
                index = index_1D(i + 1, j + 1);
                if(dru < _sdf_[index])
                    _sdf_[index] = dru;
            }
        }
    }

    _sweep_sdf(_sdf_sweep_count_, true);
    _sweep_sdf(_sdf_sweep_count_, false);
}

//Propagates the signed distance field given an initialized signed distance field. This is
//part of the signed distance field computation.
void MacGrid::_sweep_sdf(size_t iterations, bool neg)
{
    size_t index;

    for(size_t i = 0; i < _total_cells_; i++)
    {
        if(_sdf_[i] == LARGE_NUMBER)
            _layer_[i] = -1;
        else
            _layer_[i] = 0;
    }


    for(size_t j = 0; j < _height_; j++)
    {
        for(size_t i = 0; i < _width_; i++)
        {
            index = index_1D(i, j);
            if(_layer_[index] == 0)
            {
                if(_has_neighbor_in_layer(i, j, -1) == true)
                    _layer_[index] = 1;
            }
        }   
    }
    
    int layer_to_update = 0;
    if(neg == false)
        layer_to_update = -1;

    for(int counter = 1; counter < iterations + 1; ++counter)
    {
        for(size_t j = 0; j < _height_; j++)
        {
            for(size_t i = 0; i < _width_; i++)
            {
                index = index_1D(i, j);
                if(_layer_[index] == layer_to_update)
                {
                    if(_has_neighbor_in_layer(i, j, counter) == true)
                    {
                        //update signed distance at this cell
                        _layer_[index] = counter + 1;
                        if(neg == true)
                            _sdf_[index] = _max_neighbor_sdf_in_layer(i, j, counter)
                                           - _voxel_size_;
                        else
                            _sdf_[index] = _min_neighbor_sdf_in_layer(i, j, counter)
                                           + _voxel_size_;
                    }
                }
            }
        }
    }
}

//Checks if the specified voxel has a neighbor in the desired layer.
bool MacGrid::_has_neighbor_in_layer(size_t i, size_t j, int mlayer) const
{
    if(i > 0)
        if(layer(i - 1, j) == mlayer)
            return true;

    if(i < _width_ - 1)
        if(layer(i + 1, j) == mlayer)
            return true;

    if(j > 0)
        if(layer(i, j - 1) == mlayer)
            return true;

    if(j < _height_ - 1)
        if(layer(i, j + 1) == mlayer)
            return true;

    return false;
}

//Returns the minimum signed distance value in a neighbor to voxel i, j in layer mlayer.
double MacGrid::_min_neighbor_sdf_in_layer(size_t i, size_t j, int mlayer) const
{
    double min = LARGE_NUMBER;
    double d;

    if(i > 0)
    {
        if(layer(i - 1, j) == mlayer)
        {
            d = sdf(i - 1, j);
            if(d < min)
                min = d;
        }
    }

    if(i < _width_ - 1)
    {
        if(layer(i + 1, j) == mlayer)
        {
            d = sdf(i + 1, j);
            if(d < min)
                min = d;
        }
    }

    if(j > 0)
    {
        if(layer(i, j - 1) == mlayer)
        {
            d = sdf(i, j - 1);
            if(d < min)
                min = d;
        }
    }

    if(j < _height_ - 1)
    {
        if(layer(i, j + 1) == mlayer)
        {
            d = sdf(i, j + 1);
            if(d < min)
                min = d;
        }
    }

    return min;
}

//Returns the maximum signed distance value in a neighbor to voxel i, j in layer mlayer.
double MacGrid::_max_neighbor_sdf_in_layer(size_t i, size_t j, int mlayer) const
{
    double max = LARGE_NEG_NUMBER;
    double d;

    if(i > 0)
    {
        if(layer(i - 1, j) == mlayer)
        {
            d = sdf(i - 1, j);
            if(d > max)
                max = d;
        }
    }

    if(i < _width_ - 1)
    {
        if(layer(i + 1, j) == mlayer)
        {
            d = sdf(i + 1, j);
            if(d > max)
                max = d;
        }
    }

    if(j > 0)
    {
        if(layer(i, j - 1) == mlayer)
        {
            d = sdf(i, j - 1);
            if(d > max)
                max = d;
        }
    }

    if(j < _height_ - 1)
    {
        if(layer(i, j + 1) == mlayer)
        {
            d = sdf(i, j + 1);
            if(d > max)
                max = d;
        }
    }

    return max;
}
