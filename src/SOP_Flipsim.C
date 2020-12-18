/*
 * Copyright (c) 2018
 *  Side Effects Software Inc.  All rights reserved.
 *
 * Redistribution and use of Houdini Development Kit samples in source and
 * binary forms, with or without modification, are permitted provided that the
 * following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. The name of Side Effects Software may not be used to endorse or
 *    promote products derived from this software without specific prior
 *    written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SIDE EFFECTS SOFTWARE `AS IS' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN
 * NO EVENT SHALL SIDE EFFECTS SOFTWARE BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *----------------------------------------------------------------------------
 *
 * SOP_Flipsim.C
 *
 *  Created on: September 5, 2018
 *      Author: Sean Flynn
 */

#include <GU/GU_Detail.h>
#include <GU/GU_RayIntersect.h>
#include <GEO/GEO_PrimPart.h>
#include <GEO/GEO_PrimPoly.h>
#include <PRM/PRM_SpareData.h>
#include <OP/OP_AutoLockInputs.h>
#include <OP/OP_Director.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <PRM/PRM_Include.h>
#include <UT/UT_DSOVersion.h>
#include <UT/UT_Interrupt.h>
#include <UT/UT_Vector3.h>
#include <UT/UT_Vector4.h>

#include "SOP_Flipsim.h"

using namespace std;

//constants
const double EPSILON = .00001;
const int LARGE_NUMBER = 100000000;
const size_t MAX_PARTICLES = 1000000000;
const double VECTOR_VIZ_SCALE = 10.0;

//Creates the SOP node in Houdini
void newSopOperator(OP_OperatorTable *table)
{
    table->addOperator(new OP_Operator(
        "hdk_flipsim",
        "flip-fluid-sim",
        SOP_Flipsim::myConstructor,
        SOP_Flipsim::myTemplateList,
        (int)0,      // Min required sources
        (int)1,      // Maximum sources
        0));
}

// The names here have to match the inline evaluation functions
static PRM_Name names[] =
{
    PRM_Name("reset_on_one", "Reset Simulation on Frame 1"),
    PRM_Name("grid_dimensions", "Grid Dimensions"),
    PRM_Name("particle_radius", "Particle Radius"),
    PRM_Name("voxel_size", "Voxel Size"),
    PRM_Name("gravity", "Gravity"),
    PRM_Name("src_vel", "Initial Source Velocity"),
    PRM_Name("st_const", "Surface Tension Constant"),
    PRM_Name("flip_ratio", "FLIP to PIC Ratio"),
    PRM_Name("vis_cell_type", "Cell Type to Visualize"),
    PRM_Name("level_set_thresh", "Level Set Viz Threshold"),
    PRM_Name("show_particles", "Show Fluid Particles"),
    PRM_Name("show_particle_circles", "Show Fluid Particles As Circles"),
    PRM_Name("show_particle_velocities", "Show Fluid Particle Velocities"),
    PRM_Name("show_cells", "Show Grid Cells"),
    PRM_Name("vis_vel_comp", "Visualize Velocity Components"),
    PRM_Name("interp_type", "Visualize Interpolation Type"),
    PRM_Name("interp_particles_sd_bounds", "Min/Max SDF Interpolation Particles Viz"),
    PRM_Name("vis_interp_particles", "Visualize Interpolation Particles"),
    PRM_Name("vis_interp_vectors", "Visualize Interpolation Vectors"),
};

//default UI parameters
static PRM_Default  def_dimensions[] =
{
    PRM_Default(64),
    PRM_Default(32),
    PRM_Default(1)
};
static PRM_Default  def_gravity[] =
{
    PRM_Default(0.0),
    PRM_Default(-9.8 / 24.0),
    PRM_Default(0.0)
};
static PRM_Default  def_src_vel[] = {
   PRM_Default(0.0),
   PRM_Default(0.0),
   PRM_Default(0.0)
};
static PRM_Default  def_st_const(0.2);
static PRM_Default  def_flip_ratio(0.85);
static PRM_Default  def_particle_radius(0.71);
static PRM_Default  def_voxel_size(1.0);
static PRM_Default  def_level_set_thresh(.5);
static PRM_Default  def_sd_interp_bounds[] =
{
    PRM_Default(-1.0),
    PRM_Default(1.0)
};
static PRM_Range def_range(PRM_RANGE_UI, 0.0, 
                           PRM_RANGE_UI, 10.0);
static PRM_Range def_st_range(PRM_RANGE_UI, 0.0, 
                              PRM_RANGE_UI, 3.0);
static PRM_Range def_flip_ratio_range(PRM_RANGE_UI, 0.0, 
                                      PRM_RANGE_UI, 1.0);
static PRM_Name cell_types[] =
{
    PRM_Name("all", "All"),
    PRM_Name("ls", "Level Set"),
    PRM_Name("fluid", "Fluid"),
    PRM_Name("air", "Air"),
    PRM_Name("fluid_air", "Fluid and Air"),
    PRM_Name("solid", "Solid"),
    PRM_Name("unused", "Unused"),
    PRM_Name(0)
};
static PRM_ChoiceList cell_types_menu(PRM_CHOICELIST_SINGLE, cell_types);
static PRM_Name interp_names[] =
{
    PRM_Name("ux", "UX"),
    PRM_Name("uy", "UY"),
    PRM_Name("uz", "UZ"),
    PRM_Name("u", "U"),
    PRM_Name("sd", "SD"),
    PRM_Name("p", "P"),
    PRM_Name(0)
};
static PRM_ChoiceList interp_menu(PRM_CHOICELIST_SINGLE, interp_names);

//list of templates for UI parameters
PRM_Template SOP_Flipsim::myTemplateList[] =
{
    PRM_Template(PRM_TOGGLE, 1, &names[0]),
    PRM_Template(PRM_INT, 2, &names[1], def_dimensions),
    PRM_Template(PRM_FLT, 1, &names[2], &def_particle_radius, 0, &def_range),
    PRM_Template(PRM_FLT, 1, &names[3], &def_voxel_size, 0, &def_range),
    PRM_Template(PRM_FLT, 2, &names[4], def_gravity),
    PRM_Template(PRM_FLT, 2, &names[5], def_src_vel),
    PRM_Template(PRM_FLT, 1, &names[6], &def_st_const, 0, &def_st_range),
    PRM_Template(PRM_FLT, 1, &names[7], &def_flip_ratio, 0, &def_flip_ratio_range),
    PRM_Template(PRM_ORD, PRM_Template::PRM_EXPORT_MAX, 1,
                                &names[8], 0, &cell_types_menu),
    PRM_Template(PRM_FLT, 1, &names[9], &def_level_set_thresh, 0, &def_range),
    PRM_Template(PRM_TOGGLE, 1, &names[10], PRMoneDefaults),
    PRM_Template(PRM_TOGGLE, 1, &names[11]),
    PRM_Template(PRM_TOGGLE, 1, &names[12]),
    PRM_Template(PRM_TOGGLE, 1, &names[13]),
    PRM_Template(PRM_TOGGLE, 1, &names[14]),
    PRM_Template(PRM_ORD, PRM_Template::PRM_EXPORT_MAX, 1,
                                &names[15], 0, &interp_menu),
    PRM_Template(PRM_FLT, 2, &names[16], def_sd_interp_bounds),
    PRM_Template(PRM_TOGGLE, 1, &names[17]),
    PRM_Template(PRM_TOGGLE, 1, &names[18]),
    PRM_Template(),
};

int *SOP_Flipsim::myOffsets = 0;

//Calls the SOP constructor
OP_Node * SOP_Flipsim::myConstructor(OP_Network *net, const char *name, OP_Operator *op)
{
    return new SOP_Flipsim(net, name, op);
}

//Constructs a SOP_Flipsim object.
SOP_Flipsim::SOP_Flipsim(OP_Network *net, const char *name, OP_Operator *op)
    : SOP_Node(net, name, op)
    , _particle_system_(nullptr), _simulator_(nullptr), _poly_pt_group_(nullptr)
{
    _prev_source_pos_[0] = 0.0;
    _prev_source_pos_[1] = 0.0;
    _prev_source_pos_[2] = 0.0;

    // Make sure that our offsets are allocated.  Here we allow up to 32
    // parameters, no harm in over allocating.  The definition for this
    // function is in OP/OP_Parameters.h
    if (!myOffsets)
        myOffsets = allocIndirect(32);
}

//Creates the label for the inputs to this node
const char * SOP_Flipsim::inputLabel(unsigned idx) const
{
    switch (idx)
    {
        case 0: return "Source fluid particles";
        default: return "default";
    }
}

//Deletes this SOP_Flipsim and any associated data.
SOP_Flipsim::~SOP_Flipsim()
{
    delete _simulator_;
}

//Initializes this SOP_Flipsim.
void SOP_Flipsim::initSystem()
{
    if(!gdp)
       gdp = new GU_Detail;

    // Check to see if we really need to reset everything
    if(!_simulator_ || !_particle_system_ || _particle_system_->getNumParticles() > 0 ||
       gdp->getPointMap().indexSize() > 0)
    {
        gdp->clearAndDestroy();
        _particle_system_ = (GEO_PrimParticle *)gdp->appendPrimitive(GEO_PRIMPART);
        _particle_system_->clearAndDestroy();

        _poly_pt_group_ = gdp->newPointGroup("__test_group__", false);

        if(!_simulator_)
            _reset_simulator();
    }
}

//Cooks the geometry for this SOP.
OP_ERROR SOP_Flipsim::cookMySop(OP_Context &context)
{
    fpreal t = context.getTime();
    size_t frame = (size_t)(t * 24.0) + 1;
    OP_Node::flags().setTimeDep(true);

    initSystem();
    
    if(frame == 1 && RESETONONE())
        _reset_simulator();

    //get parameters from ui and set them in the simulator
    _simulator_->set_st_const(STCONST());
    UT_Vector3 g(XGRAV(), YGRAV(), ZGRAV());
    _simulator_->set_gravity(g);
    _simulator_->set_flip_ratio(FLIPRATIO());

    OP_AutoLockInputs inputs(this);
    if (inputs.lock(context) >= UT_ERROR_ABORT)
        return error();

    //get the source particles if any
    const GU_Detail *points = inputGeo(0);
    if(points)
    {
        GA_ROHandleV3 Phandle(points->findAttribute(GA_ATTRIB_POINT, "P"));
        GA_Offset ptoff;
        UT_Vector3 pos, v;
        v[0] = XSRCVEL();
        v[1] = YSRCVEL();
        v[2] = ZSRCVEL();

        GA_FOR_ALL_PTOFF(points, ptoff)
        {
            // get current pont position
            pos = Phandle.get(ptoff);
            _simulator_->add_particle(pos[0], pos[1], .5, v[0], v[1], v[2], frame - 1);
        }
    }

    _simulator_->simulate_flip_to_frame(frame);

    if(SHOWPARTICLES())
        _draw_fluid_particles(frame);

    if(VISINTERP())
    {
        int type_index = INTERPTYPE();
        ScalarType type;
        if(type_index == 0)
            type = UX;
        else if(type_index == 1)
            type = UY;
        else if(type_index == 2)
            type = UZ;
        else if(type_index == 3)
            type = U;
        else if(type_index == 4)
            type = SD;
        else if(type_index == 5)
            type = PRESSURE;

        _draw_interp_particles(CELLTYPE(), type, frame);
    }

    _delete_polygons();

    if(SHOWPARTICLES())
    {
        if(SHOWPARTICLECIRCLES())
            _draw_fluid_circles(frame);
        if(SHOWPARTICLEVELOCITIES())
            _draw_fluid_velocities(frame);
    }

    if(VISINTERPVEC())
    {
        int type_index = INTERPTYPE();
        ScalarType type;
        if(type_index == 0)
            type = UX;
        else if(type_index == 1)
            type = UY;
        else if(type_index == 2)
            type = UZ;
        else if(type_index == 3)
            type = U;
        else if(type_index == 4)
            type = SD;

        _draw_interp_vectors(CELLTYPE(), type, frame);
    }

    _draw_viz(frame);
    return error();
}

//Resets the simulator that this SOP_Flipsim visualizes.
void SOP_Flipsim::_reset_simulator(void)
{
    int x_dim = XDIM();
    int y_dim = YDIM();
    int z_dim = ZDIM();
    fpreal p_rad = PRADIUS();
    fpreal v_size = VSIZE();

    if(_simulator_)
        delete _simulator_;

    MacGrid *grid = new MacGrid(x_dim, y_dim, z_dim, v_size, p_rad);

    double st_const = STCONST();
    UT_Vector3 gravity(XGRAV(), YGRAV(), ZGRAV());
    double flip_ratio = FLIPRATIO();

    _simulator_ = new Simulator(grid, gravity, st_const, flip_ratio);

    _ip_count_ = grid->total_cells() * 84;
    for(size_t i = 0; i < _ip_count_; ++i)
    {
        _interp_particles_[i][0] = fRand(0.0, x_dim);
        _interp_particles_[i][1] = fRand(0.0, y_dim);
        _interp_particles_[i][2] = 0.0;
    }
}

//Draws the fluid particles being simulated.
void SOP_Flipsim::_draw_fluid_particles(const size_t frame)
{
    vector<UT_Vector3> particle_positions = _simulator_->get_particle_positions(frame);
    GA_RWHandleV3 colorh(gdp->findDiffuseAttribute(GA_ATTRIB_POINT));
    if(particle_positions.size() > 0)
    {
        if (!colorh.isValid())
            colorh = GA_RWHandleV3(gdp->addDiffuseAttribute(GA_ATTRIB_POINT));
    }

    UT_Vector3 pos;
    for(size_t i = 0; i < particle_positions.size(); i++)
    {
        pos = particle_positions[i];
        GA_Offset vtxoff = _particle_system_->giveBirth();
        UT_Vector3 clr(.05, .1, .95);
        colorh.set(vtxoff, clr);
        _particle_system_->setPos3(vtxoff, pos);
    }
}

//Draws circles around fluid particles.
void SOP_Flipsim::_draw_fluid_circles(const size_t frame)
{
    vector<UT_Vector3> particle_positions = _simulator_->get_particle_positions(frame);
    GA_RWHandleV3 colorh(gdp->findDiffuseAttribute(GA_ATTRIB_POINT));
    if(particle_positions.size() > 0)
    {
        if (!colorh.isValid())
            colorh = GA_RWHandleV3(gdp->addDiffuseAttribute(GA_ATTRIB_POINT));
    }

    UT_Vector3 pos;
    double radius = _simulator_->get_grid(frame)->particle_radius();

    for(size_t i = 0; i < particle_positions.size(); i++)
    {
        pos = particle_positions[i];
        _draw_circle(pos[0], pos[1], pos[2] - .02, radius, .09, .15, .75);
    }
}

//Draws fluid velocity vectors.
void SOP_Flipsim::_draw_fluid_velocities(const size_t frame)
{
    vector<UT_Vector3> particle_positions = _simulator_->get_particle_positions(frame);
    vector<UT_Vector3> particle_velocities = _simulator_->get_particle_velocities(frame);
    GA_RWHandleV3 colorh(gdp->findDiffuseAttribute(GA_ATTRIB_POINT));
    if(particle_positions.size() > 0)
    {
        if (!colorh.isValid())
            colorh = GA_RWHandleV3(gdp->addDiffuseAttribute(GA_ATTRIB_POINT));
    }

    UT_Vector3 pos, vel;
    for(size_t i = 0; i < particle_positions.size(); i++)
    {
        pos = particle_positions[i];
        vel = particle_velocities[i];
        _draw_vector(pos[0], pos[1], pos[2], vel[0] * 3, vel[1] * 3, vel[2] * 3);
    }
}

//Draws particles that show interpolated values on the MacGrid.
void SOP_Flipsim::_draw_interp_particles(const int cell_type, const ScalarType type,
                                         const size_t frame)
{
    GA_RWHandleV3 colorh(gdp->findDiffuseAttribute(GA_ATTRIB_POINT));
    if(_ip_count_ > 0)
    {
        if (!colorh.isValid())
            colorh = GA_RWHandleV3(gdp->addDiffuseAttribute(GA_ATTRIB_POINT));
    }
    else
        return;

    UT_Vector3 pos, clr, u;
    double sd, p;
    int index;
    double level_set_thresh = LEVELSETTHRESH();
    double min_sd = INTERPMINSD();
    double max_sd = INTERPMAXSD();
    MacGrid *grid = _simulator_->get_grid(frame);

    for(size_t i = 0; i < _ip_count_; i++)
    {
        pos = _interp_particles_[i];
        GA_Offset vtxoff = _particle_system_->giveBirth();
        sd = grid->get_sd(pos[0], pos[1]);
        if(sd < min_sd || sd > max_sd)
            continue;

        index = grid->index_1D((size_t)floor(pos[0]), (size_t)floor(pos[1]));
        if(index == -1)
            continue;

        if(cell_type == 1 && abs(grid->sdf(index)) > level_set_thresh)
            continue;
        else if(cell_type == 2 && grid->type(index) != FLUID)
            continue;
        else if(cell_type == 3 && grid->type(index) != AIR)
            continue;
        else if(cell_type == 4 && (grid->type(index) != AIR && grid->type(index) != FLUID))
                continue;
        else if(cell_type == 5 && grid->type(index) != SOLID)
            continue;
        else if(cell_type == 6 && grid->type(index) != UNUSED)
            continue;

        if(type == UX || type == UY || type == UZ || type == U)
        {
            grid->get_velocity(pos[0], pos[1], pos[2], CUR, u);

            clr[0] = 0;
            clr[1] = 0;
            clr[2] = 0;

            if(type == UX || type == U)
                clr[0] = abs(u[0]) / 2.0;
            if(type == UY || type == U)
                clr[1] = abs(u[1]) / 2.0;
            if(type == UZ || type == U)
                clr[2] = abs(u[2]) / 2.0;
        }
        else if(type == SD)
        {
            sd = grid->get_sd(pos[0], pos[1]);
            if(sd >= 0)
            {
                clr[0] = sd;
                clr[1] = sd;
                clr[2] = sd;
            }
            else
            {
                clr[0] = 0;
                clr[1] = 0;
                clr[2] = abs(sd);
            }
        }
        else if(type == PRESSURE)
        {
            p = grid->get_pressure(pos[0], pos[1]);
            if(p > 0)
            {
                //clr[0] = abs(p);
                clr[0] = 1.0;
                clr[1] = 0.0;
                clr[2] = 0.0;
            }
            else
            {
                clr[0] = abs(p);
                clr[1] = abs(p);
                clr[2] = abs(p);
            }
        }
        colorh.set(vtxoff, clr);
        _particle_system_->setPos3(vtxoff, pos);
    }
}

//Draws vectors that show interpolated values on the MacGrid.
void SOP_Flipsim::_draw_interp_vectors(const int cell_type, const ScalarType type,
                                       const size_t frame)
{
    UT_Vector3 pos, clr, u;
    double sd;
    int index;
    double level_set_thresh = LEVELSETTHRESH();
    double min_sd = INTERPMINSD();
    double max_sd = INTERPMAXSD();
    MacGrid *grid = _simulator_->get_grid(frame);
    
    for(size_t i = 0; i < _ip_count_; i++)
    {
        if(i % 4)
            continue;

        if(type == UX || type == UY || type == UZ || type == U)
        {
            pos = _interp_particles_[i];
            sd = grid->get_sd(pos[0], pos[1]);
            if(sd < min_sd || sd > max_sd)
                continue;

            index = grid->index_1D((size_t)floor(pos[0]), (size_t)floor(pos[1]));
            if(index == -1)
                continue;

            if(cell_type == 1 && abs(grid->sdf(index)) > level_set_thresh)
                continue;
            else if(cell_type == 2 && grid->type(index) != FLUID)
                continue;
            else if(cell_type == 3 && grid->type(index) != AIR)
                continue;
            else if(cell_type == 4 && (grid->type(index) != AIR && grid->type(index) != FLUID))
                continue;
            else if(cell_type == 5 && grid->type(index) != SOLID)
                continue;
            else if(cell_type == 6 && grid->type(index) != UNUSED)
                continue;

            grid->get_velocity(pos[0], pos[1], pos[2], CUR, u);

            if(type == UX)
                _draw_vector(pos[0], pos[1], pos[2], u[0] / 2.0, 0.0, 0.0);
            else if(type == UY)
                _draw_vector(pos[0], pos[1], pos[2], 0.0, u[1] / 2.0, 0.0);
            else if(type == UZ)
                _draw_vector(pos[0], pos[1], pos[2], 0.0, 0.0, u[2] / 2.0);
            else
                _draw_vector(pos[0], pos[1], pos[2], u[0] * 2.0, u[1] * 2.0, u[2] * 2.0);
        }
    }
}

//Draws the visualizations for this FLIP simulation.
void SOP_Flipsim::_draw_viz(const size_t frame)
{
    int show_cells = SHOWGRIDCELLS();
    if(!show_cells)
        return;

    int cell_type = CELLTYPE();
    int vis_vel_comp = VISVELCOMP();

    double level_set_thresh = LEVELSETTHRESH();
    double ux, uy;
    MacGrid *grid = _simulator_->get_grid(frame);

    for(size_t j = 0; j < grid->height(); j++)
    {
        for(size_t i = 0; i < grid->width(); i++)
        {
            if(cell_type == 0 || cell_type == 2 || cell_type == 4)//Fluid
            {
                if(grid->type(i, j) == FLUID)
                {
                    _draw_cube(i, j, 0, 1, 20, 20, 160, false);

                    if(vis_vel_comp)
                    {
                        ux = grid->ux(i, j);
                        uy = grid->uy(i, j);

                        _draw_vector(i, j + .5, 0.05, ux * VECTOR_VIZ_SCALE, 0, 0);
                        _draw_vector(i + .5, j, 0.05, 0, uy * VECTOR_VIZ_SCALE, 0);
                    }
                }
            }
            if(cell_type == 0 || cell_type == 3 || cell_type == 4)//Air
            {
                if(grid->type(i, j) == AIR)
                {
                    _draw_cube(i, j, 0, 1, 20, 20, 160, false);

                    if(vis_vel_comp)
                    {
                        ux = grid->ux(i, j);
                        uy = grid->uy(i, j);

                        _draw_vector(i, j + .5, 0.05, ux * VECTOR_VIZ_SCALE, 0, 0);
                        _draw_vector(i + .5, j, 0.05, 0, uy * VECTOR_VIZ_SCALE, 0);
                    }
                }
            }
            if(cell_type == 0 || cell_type == 5)//Solid
            {
                if(grid->type(i, j) == SOLID)
                {
                    _draw_cube(i, j, 0, 1, 240, 40, 40, false);

                    if(vis_vel_comp)
                    {
                        ux = grid->ux(i, j);
                        uy = grid->uy(i, j);

                        _draw_vector(i, j + .5, 0.05, ux * VECTOR_VIZ_SCALE, 0, 0);
                        _draw_vector(i + .5, j, 0.05, 0, uy * VECTOR_VIZ_SCALE, 0);
                    }
                }
            }
            if(cell_type == 0 || cell_type == 6)//Unused
            {
                if(grid->type(i, j) == UNUSED)
                {
                    _draw_cube(i, j, 0, 1, 10, 10, 10, false);

                    if(vis_vel_comp)
                    {
                        ux = grid->ux(i, j);
                        uy = grid->uy(i, j);

                        _draw_vector(i, j + .5, 0.05, ux * VECTOR_VIZ_SCALE, 0, 0);
                        _draw_vector(i + .5, j, 0.05, 0, uy * VECTOR_VIZ_SCALE, 0);
                    }
                }
            }
            if(cell_type == 0 || cell_type == 1)//Level set
            {
                if(abs(grid->sdf(i, j)) <= level_set_thresh)
                {
                    _draw_cube(i, j, 0, 1, 10, 10, 100, false);

                    if(vis_vel_comp)
                    {
                        ux = grid->ux(i, j);
                        uy = grid->uy(i, j);

                        _draw_vector(i, j + .5, 0.05, ux * VECTOR_VIZ_SCALE, 0, 0);
                        _draw_vector(i + .5, j, 0.05, 0, uy * VECTOR_VIZ_SCALE, 0);
                    }
                }
            }
        }
    }
}

//Draws a circle with the given parameters.
GEO_PrimPoly * SOP_Flipsim::_draw_circle(const fpreal x, const fpreal y, const fpreal z,
                                         const fpreal radius, const fpreal r, const fpreal g,
                                         const fpreal b)
{
    GA_RWHandleV3 colorh(gdp->findDiffuseAttribute(GA_ATTRIB_POINT));
    if (!colorh.isValid())
        colorh = GA_RWHandleV3(gdp->addDiffuseAttribute(GA_ATTRIB_POINT));

    //a velocity magnitude of 1 would be fully colored
    UT_Vector3 clr(r, g, b);

    //TODO--replace this so it's not hard coded.
    double sqrt_3_2 = sqrt(3.0) / 2.0;
    double sqrt_2_2 = sqrt(2.0) / 2.0;
    double x_pos[17] = {1.0, sqrt_3_2, sqrt_2_2, .5, 0.0, -.5, -1.0 * sqrt_2_2, -1.0 * sqrt_3_2,
                        -1.0, -1.0 * sqrt_3_2, -1.0 * sqrt_2_2, -.5, 0.0, .5, sqrt_2_2, sqrt_3_2,
                        1.0};
    double y_pos[17] = {0.0, .5, sqrt_2_2, sqrt_3_2, 1.0, sqrt_3_2, sqrt_2_2, .5, 0.0,
                        -.5, -1.0 * sqrt_2_2, -1.0 * sqrt_3_2, -1.0, -1.0 * sqrt_3_2,
                        -1.0 * sqrt_2_2, -.5, 0.0};

    GEO_PrimPoly *poly = GEO_PrimPoly::build(gdp, 17, GU_POLY_CLOSED);
    UT_Vector3 pos;
    int index = 0;
    for(size_t i = 0; i < 17; i++)
    {
        pos[0] = x + x_pos[i] * radius;
        pos[1] = y + y_pos[i] * radius;
        pos[2] = z;
        GA_Offset pptoff = poly->getPointOffset(index++);
        gdp->setPos3(pptoff, pos);
        colorh.set(pptoff, clr);
        _poly_pt_group_->addOffset(pptoff);
    }

    return poly;
}

//Draws a vector with the given parameters.
GEO_PrimPoly * SOP_Flipsim::_draw_vector(const fpreal px, const fpreal py, const fpreal pz,
                                         const fpreal vx, const fpreal vy, const fpreal vz)
{
    GA_RWHandleV3 colorh(gdp->findDiffuseAttribute(GA_ATTRIB_POINT));
    if (!colorh.isValid())
        colorh = GA_RWHandleV3(gdp->addDiffuseAttribute(GA_ATTRIB_POINT));

    //a velocity magnitude of 1 would be fully colored
    UT_Vector3 clr(abs(vx), abs(vy), abs(vz));

    GEO_PrimPoly *poly = GEO_PrimPoly::build(gdp, 5, GU_POLY_OPEN);
    UT_Vector3 pos;
    int index = 0;

    pos[0] = px;
    pos[1] = py;
    pos[2] = pz;
    GA_Offset pptoff = poly->getPointOffset(index++);
    gdp->setPos3(pptoff, pos);
    colorh.set(pptoff, clr);
    _poly_pt_group_->addOffset(pptoff);

    pos[0] = px + vx;
    pos[1] = py + vy;
    pos[2] = pz + vz;
    pptoff = poly->getPointOffset(index++);
    gdp->setPos3(pptoff, pos);
    colorh.set(pptoff, clr);
    _poly_pt_group_->addOffset(pptoff);


    pos[0] = (px + (vx * .9)) + (-1 * vy / 10.0);
    pos[1] = (py + (vy * .9)) + (vx / 10.0);
    pos[2] = pz;
    pptoff = poly->getPointOffset(index++);
    gdp->setPos3(pptoff, pos);
    colorh.set(pptoff, clr);
    _poly_pt_group_->addOffset(pptoff);

    pos[0] = px + vx;
    pos[1] = py + vy;
    pos[2] = pz + vz;
    pptoff = poly->getPointOffset(index++);
    gdp->setPos3(pptoff, pos);
    colorh.set(pptoff, clr);
    _poly_pt_group_->addOffset(pptoff);

    pos[0] = (px + (vx * .9)) - (-1 * vy / 10.0);
    pos[1] = (py + (vy * .9)) - (vx / 10.0);
    pos[2] = pz;
    pptoff = poly->getPointOffset(index++);
    gdp->setPos3(pptoff, pos);
    colorh.set(pptoff, clr);
    _poly_pt_group_->addOffset(pptoff);

    return poly;
}

//Draws a cube with the given parameters.
GEO_PrimPoly * SOP_Flipsim::_draw_cube(const fpreal x, const fpreal y, const fpreal z,
                                       const fpreal size, const fpreal r, const fpreal g,
                                       const fpreal b, const bool filled)
{
    return _draw_box(x, y, z, size, size, size, r, g, b, filled);
}

//Draws a box with the given parameters.
GEO_PrimPoly * SOP_Flipsim::_draw_box(const fpreal x, const fpreal y, const fpreal z,
                                      const fpreal sx, const fpreal sy, const fpreal sz,
                                      const fpreal r, const fpreal g, const fpreal b,
                                      const bool filled)
{
    GA_RWHandleV3 colorh(gdp->findDiffuseAttribute(GA_ATTRIB_POINT));
    if (!colorh.isValid())
        colorh = GA_RWHandleV3(gdp->addDiffuseAttribute(GA_ATTRIB_POINT));

    UT_Vector3 clr(r, g, b);
    clr *= 1.0/255.0;

    GEO_PrimPoly *poly;
    if(filled == true)
        poly = GEO_PrimPoly::build(gdp, 16, GU_POLY_CLOSED);
    else
        poly = GEO_PrimPoly::build(gdp, 16, GU_POLY_OPEN);

    UT_Vector3 pos;
    int index = 0;

    pos[0] = x;
    pos[1] = y;
    pos[2] = z;
    GA_Offset pptoff = poly->getPointOffset(index++);
    gdp->setPos3(pptoff, pos);
    colorh.set(pptoff, clr);
    _poly_pt_group_->addOffset(pptoff);

    pos[0] = x + sx;
    pos[1] = y;
    pos[2] = z;
    pptoff = poly->getPointOffset(index++);
    gdp->setPos3(pptoff, pos);
    colorh.set(pptoff, clr);
    _poly_pt_group_->addOffset(pptoff);

    pos[0] = x + sx;
    pos[1] = y + sy;
    pos[2] = z;
    pptoff = poly->getPointOffset(index++);
    gdp->setPos3(pptoff, pos);
    colorh.set(pptoff, clr);
    _poly_pt_group_->addOffset(pptoff);

    pos[0] = x;
    pos[1] = y + sy;
    pos[2] = z;
    pptoff = poly->getPointOffset(index++);
    gdp->setPos3(pptoff, pos);
    colorh.set(pptoff, clr);
    _poly_pt_group_->addOffset(pptoff);

    pos[0] = x;
    pos[1] = y;
    pos[2] = z;
    pptoff = poly->getPointOffset(index++);
    gdp->setPos3(pptoff, pos);
    colorh.set(pptoff, clr);
    _poly_pt_group_->addOffset(pptoff);

    pos[0] = x;
    pos[1] = y;
    pos[2] = z + sz;
    pptoff = poly->getPointOffset(index++);
    gdp->setPos3(pptoff, pos);
    colorh.set(pptoff, clr);
    _poly_pt_group_->addOffset(pptoff);

    pos[0] = x + sx;
    pos[1] = y;
    pos[2] = z + sz;
    pptoff = poly->getPointOffset(index++);
    gdp->setPos3(pptoff, pos);
    colorh.set(pptoff, clr);
    _poly_pt_group_->addOffset(pptoff);

    pos[0] = x + sx;
    pos[1] = y;
    pos[2] = z;
    pptoff = poly->getPointOffset(index++);
    gdp->setPos3(pptoff, pos);
    colorh.set(pptoff, clr);
    _poly_pt_group_->addOffset(pptoff);

    pos[0] = x + sx;
    pos[1] = y;
    pos[2] = z + sz;
    pptoff = poly->getPointOffset(index++);
    gdp->setPos3(pptoff, pos);
    colorh.set(pptoff, clr);
    _poly_pt_group_->addOffset(pptoff);

    pos[0] = x + sx;
    pos[1] = y + sy;
    pos[2] = z + sz;
    pptoff = poly->getPointOffset(index++);
    gdp->setPos3(pptoff, pos);
    colorh.set(pptoff, clr);
    _poly_pt_group_->addOffset(pptoff);

    pos[0] = x + sx;
    pos[1] = y + sy;
    pos[2] = z;
    pptoff = poly->getPointOffset(index++);
    gdp->setPos3(pptoff, pos);
    colorh.set(pptoff, clr);
    _poly_pt_group_->addOffset(pptoff);

    pos[0] = x + sx;
    pos[1] = y + sy;
    pos[2] = z + sz;
    pptoff = poly->getPointOffset(index++);
    gdp->setPos3(pptoff, pos);
    colorh.set(pptoff, clr);
    _poly_pt_group_->addOffset(pptoff);

    pos[0] = x;
    pos[1] = y + sy;
    pos[2] = z + sz;
    pptoff = poly->getPointOffset(index++);
    gdp->setPos3(pptoff, pos);
    colorh.set(pptoff, clr);
    _poly_pt_group_->addOffset(pptoff);

    pos[0] = x;
    pos[1] = y + sy;
    pos[2] = z;
    pptoff = poly->getPointOffset(index++);
    gdp->setPos3(pptoff, pos);
    colorh.set(pptoff, clr);
    _poly_pt_group_->addOffset(pptoff);

    pos[0] = x;
    pos[1] = y + sy;
    pos[2] = z + sz;
    pptoff = poly->getPointOffset(index++);
    gdp->setPos3(pptoff, pos);
    colorh.set(pptoff, clr);
    _poly_pt_group_->addOffset(pptoff);

    pos[0] = x;
    pos[1] = y;
    pos[2] = z + sz;
    pptoff = poly->getPointOffset(index++);
    gdp->setPos3(pptoff, pos);
    colorh.set(pptoff, clr);
    _poly_pt_group_->addOffset(pptoff);

    return poly;
}

//Deletes the polygons in this SOP_Node
void SOP_Flipsim::_delete_polygons(void)
{
    GA_Primitive* prim;
    GA_Size primCount = gdp->primitives().entries();

    for(GA_Offset i = 1; i < primCount; i++)
    {
        // Get the first primitive in the detail (this assumes it exists).
        prim = gdp->getPrimitiveList().get(i);
        gdp->destroyPrimitive(*prim);
    }

    gdp->destroyUnusedPoints();
    gdp->destroyDegeneratePrimitives();
}
