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
 * SOP_Flipsim.h
 *
 *  Created on: September 5, 2018
 *      Author: Sean Flynn
 */

#ifndef __SOP_Flipsim_h__
#define __SOP_Flipsim_h__

#include <memory>

#include <SOP/SOP_Node.h>

#include "macgrid.h"
#include "simulator.h"

class GEO_ParticleVertex;
class GEO_PrimParticle;
class GU_RayIntersect;

/**
* This class defines a FLIP fluid simulator SOP node in Houdini. It is responsible for the UI
* and all the visualizations for the Simulator and MacGrid classes. 
*/
class SOP_Flipsim : public SOP_Node
{
public:

    /**
    * Constructs a SOP_Flipsim object.
    *
    * @param net Needed by HDK to initialize SOP_Node object
    * @param name Needed by HDK to initialize SOP_Node object
    * @param op Needed by HDK to initialize SOP_Node object
    */
    SOP_Flipsim(OP_Network *net, const char *name, OP_Operator *op);

    //code provided by SOP_Node sample code
    static PRM_Template myTemplateList[];
    static OP_Node *myConstructor(OP_Network*, const char *, OP_Operator *);

protected:

    /**
    * Creates the label for the inputs to this node
    *
    * @param idx The index of the input
    */
    virtual const char *inputLabel(unsigned idx) const;

    /**
    * Initializes this SOP_Flipsim.
    */
    void initSystem();

    /**
    * Cooks the geometry for this SOP.
    *
    * @param context The current context of this SOP_Node
    */
    virtual OP_ERROR cookMySop(OP_Context &context);

private:

    //variable from sample SOP_Node
    static int *myOffsets;


    //variables needed to create particles and polygons
    GEO_PrimParticle *_particle_system_;
    std::unique_ptr<Simulator> _simulator_;
    UT_Vector3 _interp_particles_[10000000];
    size_t _ip_count_;
    UT_Vector3 _prev_source_pos_;
    GA_PointGroup *_poly_pt_group_;

    /**
    * Functions to get the values from UI elements
    */
    int RESETONONE() { return evalInt  ("reset_on_one", 0, 0); }
    int XDIM()   { return evalInt("grid_dimensions", 0, 0); }
    int YDIM()   { return evalInt("grid_dimensions", 1, 0); }
    int ZDIM()   { return 1; }
    fpreal PRADIUS() { return evalFloat("particle_radius", 0, 0); }
    fpreal VSIZE() { return evalFloat("voxel_size", 0, 0); }
    fpreal XGRAV()   { return evalFloat("gravity", 0, 0); }
    fpreal YGRAV()   { return evalFloat("gravity", 1, 0); }
    fpreal ZGRAV()   { return 0.0; }
    fpreal XSRCVEL()   { return evalFloat("src_vel", 0, 0); }
    fpreal YSRCVEL()   { return evalFloat("src_vel", 1, 0); }
    fpreal ZSRCVEL()   { return 0.0; }
    fpreal STCONST() { return evalFloat("st_const", 0, 0); }
    fpreal FLIPRATIO() { return evalFloat("flip_ratio", 0, 0); }
    int CELLTYPE() { return evalInt("vis_cell_type", 0, 0); }
    fpreal LEVELSETTHRESH() { return evalFloat("level_set_thresh", 0, 0); }
    int SHOWPARTICLES() { return evalInt  ("show_particles", 0, 0); }
    int SHOWPARTICLECIRCLES() { return evalInt  ("show_particle_circles", 0, 0); }
    int SHOWPARTICLEVELOCITIES() { return evalInt  ("show_particle_velocities", 0, 0); }
    int SHOWGRIDCELLS() { return evalInt  ("show_cells", 0, 0); }
    int VISVELCOMP() { return evalInt  ("vis_vel_comp", 0, 0); }
    int INTERPTYPE() { return evalInt("interp_type", 0, 0); }
    fpreal INTERPMINSD()   { return evalFloat("interp_particles_sd_bounds", 0, 0); }
    fpreal INTERPMAXSD()   { return evalFloat("interp_particles_sd_bounds", 1, 0); }
    int VISINTERP() { return evalInt  ("vis_interp_particles", 0, 0); }
    int VISINTERPVEC() { return evalInt  ("vis_interp_vectors", 0, 0); }

    /**
    * Resets the simulator that this SOP_Flipsim visualizes.
    */
    void _reset_simulator(void);

    /**
    * Draws the fluid particles being simulated.
    */
    void _draw_fluid_particles(size_t frame);

    /**
    * Draws circles around fluid particles.
    */
    void _draw_fluid_circles(size_t frame);

    /**
    * Draws fluid velocity vectors.
    */
     void _draw_fluid_velocities(size_t frame);

    /**
    * Draws particles that show interpolated values on the MacGrid.
    */
    void _draw_interp_particles(int cell_type, ScalarType type, size_t frame);

    /**
    * Draws vectors that show interpolated values on the MacGrid.
    */
    void _draw_interp_vectors(int cell_type, ScalarType type, size_t frame);

    /**
    * Draws the visualizations for this FLIP simulation
    */
    void _draw_viz(size_t frame);

    /**
    * Draws a circle with the given parameters.
    */
    GEO_PrimPoly * _draw_circle(fpreal x, fpreal y, fpreal z,
                                fpreal radius, fpreal r, fpreal g,
                                fpreal b);

    /**
    * Draws a cube with the given parameters.
    */
    GEO_PrimPoly * _draw_cube(fpreal x, fpreal y, fpreal z, fpreal size,
                              fpreal r, fpreal g, fpreal b, bool filled);

    /**
    * Draws a box with the given parameters.
    */
    GEO_PrimPoly * _draw_box(fpreal x, fpreal y, fpreal z,
                             fpreal sx, fpreal sy, fpreal sz,
                             fpreal r, fpreal g, fpreal b, bool filled);

    /**
    * Draws a vector with the given parameters.
    */
    GEO_PrimPoly * _draw_vector(fpreal px, fpreal py, fpreal pz,
                                fpreal vx, fpreal vy, fpreal vz);

    /**
    * Deletes the polygons in this SOP_Node
    */
    void _delete_polygons(void);
};

#endif
