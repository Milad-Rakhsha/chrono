// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Milad Rakhsha
// =============================================================================

#include <cstdio>
#include <vector>
#include <cmath>
#include <memory>
#include "chrono_parallel/physics/ChSystemParallel.h"

#include "chrono/ChConfig.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsInputOutput.h"
#include "chrono/utils/ChUtilsGeometry.h"
#include "chrono/utils/ChUtilsGenerators.h"
#include "chrono_thirdparty/filesystem/path.h"
#include "chrono_thirdparty/filesystem/resolver.h"
#include "Utils.h"

std::string out_folder = "Box_Spheres/";
#ifdef CHRONO_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

using namespace chrono;
using namespace chrono::collision;
std::shared_ptr<ChBody> box;

double box_mass = 1.0;
double friction = 0.1f;
ChVector<> hdim(2, 2, 0.1);

double sphere_radius = 0.1;

// -----------------------------------------------------------------------------
// Create the falling spherical objects in a uniform rectangular grid.
// -----------------------------------------------------------------------------
void CreateModel(ChSystemParallel* sys) {
    double envelope = sphere_radius * 0.05;
    double sphere_mass = 1.0;
    double sphere_radius = 0.1;
    double sphere_z = -sphere_radius;
    double box_z = sphere_z + sphere_radius + hdim.z() + envelope;

    // Create the middle ball material
    // Material properties (same on bin and balls)
    float Y = 1e8;
    float cr = 0.0f;
    auto mat = chrono_types::make_shared<ChMaterialSurfaceSMC>();
    mat->SetFriction(friction);

    mat->SetYoungModulus(Y);
    mat->SetPoissonRatio(0.4);
    //    sys->GetSettings()->solver.use_material_properties = false;
    //    mat->SetKn(1e10);
    //    mat->SetGn(10);
    //    mat->SetKt(1e10);
    //    mat->SetGt(10);

    //    sys->GetSettings()->solver.tangential_displ_mode == ChSystemSMC::TangentialDisplacementModel::MultiStep;
    mat->SetRestitution(cr);
    mat->SetAdhesion(0);  // Magnitude of the adhesion in Constant adhesion model
    int ballId = 0;

    ChVector<> box_pos = ChVector<>(0, 0, box_z);
    box = chrono_types::make_shared<ChBody>(chrono_types::make_shared<ChCollisionModelParallel>(), ChMaterialSurface::SMC);
    box->SetMaterialSurface(mat);
    box->SetIdentifier(ballId);
    box->SetMass(box_mass);
    box->SetInertiaXX(utils::CalcBoxGyration(hdim).diagonal());
    box->SetPos(box_pos);
    box->SetRot(ChQuaternion<>(1, 0, 0, 0));
    box->SetBodyFixed(false);
    box->SetCollide(true);
    box->GetCollisionModel()->ClearModel();
    utils::AddBoxGeometry(box.get(), hdim);
    box->GetCollisionModel()->BuildModel();
    sys->AddBody(box);
    ballId++;

    double mass = 1.0;
    ChVector<> inertia = (2.0 / 5.0) * mass * sphere_radius * sphere_radius * ChVector<>(1, 1, 1);

    for (int x = -1; x <= 1; x++) {
        for (int y = -1; y <= 1; y++) {
            if ((x * y == 0) && (x + y != 0))
                continue;
            ChVector<> pos = ChVector<>(x, y, sphere_z);
            auto ball = chrono_types::make_shared<ChBody>(chrono_types::make_shared<ChCollisionModelParallel>(), ChMaterialSurface::SMC);
            ball->SetMaterialSurface(mat);
            ball->SetIdentifier(ballId);
            ball->SetMass(mass);
            ball->SetInertiaXX(inertia);
            ball->SetPos(pos);
            ball->SetRot(ChQuaternion<>(1, 0, 0, 0));
            ball->SetBodyFixed(true);
            ball->SetCollide(true);
            ball->GetCollisionModel()->ClearModel();
            utils::AddSphereGeometry(ball.get(), sphere_radius);
            ball->GetCollisionModel()->BuildModel();
            sys->AddBody(ball);
            ballId++;
        }
    }
}

// -----------------------------------------------------------------------------
// Create the system, specify simulation parameters, and run simulation loop.
// -----------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    // Simulation parameters
    // ---------------------
    double gravity = 10;
    double time_step = 1e-4;
    double time_end = 0.5;
    double out_fps = 100;

    double push = atof(argv[1]);

    if (!filesystem::create_directory(filesystem::path(out_folder))) {
        std::cout << "Error creating directory " << out_folder << std::endl;
        return 1;
    }

    const std::string removeFiles = (std::string("rm ") + out_folder + "/* ");
    std::cout << removeFiles << std::endl;
    system(removeFiles.c_str());

    int threads = 8;

    // Create system
    // -------------
    ChSystemParallelSMC msystem;
    // Set number of threads.
    int max_threads = CHOMPfunctions::GetNumProcs();
    if (threads > max_threads)
        threads = max_threads;
    
    CHOMPfunctions::SetNumThreads(threads);

    // Set gravitational acceleration
    msystem.Set_G_acc(ChVector<>(0, 0, -gravity));
    msystem.GetSettings()->collision.narrowphase_algorithm = NarrowPhaseType::NARROWPHASE_HYBRID_MPR;
    msystem.GetSettings()->collision.collision_envelope = 1e-4;
    msystem.GetSettings()->collision.bins_per_axis = vec3(10, 10, 10);
    // Create the fixed and moving bodies
    // ----------------------------------

    CreateModel(&msystem);

    // Perform the simulation
    // ----------------------
    int num_steps = std::ceil(time_end / time_step);
    int out_steps = std::ceil((1.0 / time_step) / out_fps);
    int out_frame = 0;
    double time = 0;

    // Number of steps
    int sim_frame = 0;
    int next_out_frame = 0;

#if 0
    opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
    gl_window.Initialize(1280, 720, "ballsDVI", &msystem);
    gl_window.SetCamera(ChVector<>(0, -3, 0), ChVector<>(0, -2, 0), ChVector<>(0, 0, 1));
    gl_window.Pause();
    // Uncomment the following two lines for the OpenGL manager to automatically
    // run the simulation in an infinite loop.
    // gl_window.StartDrawLoop(time_step);
    // return 0;

    while (time < 10) {
        if (gl_window.Active()) {
            gl_window.DoStepDynamics(time_step);
            gl_window.Render();
            time += time_step;
        }
    }
#else
    // Run simulation for specified time
    for (int i = 0; i < num_steps; i++) {
        double t = msystem.GetChTime();
        double f = push * friction * box_mass * gravity;

        if (i > num_steps / 2) {
            box->Empty_forces_accumulators();
            box->Accumulate_force(ChVector<>(i / num_steps * f, 0, 0), ChVector<>(0, 0, -hdim.z()), true);
        }
        msystem.DoStepDynamics(time_step);
        // If enabled, output data for PovRay postprocessing.
        if (i == next_out_frame) {
            std::cout << time << ", Nc= " << msystem.data_manager->host_data.bids_rigid_rigid.size() << std::endl;
            printf("tangential force =%f\n", f);

            writeCSV(&msystem, out_frame, out_folder);
            next_out_frame += out_steps;
            std::ofstream ofile(out_folder + "F_SCM_" + std::to_string(out_frame) + ".txt");

            custom_vector<vec2>& pairs = msystem.data_manager->host_data.bids_rigid_rigid;
            custom_vector<real3>& gamma_N = msystem.data_manager->host_data.contact_force_N;
            custom_vector<real3>& gamma_T = msystem.data_manager->host_data.contact_force_T;
            ofile << "bi,bj,Fn,Ft\n";
            int N = pairs.size();
            std::cout << "num contacts: " << N << std::endl;
            for (int i = 0; i < N; i++)
                ofile << pairs[i].x << "," << pairs[i].y << "," << Length(gamma_N[i]) << "," << Length(gamma_T[i])
                      << "\n";
            out_frame++;

            ofile.close();
        }

        time += time_step;
    }

#endif

    return 0;
}
