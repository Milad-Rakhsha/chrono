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

#include "chrono_parallel/physics/Ch3DOFContainer.h"

#ifdef CHRONO_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

using namespace chrono;
using namespace chrono::collision;
std::shared_ptr<ChBody> box;

double box_mass = 1.0;
double friction = 0.5f;
std::string out_folder = "CannonballSMC";
int num_ball_x = 5;

double sphere_radius = 0.1;
// -----------------------------------------------------------------------------
// Create the falling spherical objects in a uniform rectangular grid.
// -----------------------------------------------------------------------------
void CreateModel(ChSystemParallel* sys) {
    double sphere_mass = 1.0;
    double envelope = sphere_radius * 0.02;

    // Create the middle ball material
    // Material properties (same on bin and balls)
    float Y = 1e10f;
    float cr = 0.0f;
    auto mat = chrono_types::make_shared<ChMaterialSurfaceSMC>();
    mat->SetFriction(friction);
    mat->SetYoungModulus(Y);
    mat->SetRestitution(cr);
    mat->SetAdhesion(0);  // Magnitude of the adhesion in Constant adhesion model

    int ballId = 0;
    double mass = 1.0;
    ChVector<> inertia = (2.0 / 5.0) * mass * sphere_radius * sphere_radius * ChVector<>(1, 1, 1);

    int _num_ball_ = num_ball_x;
    double sphere_z = sphere_radius;
    double shift = 0;
    while (_num_ball_ >= 0) {
        for (int x = 0; x <= _num_ball_; x++) {
            for (int y = 0; y <= _num_ball_; y++) {
                ChVector<> pos = ChVector<>(x * sphere_radius * 2 + shift, y * sphere_radius * 2 + shift, sphere_z);
                auto ball =
                    chrono_types::make_shared<ChBody>(chrono_types::make_shared<ChCollisionModelParallel>(), ChMaterialSurface::SMC);
                ball->SetMaterialSurface(mat);
                ball->SetIdentifier(ballId);
                ball->SetMass(mass);
                ball->SetInertiaXX(inertia);
                ball->SetPos(pos);
                ball->SetRot(ChQuaternion<>(1, 0, 0, 0));
                ball->SetBodyFixed(_num_ball_ == num_ball_x);
                ball->SetCollide(true);
                ball->GetCollisionModel()->ClearModel();
                utils::AddSphereGeometry(ball.get(), sphere_radius);
                ball->GetCollisionModel()->BuildModel();
                sys->AddBody(ball);
                ballId++;
            }
        }
        printf("%d ball were added.\n", ballId);
        _num_ball_--;
        shift += sphere_radius;
        sphere_z += sqrt(2.0) * sphere_radius + envelope;
    }
}

// -----------------------------------------------------------------------------
// Create the system, specify simulation parameters, and run simulation loop.
// -----------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    out_folder = out_folder + argv[1] + "/";
    // Simulation parameters
    // ---------------------
    double gravity = 10;
    double time_step = 2e-5;
    double time_end = 0.20;
    double out_fps = 100;

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
    int total_its = 0;
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
    double tmp = num_ball_x * sphere_radius;
    gl_window.SetCamera(ChVector<>(tmp, -tmp, tmp), ChVector<>(tmp, tmp, 0), ChVector<>(0, 0, 1));
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
            //            msystem.CalculateContactForces();
            //            real3 frc = msystem.GetBodyContactForce(0);
            //            std::cout << frc.x << "  " << frc.y << "  " << frc.z << std::endl;
        }
    }
#else

    // Run simulation for specified time
    for (int i = 0; i < num_steps; i++) {
        msystem.DoStepDynamics(time_step);

        // If enabled, output data for PovRay postprocessing.
        if (i == next_out_frame) {
            std::cout << time << ", Nc= " << msystem.data_manager->host_data.bids_rigid_rigid.size() << std::endl;

            writeCSV(&msystem, out_frame, out_folder);
            out_frame++;
            next_out_frame += out_steps;
            std::ofstream ofile(out_folder + "F_SCM_" + std::to_string(out_frame) + ".txt");

            custom_vector<vec2>& pairs = msystem.data_manager->host_data.bids_rigid_rigid;
            custom_vector<real3>& gamma_N = msystem.data_manager->host_data.contact_force_N;
            custom_vector<real3>& gamma_T = msystem.data_manager->host_data.contact_force_T;
            custom_vector<real3>& body_force = msystem.data_manager->host_data.ct_body_force;

            ofile << "bi,bj,Fn,Ft\n";
            int N = pairs.size();
            std::cout << "num contacts: " << N << std::endl;
            for (int i = 0; i < N; i++)
                ofile << pairs[i].x << "," << pairs[i].y << "," << Length(gamma_N[i]) << "," << Length(gamma_T[i])
                      << "\n";
            ofile.close();

            std::ofstream obfile(out_folder + "F_SCM_body_" + std::to_string(out_frame) + ".txt");
            int Nb = body_force.size();
            obfile << "i,fx,fy,fz\n";
            for (int i = 0; i < Nb; i++)
                obfile << i << "," << body_force[i].x << "," << body_force[i].y << "," << body_force[i].z << "\n";

            obfile.close();
        }

        time += time_step;
    }

#endif

    return 0;
}
