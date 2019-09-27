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

#include "chrono_parallel/physics/Ch3DOFContainer.h"

#ifdef CHRONO_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

using namespace chrono;
using namespace chrono::collision;
std::shared_ptr<ChBody> box;

bool perfrom_regularization = false;
bool use_compliance = false;
bool regularize_tangential = true;
bool randomize_initials = false;
double reg_alpha0 = 10;
double reg_alpha_tau = 1e-6;
double reg_alpha_final = 9.;
double box_mass = 1.0;
double friction = 0.5f;
std::string out_folder = "CannonballNSC_reg/";

int num_ball_x = 20;
double sphere_radius = 0.1;
real tolerance = 1e-3;

// -----------------------------------------------------------------------------
// save csv data file
// -----------------------------------------------------------------------------
void writeCSV(ChSystemParallel* msystem, int out_frame) {
    char filename2[100];
    sprintf(filename2, "%s/NSC_%04d.csv", out_folder.c_str(), out_frame + 1);

    const std::string& delim = ",";
    utils::CSV_writer csv(delim);
    int numMarkers = msystem->data_manager->host_data.pos_rigid.size();
    csv << "t,x,y,z,vx,vy,vz,|U|" << std::endl;
    for (int i = 0; i < numMarkers; i++) {
        real3 pos3 = msystem->data_manager->host_data.pos_rigid[i];
        real vx = msystem->data_manager->host_data.v[6 * i];
        real vy = msystem->data_manager->host_data.v[6 * i + 1];
        real vz = msystem->data_manager->host_data.v[6 * i + 2];
        real u = sqrt(vx * vx + vy * vy + vz * vz);
        csv << msystem->GetChTime() << pos3.x << pos3.y << pos3.z << vx << vy << vz << std::endl;
    }

    csv.write_to_file(filename2);
}

// -----------------------------------------------------------------------------
// Create the falling spherical objects in a uniform rectangular grid.
// -----------------------------------------------------------------------------
void CreateModel(ChSystemParallel* sys) {
    double sphere_mass = 1.0;
    double envelope = sphere_radius * 0.05;

    // Create the middle ball material
    std::shared_ptr<ChMaterialSurfaceNSC> mat = std::make_shared<ChMaterialSurfaceNSC>();
    mat->SetFriction(friction);

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
                auto ball = std::make_shared<ChBody>(std::make_shared<ChCollisionModelParallel>());
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
    double box_x = 2 * shift * 1.5;
    ChVector<double> box_dim(box_x, box_x, 0.2);
    ChVector<double> box_pos(shift, shift, -box_dim.z() / 2);
    box = std::make_shared<ChBody>(std::make_shared<ChCollisionModelParallel>());
    box->SetMaterialSurface(mat);
    box->SetIdentifier(ballId);
    box->SetMass(box_mass);
    box->SetInertiaXX(utils::CalcBoxGyration(box_dim / 2).Get_Diag());
    box->SetPos(box_pos);
    box->SetRot(QUNIT);
    box->SetBodyFixed(true);
    box->SetCollide(true);
    box->GetCollisionModel()->ClearModel();
    utils::AddBoxGeometry(box.get(), box_dim / 2);
    box->GetCollisionModel()->BuildModel();
    sys->AddBody(box);
}

// -----------------------------------------------------------------------------
// Create the system, specify simulation parameters, and run simulation loop.
// -----------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    if (!filesystem::create_directory(filesystem::path(out_folder))) {
        std::cout << "Error creating directory " << out_folder << std::endl;
        return 1;
    }

    const std::string removeFiles = (std::string("rm ") + out_folder + "/* ");
    std::cout << removeFiles << std::endl;
    system(removeFiles.c_str());

    int threads = 8;
    int solver = 1;
    int max_iteration = 200;
    bool enable_alpha_init;
    bool enable_cache_step;

    if (argc == 2) {
        solver = atoi(argv[1]);
    }
    if (argc == 5) {
        solver = atoi(argv[1]);
        perfrom_regularization = atoi(argv[2]);
        use_compliance = atoi(argv[3]);
        regularize_tangential = atoi(argv[4]);
    }
    if (argc == 7) {
        solver = atoi(argv[1]);
        perfrom_regularization = atoi(argv[2]);
        use_compliance = atoi(argv[3]);
        regularize_tangential = atoi(argv[4]);
        reg_alpha0 = atof(argv[5]);
        randomize_initials = atoi(argv[6]);
    }

    std::cout << "solver type " << solver << ", perfrom_regularization " << perfrom_regularization << std::endl;

    // Simulation parameters
    // ---------------------

    double gravity = 10;
    double time_step = 5e-3;
    double time_end = 0.5;
    double out_fps = 100;

    // Create system
    // -------------

    ChSystemParallelNSC msystem;

    // Set number of threads.
    int max_threads = CHOMPfunctions::GetNumProcs();
    if (threads > max_threads)
        threads = max_threads;
    msystem.SetParallelThreadNumber(threads);
    CHOMPfunctions::SetNumThreads(threads);

    // Set gravitational acceleration
    msystem.Set_G_acc(ChVector<>(0, 0, -gravity));

    // Set solver parameters

    msystem.GetSettings()->solver.perfrom_regularization = perfrom_regularization;
    msystem.GetSettings()->solver.use_compliance_info = use_compliance;
    msystem.GetSettings()->solver.regularize_tangential = regularize_tangential;
    msystem.GetSettings()->solver.randomize_initial_guess = randomize_initials;
    msystem.GetSettings()->solver.reg_alpha0 = reg_alpha0;
    msystem.GetSettings()->solver.reg_alpha_tau = reg_alpha_tau;
    msystem.GetSettings()->solver.reg_alpha_final = reg_alpha_final;

    msystem.GetSettings()->solver.solver_mode = SolverMode::SLIDING;
    msystem.GetSettings()->solver.max_iteration_normal = max_iteration;
    msystem.GetSettings()->solver.max_iteration_sliding = max_iteration;
    msystem.GetSettings()->solver.max_iteration_spinning = 0;
    msystem.GetSettings()->solver.max_iteration_bilateral = 0;
    msystem.GetSettings()->solver.tolerance = tolerance;
    msystem.GetSettings()->solver.tol_speed = tolerance;

    msystem.GetSettings()->solver.alpha = 0;
    msystem.GetSettings()->solver.use_power_iteration = enable_alpha_init;
    msystem.GetSettings()->solver.cache_step_length = enable_cache_step;
    msystem.GetSettings()->solver.contact_recovery_speed = 0.0;
    if (solver == 0) {
        msystem.ChangeSolverType(SolverType::JACOBI);
    } else if (solver == 1) {
        msystem.ChangeSolverType(SolverType::GAUSS_SEIDEL);
    } else if (solver == 2) {
        msystem.ChangeSolverType(SolverType::APGD);
    } else if (solver == 3) {
        msystem.ChangeSolverType(SolverType::APGDREF);
    } else if (solver == 4) {
        msystem.ChangeSolverType(SolverType::BB);
    } else if (solver == 5) {
        msystem.ChangeSolverType(SolverType::SPGQP);
    }

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
        std::cout << time << " " << msystem.data_manager->measures.solver.residual << " "
                  << msystem.data_manager->measures.solver.total_iteration << " "
                  << msystem.data_manager->system_timer.GetTime("ChLcpSolverParallel_Solve") << std::endl;
        total_its += msystem.data_manager->measures.solver.total_iteration;

        // If enabled, output data for PovRay postprocessing.
        if (i == next_out_frame) {
            writeCSV(&msystem, out_frame);
            out_frame++;
            next_out_frame += out_steps;
            std::ofstream ofile(out_folder + "F_NSC_" + std::to_string(out_frame) + ".txt");
            DynamicVector<real>& gamma = msystem.data_manager->host_data.gamma;
            int N = gamma.size() / 3;
            for (int i = 0; i < N; i++)
                ofile << gamma[i] / time_step << "," << gamma[N + 2 * i] / time_step << ","
                      << gamma[N + 2 * i + 1] / time_step << "\n";

            ofile.close();
        }

        time += time_step;
    }

#endif

    std::cout << "==============================" << std::endl
              << "solver = " << solver << std::endl
              << "perfrom_regularization = " << perfrom_regularization << std::endl
              << "randomize_initials = " << randomize_initials << std::endl
              << "use_compliance = " << use_compliance << std::endl
              << "regularize_tangential = " << regularize_tangential << std::endl
              << "reg_alpha0 = " << reg_alpha0 << std::endl
              << "Average itrs / step = " << total_its / num_steps << std::endl
              << "==============================" << std::endl;

    return 0;
}
