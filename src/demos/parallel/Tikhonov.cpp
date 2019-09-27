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

#include "chrono_parallel/physics/Ch3DOFContainer.h"

#ifdef CHRONO_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

using namespace chrono;
using namespace chrono::collision;
std::shared_ptr<ChBody> box;

bool perfrom_regularization = true;
bool use_compliance = false;
bool regularize_tangential = true;
bool randomize_initials = false;
double reg_alpha0 = 10;
double reg_alpha_tau = 1e-10;
double reg_alpha_final = 9.;
double box_mass = 1.0;
double friction = 0.5f;
const char* out_folder = "../Tikhonov_DVI/";
ChVector<> box_dim;

void OutputBlazeMatrix(CompressedMatrix<real> src, std::string filename) {
    const char* numformat = "%.16g";
    ChStreamOutAsciiFile stream(filename.c_str());
    stream.SetNumFormat(numformat);
    for (int i = 0; i < src.rows(); ++i) {
        for (CompressedMatrix<real>::Iterator it = src.begin(i); it != src.end(i); ++it) {
            stream << it->value() << ",";
        }
        stream << "\n";
    }
}

void OutputBlazeVector(DynamicVector<real> src, std::string filename) {
    const char* numformat = "%.16g";
    ChStreamOutAsciiFile stream(filename.c_str());
    stream.SetNumFormat(numformat);

    for (int i = 0; i < src.size(); i++)
        stream << src[i] << ",";
}

// -----------------------------------------------------------------------------
// Create the falling spherical objects in a uniform rectangular grid.
// -----------------------------------------------------------------------------
void CreateModel(ChSystemParallel* sys) {
    double sphere_mass = 1.0;
    double sphere_radius = 0.1;
    double sphere_z = -sphere_radius;
    ChVector<> hdim(2, 2, 0.2);
    double envelope = sphere_radius * 0.05;
    double box_z = sphere_z + sphere_radius + hdim.z() / 2 + envelope;
    box_dim = ChVector<double>(2, 2, 0.2);

    // Create the box material
    std::shared_ptr<ChMaterialSurfaceNSC> mat_box = std::make_shared<ChMaterialSurfaceNSC>();
    mat_box->SetFriction(friction);

    // Create the middle ball material
    std::shared_ptr<ChMaterialSurfaceNSC> mat = std::make_shared<ChMaterialSurfaceNSC>();
    mat->SetFriction(friction);
    if (use_compliance) {
        mat->SetCompliance(6.0e-5);
        mat->SetComplianceT(6.0e-5);
    }

    // Create the corner balls' material
    std::shared_ptr<ChMaterialSurfaceNSC> mat2 = std::make_shared<ChMaterialSurfaceNSC>();
    mat2->SetFriction(friction);
    if (use_compliance) {
        mat2->SetCompliance(1.0e-5);
        mat2->SetComplianceT(1.0e-5);
    }

    ChVector<> box_pos = ChVector<>(0, 0, box_z);
    box = std::make_shared<ChBody>(std::make_shared<ChCollisionModelParallel>());
    box->SetMaterialSurface(mat_box);
    box->SetIdentifier(0);
    box->SetMass(box_mass);
    box->SetInertiaXX(utils::CalcBoxGyration(hdim / 2).Get_Diag());
    box->SetPos(box_pos);
    box->SetRot(QUNIT);
    box->SetBodyFixed(false);
    box->SetCollide(true);
    box->GetCollisionModel()->ClearModel();
    utils::AddBoxGeometry(box.get(), hdim / 2);
    box->GetCollisionModel()->BuildModel();
    sys->AddBody(box);

    int ballId = 1;
    double mass = 1.0;
    ChVector<> inertia = (2.0 / 5.0) * mass * sphere_radius * sphere_radius * ChVector<>(1, 1, 1);

    for (int x = -1; x <= 1; x++) {
        for (int y = -1; y <= 1; y++) {
            if ((x * y == 0) && (x + y != 0))
                continue;
            ChVector<> pos = ChVector<>(x, y, sphere_z);
            auto ball = std::make_shared<ChBody>(std::make_shared<ChCollisionModelParallel>());
            // Use a different material type for the middle ball
            if (x == 0 && y == 0)
                ball->SetMaterialSurface(mat2);
            else
                ball->SetMaterialSurface(mat);

            printf("ball %d was added.\n", ballId);
            ball->SetIdentifier(ballId++);
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
        }
    }
}

// -----------------------------------------------------------------------------
// Create the system, specify simulation parameters, and run simulation loop.
// -----------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    int threads = 8;

    int solver = 1;  // 0:BB, 1:Jacobi, 2: APGD
    int max_iteration = 10000;
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
    double time_step = 1e-3;
    double time_end = 0.2;
    double out_fps = 50;

    real tolerance = 1e-8;

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
    msystem.GetSettings()->solver.tolerance = tolerance / time_step;
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

    msystem.GetSettings()->collision.collision_envelope = 1e-3;
    msystem.GetSettings()->collision.bins_per_axis = vec3(10, 10, 10);

    // Create the fixed and moving bodies
    // ----------------------------------

    CreateModel(&msystem);

    // Perform the simulation
    // ----------------------
    int total_its = 0;
    int num_steps = std::ceil(time_end / time_step);
    int out_steps = std::ceil((1 / time_step) / out_fps);
    int out_frame = 0;
    double time = 0;
#if 0
    opengl::ChOpenGLWindow& gl_window = opengl::ChOpenGLWindow::getInstance();
    gl_window.Initialize(1280, 720, "ballsDVI", &msystem);
    gl_window.SetCamera(ChVector<>(0, -10, 0), ChVector<>(0, 0, 0), ChVector<>(0, 0, 1));
    gl_window.Pause();
    // Uncomment the following two lines for the OpenGL manager to automatically
    // run the simulation in an infinite loop.
    // gl_window.StartDrawLoop(time_step);
    // return 0;

    while (time < 20) {
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

    std::ofstream ofile("residual" + std::to_string(solver) + ".txt");
    // Run simulation for specified time

    for (int i = 0; i < num_steps; i++) {
        if (i > num_steps / 2) {
            box->Empty_forces_accumulators();
            //            box->Accumulate_force(ChVector<>(friction * box_mass * gravity / 2, 0, 0), box->GetPos(),
            //            true);
            box->Accumulate_force(ChVector<>(friction * box_mass * gravity / 2, 0, 0),
                                  ChVector<>(0, 0, -box_dim.z() / 2), true);
        }
        msystem.DoStepDynamics(time_step);
        std::cout << time << " " << msystem.data_manager->measures.solver.residual << " "
                  << msystem.data_manager->measures.solver.total_iteration << " "
                  << msystem.data_manager->system_timer.GetTime("ChLcpSolverParallel_Solve") << std::endl;
        total_its += msystem.data_manager->measures.solver.total_iteration;

        DynamicVector<real>& gamma = msystem.data_manager->host_data.gamma;
        int N = gamma.size() / 3;
        for (int i = 0; i < N; i++)
            std::cout << gamma[i] / time_step << "," << gamma[N + 2 * i] / time_step << ","
                      << gamma[N + 2 * i + 1] / time_step << std::endl;
        std::cout << "------------------" << std::endl;

        time += time_step;
    }
    ofile.close();

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
