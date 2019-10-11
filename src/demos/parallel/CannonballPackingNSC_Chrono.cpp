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

#include "chrono/ChConfig.h"
#include "chrono/utils/ChUtilsCreators.h"
#include "chrono/utils/ChUtilsInputOutput.h"
#include "chrono/physics/ChSystemNSC.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/physics/ChContactContainer.h"
#include "chrono/physics/ChContactContainerNSC.h"
#include "chrono/solver/ChSystemDescriptor.h"

#include "chrono/utils/ChUtilsGeometry.h"
#include "chrono/utils/ChUtilsGenerators.h"
#include "chrono_thirdparty/filesystem/path.h"
#include "chrono_thirdparty/filesystem/resolver.h"

#ifdef CHRONO_OPENGL
#include "chrono_opengl/ChOpenGLWindow.h"
#endif

using namespace chrono;
using namespace chrono::collision;
using namespace std;

shared_ptr<ChBody> box;

double friction = 0.f;
string out_folder = "NSC_Chrono_/";

int num_ball_x = 5;
double sphere_radius = 0.1;
double tolerance = 0.0;

// Custom contact container -- get access to the contact lists in the base
// class.
class MyContactContainer : public ChContactContainerNSC {
  public:
    vector<shared_ptr<ChBody>> body_list;
    MyContactContainer(vector<shared_ptr<ChBody>> bodyList) { body_list = bodyList; }
    MyContactContainer() {}
    // Traverse the list contactlist_6_6
    void DoStuffWithContainer(int out_frame) {
        // Look at ChContactContainerNSC and there are other lists (than
        // contactlist_6_6) that you might want to go over
        auto iter = contactlist_6_6.begin();
        int num_contact = 0;
        ofstream ofile(out_folder + "F_NSC_" + to_string(out_frame) + ".txt");
        ofile << "bi,bj,Fn,Ft,phi\n";
        while (iter != contactlist_6_6.end()) {
            ChContactable* objA = (*iter)->GetObjA();
            ChContactable* objB = (*iter)->GetObjB();
            ChVector<> p1 = (*iter)->GetContactP1();
            ChVector<> p2 = (*iter)->GetContactP2();
            ChVector<> F = (*iter)->GetContactForce();
            ChVector<> N = (*iter)->GetContactNormal();
            ChCoordsys<> CS = (*iter)->GetContactCoords();
            double CD = (*iter)->GetContactDistance();
            int Nc = body_list.size() / 3;
            double Tan = std::sqrt(F.y() * F.y() + F.z() * F.z());
            ofile << objA->GetPhysicsItem()->GetIdentifier() << "," << objB->GetPhysicsItem()->GetIdentifier() << ","
                  << F.x() << "," << Tan << "," << CD << "\n";
            num_contact++;
            ++iter;
        }
        ofile.close();
    }
};

// -----------------------------------------------------------------------------
// save csv data file
// -----------------------------------------------------------------------------
void writeCSV(ChSystemNSC& msystem, int out_frame) {
    char filename2[100];
    sprintf(filename2, "%s/Chrono_%04d.csv", out_folder.c_str(), out_frame + 1);

    const string& delim = ",";
    utils::CSV_writer csv(delim);
    vector<shared_ptr<ChBody>> bodyList = msystem.Get_bodylist();
    int numBodies = bodyList.size();
    csv << "x,y,z,vx,vy,vz,|U|" << endl;
    for (int i = 0; i < numBodies; i++) {
        ChVector<> pos = bodyList[i]->GetPos();
        ChVector<> vel = bodyList[i]->GetPos_dt();
        csv << pos.x() << pos.y() << pos.z() << vel.x() << vel.y() << vel.z() << vel.Length() << endl;
    }
    csv.write_to_file(filename2);
}

// -----------------------------------------------------------------------------
// Create the falling spherical objects in a uniform rectangular grid.
// -----------------------------------------------------------------------------
void CreateModel(ChSystemNSC& sys) {
    double sphere_mass = 1.0;
    double envelope = sphere_radius * 0.02;

    // Create the middle ball material
    shared_ptr<ChMaterialSurfaceNSC> mat = make_shared<ChMaterialSurfaceNSC>();
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
                auto msphereBody =
                    make_shared<ChBodyEasySphere>(0.1,                                              // radius size
                                                  1 / (4. / 3. * CH_C_PI * pow(sphere_radius, 3)),  // density
                                                  true,                                             // collide enable?
                                                  true);                                            // visualization?
                msphereBody->SetBodyFixed(_num_ball_ == num_ball_x);
                msphereBody->SetIdentifier(ballId);
                msphereBody->SetPos(pos);
                msphereBody->GetMaterialSurfaceNSC()->SetFriction(friction);
                sys.Add(msphereBody);
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
    // Simulation parameters
    // ---------------------

    double gravity = 10;
    double time_step = 1e-3;
    double time_end = 0.2;
    double out_fps = 100;
    int threads = 1;
    int solver = 1;
    int max_iteration = 2000;
    bool enable_alpha_init;
    bool enable_cache_step;
    if (argc == 3) {
        solver = atoi(argv[1]);
        out_folder = argv[2];
        cout << out_folder << endl;
    }
    if (!filesystem::create_directory(filesystem::path(out_folder))) {
        cout << "Error creating directory " << out_folder << endl;
        return 1;
    }

    out_folder = out_folder + "/";
    const string removeFiles = (string("rm ") + out_folder + "/* ");
    cout << removeFiles << endl;
    system(removeFiles.c_str());

    // Create system
    // -------------

    ChSystemNSC msystem;

    // Set number of threads.
    int max_threads = CHOMPfunctions::GetNumProcs();
    if (threads > max_threads)
        threads = max_threads;
    msystem.SetParallelThreadNumber(threads);
    CHOMPfunctions::SetNumThreads(threads);

    // Set gravitational acceleration
    msystem.Set_G_acc(ChVector<>(0, 0, -gravity));

    if (solver == 0) {
        msystem.SetSolverType(ChSolver::Type::JACOBI);
    } else if (solver == 1) {
        msystem.SetSolverType(ChSolver::Type::PCG);
    } else if (solver == 2) {
        msystem.SetSolverType(ChSolver::Type::APGD);
    } else if (solver == 3) {
        msystem.SetSolverType(ChSolver::Type::APGD);
    } else if (solver == 4) {
        msystem.SetSolverType(ChSolver::Type::BARZILAIBORWEIN);
    }

    // Create the fixed and moving bodies
    // ----------------------------------

    CreateModel(msystem);

    // Perform the simulation
    // ----------------------
    int total_its = 0;
    int num_steps = ceil(time_end / time_step);
    int out_steps = ceil((1.0 / time_step) / out_fps);
    int out_frame = 0;
    double time = 0;

    // Number of steps
    int sim_frame = 0;
    int next_out_frame = 0;

    // mphysicalSystem.SetUseSleeping(true);
    msystem.SetMaxItersSolverSpeed(max_iteration);
    msystem.SetSolverWarmStarting(false);
    msystem.SetMaxPenetrationRecoverySpeed(0.0);

    vector<shared_ptr<ChBody>> blist = msystem.Get_bodylist();
    auto container = make_shared<MyContactContainer>(blist);
    msystem.SetContactContainer(container);
    // Run simulation for specified time
    for (int i = 0; i < num_steps; i++) {
        msystem.DoStepDynamics(time_step);
        std::cout << "T= " << time_step * i << std::endl;
        if (i == next_out_frame) {
            shared_ptr<ChSystemDescriptor> idesc = msystem.GetSystemDescriptor();
            shared_ptr<ChContactContainerNSC> iconcon =
                dynamic_pointer_cast<ChContactContainerNSC>(msystem.GetContactContainer());
            try {
                container->DoStuffWithContainer(out_frame);
            } catch (ChException myexc) {
                cout << myexc.what();
            }

            writeCSV(msystem, out_frame);
            out_frame++;
            next_out_frame += out_steps;
        }

        time += time_step;
    }

    cout << "==============================" << endl
         << "solver = " << solver << endl
         << "==============================" << endl;

    return 0;
}
