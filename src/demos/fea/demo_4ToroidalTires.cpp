// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All right reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Antonio Recuero
// =============================================================================
//
// This demo creates a toroidal geometry that may be used for simple
// tire/ground interaction (user-customizable). The tire is made up of
// ANCF shell elements and may be fully parameterized.
// Both radii of the toroidal geometry and the number of
// ANCF shell elements can be parameterized at the beginning of this file. Boolean
// variables are defined to select the addition of bodies and constraints to the system.
// Position and direction constraints attach the tire tread to the rim. Initial velocity
// and/or forces are used to move the tire forward, which remains in a perpendicular
// plane through a Plane-Plane constraint (ChLinkLockPlanePlane)
// This demo shows two ways of adding pressure to shell elements using ChLoaderPressure
// and ChLoaderUVdistributed. The latter allows for complex, user defined expressions, whereas
// the former imposes constant pressure using normal vectors defined in the element.
// =============================================================================
#include "chrono/physics/ChSystem.h"
#include "chrono/physics/ChBodyEasy.h"
#include "chrono/lcp/ChLcpIterativeMINRES.h"
#include "chrono_fea/ChElementShellANCF.h"
#include "chrono_fea/ChMesh.h"
#include "chrono_fea/ChLinkPointFrame.h"
#include "chrono_fea/ChLinkDirFrame.h"
#include "chrono/utils/ChUtilsInputOutput.h"
#include "chrono/utils/ChUtilsValidation.h"
#include "chrono/core/ChMathematics.h"
#include "chrono_mkl/ChLcpMklSolver.h"
#include "physics/ChLoadContainer.h"
#include "chrono_fea/ChVisualizationFEAmesh.h"
#include "chrono/core/ChRealtimeStep.h"
#include "chrono_irrlicht/ChBodySceneNode.h"
#include "chrono_irrlicht/ChBodySceneNodeTools.h"
#include "chrono_irrlicht/ChIrrAppInterface.h"
#include "chrono_irrlicht/ChIrrApp.h"
#include <irrlicht.h>
// Remember to use the namespace 'chrono' because all classes
// of Chrono::Engine belong to this namespace and its children...

using namespace chrono;
using namespace fea;
using namespace irr;
using namespace scene;

bool addConstRim = true;
bool addBodies = true;
bool addGroundForces = true;
bool showVisual = true;
bool addSingleLoad = false;
bool addPressureAlessandro = true;

ChSharedPtr<ChBody> BGround;
ChSharedPtr<ChBody> Body_2;       // Hub 1
ChSharedPtr<ChBody> Body_3;       // Hub 2
ChSharedPtr<ChBody> Body_4;       // Hub 3
ChSharedPtr<ChBody> Body_5;       // Hub 4
ChSharedPtr<ChBody> SimpChassis;  // Chassis body

ChSharedPtr<ChLinkPointFrame> NodePosRim;
ChSharedPtr<ChLinkDirFrame> NodeDirRim;
ChSharedPtr<ChNodeFEAxyzD> ConstrainedNode;
ChSharedPtr<ChLinkPointFrame> NodePosRim2;
ChSharedPtr<ChLinkDirFrame> NodeDirRim2;
ChSharedPtr<ChNodeFEAxyzD> ConstrainedNode2;

// Some model parameters

const int num_steps = 550;  // Number of time steps for unit test (range 1 to 4000) 550
double time_step = 0.001;   // Time step: 0.001 It works like that, but rotates at the beginning

const double ForVel = 0.0;      // Initial velocity of the tire. Applied to rim and nodes
double plate_lenght_z = 0.014;  // Thickness of the tire

// Mass of the rigid body (dumb brick)
const double RBMass = 2000;
// Location of the wheels w.r.t. rigid body's center of mass
const double Lwz = 1.1;  // Absolute value of longitudinal distance
const double Lwx = 0.7;  // Absolute value of lateral distance

// Specification of the mesh
const int numEl_Diameter = 50;  // Number of elements along diameter
const int numEl_Thread = 9;     // Number of elements along thread
const int numEl_z = 1;          // Single element along the thickness
const int N_Diameter = numEl_Diameter;
const int N_Thread = numEl_Thread + 1;
const int N_z = numEl_z + 1;
const double TorusSmallRadius = 0.195;
const double TorusRadius = 0.35;
const double Clearance = 0.00;  // Initial space between tire and ground
const double GroundLoc =
    -(TorusRadius + TorusSmallRadius + Clearance);          // Vertical position of the ground (for contact)
int TotalNumElements = numEl_Diameter * numEl_Thread;       // Number of elements in the tire
int TotalNumNodes = (numEl_Diameter) * (numEl_Thread + 1);  // Total number of nodes in the tire

double dz = plate_lenght_z / numEl_z;  // We are going to calculate distances through coordinates
double dx = 2 * CH_C_PI * (TorusRadius + TorusSmallRadius) / 2 / N_Diameter;  // Rough estimate of shell dimensions
double dy = CH_C_PI * TorusSmallRadius / numEl_Thread;

// ChLoadCustomMultiple to include basic node-Ground contact interaction
class MyLoadCustomMultiple : public ChLoadCustomMultiple {
  public:
    MyLoadCustomMultiple(std::vector<ChSharedPtr<ChLoadable>>& mloadables) : ChLoadCustomMultiple(mloadables){};

    // Compute Q=Q(x,v)
    // This is the function that you have to implement. It should return the
    // generalized Q load
    // (i.e.the force in generalized lagrangian coordinates).
    // Since here we have multiple connected ChLoadable objects (the two nodes),
    // the rule is that
    // all the vectors (load_Q, state_x, state_w) are split in the same order that
    // the loadable objects
    // are added to MyLoadCustomMultiple; in this case for instance
    // Q={Efx,Efy,Efz,Ffx,Ffy,Ffz}.
    // As this is a stiff force field, dependency from state_x and state_y must be
    // considered.
    virtual void ComputeQ(ChState* state_x,      ///< state position to evaluate Q
                          ChStateDelta* state_w  ///< state speed to evaluate Q
                          ) {
        std::vector<ChSharedPtr<ChLoadable>> NodeList;
        ChVector<> Node1_Pos;
        ChVector<> Node1_Vel;
        ChVector<> Node1_Grad;
        ChVector<> Node1_GradVel;
        this->load_Q.FillElem(0);
        const double KGround = 3e5;
        const double CGround = 8e2;
        double NormalForceNode = 0;
        double FrictionCoeff = 2.0;  // F1-like
        if (state_x && state_w) {
            for (int iii = 0; iii < TotalNumNodes; iii++) {
                Node1_Pos = state_x->ClipVector(iii * 6, 0);
                Node1_Grad = state_x->ClipVector(iii * 6 + 3, 0);
                Node1_Vel = state_w->ClipVector(iii * 6, 0);
                Node1_GradVel = state_w->ClipVector(iii * 6 + 3, 0);

                if (Node1_Pos.y < GroundLoc) {
                    NormalForceNode = KGround * abs(Node1_Pos.y - GroundLoc) - CGround * (Node1_Vel.y);
                    this->load_Q(iii * 6 + 1) = NormalForceNode;  // Fy (Vertical)
                    const double VelLimit = 0.0001;
                    // Calculation of friction forces
                    if (abs(Node1_Vel.x) > VelLimit) {
                        this->load_Q(iii * 6 + 0) =
                            -NormalForceNode * FrictionCoeff *
                            (Node1_Vel.x / sqrt((pow(Node1_Vel.x, 2) + pow(Node1_Vel.z, 2))));  // Fx (Plane x)
                    } else {
                        this->load_Q(iii * 6 + 0) =
                            -NormalForceNode * FrictionCoeff * sin(abs(Node1_Vel.x) * CH_C_PI_2 / VelLimit) *
                            (Node1_Vel.x / sqrt((pow(Node1_Vel.x, 2) + pow(Node1_Vel.z, 2))));  // Fx (Plane x)
                    }
                    if (abs(Node1_Vel.z) > VelLimit) {
                        this->load_Q(iii * 6 + 2) =
                            -NormalForceNode * FrictionCoeff *
                            (Node1_Vel.z / sqrt((pow(Node1_Vel.x, 2) + pow(Node1_Vel.z, 2))));  // Fz (Plane y)
                    } else {
                        this->load_Q(iii * 6 + 2) =
                            -NormalForceNode * FrictionCoeff * sin(abs(Node1_Vel.z) * CH_C_PI_2 / VelLimit) *
                            (Node1_Vel.z / sqrt((pow(Node1_Vel.x, 2) + pow(Node1_Vel.z, 2))));  // Fz (Plane y)
                    }
                }
                Node1_Pos.SetNull();
            }
        } else {
            // explicit integrators might call ComputeQ(0,0), null pointers mean
            // that we assume current state, without passing state_x for efficiency
        }
    }

    virtual void ComputeJacobian(ChState* state_x,       ///< state position to evaluate jacobians
                                 ChStateDelta* state_w,  ///< state speed to evaluate jacobians
                                 ChMatrix<>& mK,         ///< result dQ/dx
                                 ChMatrix<>& mR,         ///< result dQ/dv
                                 ChMatrix<>& mM)         ///< result dQ/da
    {
        double Delta = 1e-8;
        int mrows_w = this->LoadGet_ndof_w();
        // compute Q at current speed & position, x_0, v_0
        ChVectorDynamic<> Q0(mrows_w);
        this->ComputeQ(state_x, state_w);  // Q0 = Q(x, v)
        Q0 = this->load_Q;

        ChVectorDynamic<> Q1(mrows_w);
        ChVectorDynamic<> Jcolumn(mrows_w);

        // Compute K=-dQ(x,v)/dx by backward differentiation
        for (int i = 0; i < mrows_w; ++i) {
            (*state_x)(i) += Delta;            //***TODO*** use NodeIntStateIncrement
            this->ComputeQ(state_x, state_w);  // Q1 = Q(x+Dx, v)
            Q1 = this->load_Q;
            (*state_x)(i) -= Delta;                //***TODO*** use NodeIntStateIncrement
            Jcolumn = (Q1 - Q0) * (-1.0 / Delta);  // - sign because K=-dQ/dx
            this->jacobians->K.PasteMatrix(&Jcolumn, 0, i);
        }
        // Compute R=-dQ(x,v)/dv by backward differentiation
        for (int i = 0; i < mrows_w; ++i) {
            (*state_w)(i) += Delta;
            this->ComputeQ(state_x, state_w);  // Q1 = Q(x, v+Dv)
            Q1 = this->load_Q;
            (*state_w)(i) -= Delta;

            Jcolumn = (Q1 - Q0) * (-1.0 / Delta);  // - sign because R=-dQ/dv
            this->jacobians->R.PasteMatrix(&Jcolumn, 0, i);
        }
    };

    // OPTIONAL: if you want to provide an analytical jacobian, just implement the
    // following:
    //   virtual void ComputeJacobian(...)

    // Remember to set this as stiff, to enable the jacobians
    virtual bool IsStiff() { return false; }
};

int main(int argc, char* argv[]) {
    // Definition of the model
    ChSystem my_system;
    std::vector<ChSharedPtr<ChMesh>> TireMeshes(4);

    ChIrrApp application(&my_system, L"ANCF Rolling Tire", core::dimension2d<u32>(1080, 800), false);
    // Easy shortcuts to add camera, lights, logo and sky in Irrlicht scene:
    application.AddTypicalLogo();
    application.AddTypicalSky();
    application.AddTypicalLights();
    application.AddTypicalCamera(core::vector3df(0.5f, 0.5f, 1.15f),   // camera location
                                 core::vector3df(0.65f, 0.0f, 0.0f));  // "look at" location
    application.AddTypicalLights(irr::core::vector3df(30.f, -30.f, 100.f), irr::core::vector3df(30.f, 50.f, 100.f), 160,
                                 70);
    utils::CSV_writer out("\t");
    out.stream().setf(std::ios::scientific | std::ios::showpos);
    out.stream().precision(7);
    // Main loop for the definition of 4 meshes

    for (int TireNo = 0; TireNo < 4; TireNo++) {
        //  Fixing constraints, initial coordinates and velocities
        double thethaTorus = 0.0;
        double phiTorus = 0.0;
        double loc_x = 0.0;
        double loc_y = 0.0;
        double loc_z = 0.0;

        for (int j = 0; j < N_Diameter; j++) {
            for (int i = 0; i < N_Thread; i++) {
                thethaTorus = -CH_C_PI / 2 + CH_C_PI * i / (N_Thread - 1);
                thethaTorus = -CH_C_PI / 2 + CH_C_PI * i / (N_Thread - 1);
                phiTorus = 2 * CH_C_PI * j / N_Diameter;
                switch (TireNo) {
                    case 0:
                        loc_x = TorusSmallRadius * sin(thethaTorus) - Lwx;
                        loc_y = (TorusRadius + TorusSmallRadius * cos(thethaTorus)) * sin(phiTorus);
                        loc_z = (TorusRadius + TorusSmallRadius * cos(thethaTorus)) * cos(phiTorus) + Lwz;
                        break;
                    case 1:
                        loc_x = TorusSmallRadius * sin(thethaTorus) + Lwx;
                        loc_y = (TorusRadius + TorusSmallRadius * cos(thethaTorus)) * sin(phiTorus);
                        loc_z = (TorusRadius + TorusSmallRadius * cos(thethaTorus)) * cos(phiTorus) + Lwz;
                        break;
                    case 2:
                        loc_x = TorusSmallRadius * sin(thethaTorus) + Lwx;
                        loc_y = (TorusRadius + TorusSmallRadius * cos(thethaTorus)) * sin(phiTorus);
                        loc_z = (TorusRadius + TorusSmallRadius * cos(thethaTorus)) * cos(phiTorus) - Lwz;
                        break;
                    case 3:
                        loc_x = TorusSmallRadius * sin(thethaTorus) - Lwx;
                        loc_y = (TorusRadius + TorusSmallRadius * cos(thethaTorus)) * sin(phiTorus);
                        loc_z = (TorusRadius + TorusSmallRadius * cos(thethaTorus)) * cos(phiTorus) - Lwz;
                        break;
                }

                double dir_x = sin(thethaTorus);
                double dir_y = cos(thethaTorus) * sin(phiTorus);
                double dir_z = cos(thethaTorus) * cos(phiTorus);

                // Write position and gradient of nodes to a file and plot them in Matlab
                out << i + j * N_Thread << loc_x << loc_y << loc_z << dir_x << dir_y << dir_z << std::endl;
                out.write_to_file("../TorusTire.txt");

                double vel_x = 0;
                double vel_y = 0;
                double vel_z = ForVel;

                // Create the node
                ChSharedPtr<ChNodeFEAxyzD> node(
                    new ChNodeFEAxyzD(ChVector<>(loc_x, loc_y, loc_z), ChVector<>(dir_x, dir_y, dir_z)));
                // Apply initial velocities
                node.DynamicCastTo<ChNodeFEAxyz>()->SetPos_dt(ChVector<>(vel_x, vel_y, vel_z));
                node->SetMass(0);
                // Add node to mesh
                TireMeshes[TireNo]->AddNode(node);
            }
        }  // Here we finish loop over each tire

        // Create an isotropic material.
        // All elements share the same material.
        ChSharedPtr<ChMaterialShellANCF> mat(new ChMaterialShellANCF(500, 9.0e7, 0.3));

        // Create the elements
        int node0, node1, node2, node3;         // Node numbering
        for (int j = 0; j < N_Diameter; j++) {  // Start node numbering by zero
            for (int i = 0; i < numEl_Thread; i++) {
                if (j == N_Diameter - 1) {
                    node0 = i + j * (numEl_Thread + 1);
                    node1 = i;
                    node2 = i + 1;
                    node3 = i + 1 + j * (numEl_Thread + 1);
                } else {
                    node0 = i + j * (numEl_Thread + 1);
                    node1 = i + (j + 1) * (numEl_Thread + 1);
                    node2 = i + 1 + (j + 1) * (numEl_Thread + 1);
                    node3 = i + 1 + j * (numEl_Thread + 1);
                }
                // Create the element and set its nodes.
                ChSharedPtr<ChElementShellANCF> element(new ChElementShellANCF);
                element->SetNodes(TireMeshes[TireNo]->GetNode(node0).DynamicCastTo<ChNodeFEAxyzD>(),
                                  TireMeshes[TireNo]->GetNode(node1).DynamicCastTo<ChNodeFEAxyzD>(),
                                  TireMeshes[TireNo]->GetNode(node3).DynamicCastTo<ChNodeFEAxyzD>(),
                                  TireMeshes[TireNo]->GetNode(node2).DynamicCastTo<ChNodeFEAxyzD>());
                // Set element dimensions
                element->SetDimensions(dx, dy);

                // Add a single layers with a fiber angle of 0 degrees.
                element->AddLayer(dz, 0 * CH_C_DEG_TO_RAD, mat);

                // Set other element properties
                element->SetAlphaDamp(0.15);  // Structural damping for this element
                element->SetGravityOn(true);  // no gravitational forces

                // Add element to mesh
                TireMeshes[TireNo]->AddElement(element);
            }
        }  // Here the creation of elements ends
        // Switch off mesh class gravity
        TireMeshes[TireNo]->SetAutomaticGravity(false);

        // Add the mesh to the system
        my_system.Add(TireMeshes[TireNo]);

        if (addSingleLoad) {
            ChSharedPtr<ChNodeFEAxyzD> nodetip(
                TireMeshes[TireNo]->GetNode(TotalNumNodes - 1).DynamicCastTo<ChNodeFEAxyzD>());
            nodetip->SetForce(ChVector<>(0, 0, -10));
        }

    }  // Here we finish loop over all tires

    // Mark completion of system construction
    my_system.SetupInitial();

    if (addBodies) {
        // Body 1: Ground
        BGround = ChSharedPtr<ChBody>(new ChBody);
        my_system.AddBody(BGround);
        BGround->SetIdentifier(1);
        BGround->SetBodyFixed(true);
        BGround->SetCollide(false);
        BGround->SetMass(1);
        BGround->SetInertiaXX(ChVector<>(1, 1, 0.2));
        BGround->SetPos(ChVector<>(-2, 0, 0));  // Y = -1m
        ChQuaternion<> rot = Q_from_AngX(0.0);
        BGround->SetRot(rot);
        // Defining the Body 2: Rim for TireMeshes[0]
        Body_2 = ChSharedPtr<ChBody>(new ChBody);
        my_system.AddBody(Body_2);
        Body_2->SetIdentifier(2);
        Body_2->SetBodyFixed(false);
        Body_2->SetCollide(false);
        Body_2->SetMass(8);
        Body_2->SetInertiaXX(ChVector<>(0.6, 0.3, 0.3));
        Body_2->SetPos(ChVector<>(-Lwx, 0, Lwz));  // Y = -1m
        Body_2->SetRot(rot);
        // Defining the Body 3: Rim for TireMeshes[1]
        Body_3 = ChSharedPtr<ChBody>(new ChBody);
        my_system.AddBody(Body_3);
        Body_3->SetIdentifier(3);
        Body_3->SetBodyFixed(false);
        Body_3->SetCollide(false);
        Body_3->SetMass(8);
        Body_3->SetInertiaXX(ChVector<>(0.6, 0.3, 0.3));
        Body_3->SetPos(ChVector<>(Lwx, 0, Lwz));  // Y = -1m
        Body_3->SetRot(rot);
        // Defining the Body 4: Rim for TireMeshes[2]
        Body_4 = ChSharedPtr<ChBody>(new ChBody);
        my_system.AddBody(Body_4);
        Body_4->SetIdentifier(4);
        Body_4->SetBodyFixed(false);
        Body_4->SetCollide(false);
        Body_4->SetMass(8);
        Body_4->SetInertiaXX(ChVector<>(0.6, 0.3, 0.3));
        Body_4->SetPos(ChVector<>(Lwx, 0, -Lwz));  // Y = -1m
        Body_4->SetRot(rot);
        // Defining the Body 5: Rim for TireMeshes[3]
        Body_5 = ChSharedPtr<ChBody>(new ChBody);
        my_system.AddBody(Body_5);
        Body_5->SetIdentifier(4);
        Body_5->SetBodyFixed(false);
        Body_5->SetCollide(false);
        Body_5->SetMass(8);
        Body_5->SetInertiaXX(ChVector<>(0.6, 0.3, 0.3));
        Body_5->SetPos(ChVector<>(-Lwx, 0, -Lwz));  // Y = -1m
        Body_5->SetRot(rot);
        // Defining the SimpChassis
        SimpChassis = ChSharedPtr<ChBody>(new ChBody);
        my_system.AddBody(SimpChassis);
        SimpChassis->SetIdentifier(5);
        SimpChassis->SetBodyFixed(false);
        SimpChassis->SetCollide(false);
        SimpChassis->SetMass(RBMass);
        SimpChassis->SetInertiaXX(ChVector<>(900, 900, 1500));
        SimpChassis->SetPos(ChVector<>(0, 0, 0));  // Y = -1m
        SimpChassis->SetRot(rot);

        // Revolute joint between hub and rim
        ChSharedPtr<ChLinkLockRevolute> Rev1(new ChLinkLockRevolute);
        my_system.AddLink(Rev1);
        Rev1->Initialize(Body_2, SimpChassis, ChCoordsys<>(ChVector<>(-Lwx, 0, Lwz), Q_from_AngY(-CH_C_PI_2)));

        ChSharedPtr<ChLinkLockRevolute> Rev2(new ChLinkLockRevolute);
        my_system.AddLink(Rev2);
        Rev2->Initialize(Body_3, SimpChassis, ChCoordsys<>(ChVector<>(Lwx, 0, Lwz), Q_from_AngY(-CH_C_PI_2)));

        ChSharedPtr<ChLinkLockRevolute> Rev3(new ChLinkLockRevolute);
        my_system.AddLink(Rev3);
        Rev3->Initialize(Body_4, SimpChassis, ChCoordsys<>(ChVector<>(Lwx, 0, -Lwz), Q_from_AngY(-CH_C_PI_2)));

        ChSharedPtr<ChLinkLockRevolute> Rev4(new ChLinkLockRevolute);
        my_system.AddLink(Rev4);
        Rev4->Initialize(Body_5, SimpChassis, ChCoordsys<>(ChVector<>(-Lwx, 0, -Lwz), Q_from_AngY(-CH_C_PI_2)));

        my_system.Set_G_acc(ChVector<>(0, -9.81, 0));  // Hey! 4G!
    }

    // First: loads must be added to "load containers",
    // and load containers must be added to your ChSystem
    std::vector<ChSharedPtr<ChLoadContainer>> Mloadcontainer4(4);
    // Add constant pressure using ChLoaderPressure (preferred for simple, constant pressure)
    if (addPressureAlessandro) {
        for (int TireNo = 0; TireNo < 4; TireNo++) {
            for (int NoElmPre = 0; NoElmPre < TotalNumElements; NoElmPre++) {
                ChSharedPtr<ChLoad<ChLoaderPressure>> faceload(new ChLoad<ChLoaderPressure>(
                    TireMeshes[TireNo]->GetElement(NoElmPre).StaticCastTo<ChElementShellANCF>()));
                faceload->loader.SetPressure(320e3);
                faceload->loader.SetStiff(false);
                faceload->loader.SetIntegrationPoints(2);
                Mloadcontainer4[TireNo]->Add(faceload);
            }
        }
    }  // End of defining load containers for all pressure loads

    // Add constraints to the rim: ChLinkPointFrame for position constraints,
    // ChLinkDirFrame for direction constraints
    if (addConstRim) {
        for (int TireNo = 0; TireNo < 4; TireNo++) {
            int NoConstNodes = 2 * numEl_Diameter;
            int ConstNodeInx = -1;
            for (int j = 0; j < N_Diameter; j++) {
                for (int i = 0; i < N_Thread; i++) {
                    ConstNodeInx = i + j * N_Thread;  // General node index

                    if ((i + j * N_Thread) % (N_Thread) == 0) {
                        // First node to be constrained (one side of the rim)

                        switch (TireNo) {
                            case 0:
                                ConstrainedNode = ChSharedPtr<ChNodeFEAxyzD>(
                                    TireMeshes[TireNo]->GetNode(ConstNodeInx).DynamicCastTo<ChNodeFEAxyzD>());
                                NodePosRim = ChSharedPtr<ChLinkPointFrame>(new ChLinkPointFrame);
                                NodePosRim->Initialize(ConstrainedNode, Body_2);
                                my_system.Add(NodePosRim);

                                NodeDirRim = ChSharedPtr<ChLinkDirFrame>(new ChLinkDirFrame);
                                NodeDirRim->Initialize(ConstrainedNode, Body_2);
                                NodeDirRim->SetDirectionInAbsoluteCoords(ConstrainedNode->D);
                                my_system.Add(NodeDirRim);

                                // Second node to be constrained (other side of the rim)
                                ConstrainedNode2 = ChSharedPtr<ChNodeFEAxyzD>(TireMeshes[TireNo]
                                                                                  ->GetNode(ConstNodeInx + N_Thread - 1)
                                                                                  .DynamicCastTo<ChNodeFEAxyzD>());

                                NodePosRim2 = ChSharedPtr<ChLinkPointFrame>(new ChLinkPointFrame);
                                NodePosRim2->Initialize(ConstrainedNode2, Body_2);
                                my_system.Add(NodePosRim2);

                                NodeDirRim2 = ChSharedPtr<ChLinkDirFrame>(new ChLinkDirFrame);
                                NodeDirRim2->Initialize(ConstrainedNode2, Body_2);
                                NodeDirRim2->SetDirectionInAbsoluteCoords(ConstrainedNode2->D);
                                my_system.Add(NodeDirRim2);
                                break;
                            case 1:
                                ConstrainedNode = ChSharedPtr<ChNodeFEAxyzD>(
                                    TireMeshes[TireNo]->GetNode(ConstNodeInx).DynamicCastTo<ChNodeFEAxyzD>());
                                NodePosRim = ChSharedPtr<ChLinkPointFrame>(new ChLinkPointFrame);
                                NodePosRim->Initialize(ConstrainedNode, Body_3);
                                my_system.Add(NodePosRim);

                                NodeDirRim = ChSharedPtr<ChLinkDirFrame>(new ChLinkDirFrame);
                                NodeDirRim->Initialize(ConstrainedNode, Body_3);
                                NodeDirRim->SetDirectionInAbsoluteCoords(ConstrainedNode->D);
                                my_system.Add(NodeDirRim);

                                // Second node to be constrained (other side of the rim)
                                ConstrainedNode2 = ChSharedPtr<ChNodeFEAxyzD>(TireMeshes[TireNo]
                                    ->GetNode(ConstNodeInx + N_Thread - 1)
                                    .DynamicCastTo<ChNodeFEAxyzD>());

                                NodePosRim2 = ChSharedPtr<ChLinkPointFrame>(new ChLinkPointFrame);
                                NodePosRim2->Initialize(ConstrainedNode2, Body_3);
                                my_system.Add(NodePosRim2);

                                NodeDirRim2 = ChSharedPtr<ChLinkDirFrame>(new ChLinkDirFrame);
                                NodeDirRim2->Initialize(ConstrainedNode2, Body_3);
                                NodeDirRim2->SetDirectionInAbsoluteCoords(ConstrainedNode2->D);
                                my_system.Add(NodeDirRim2);
                                break;
                        case 2:
                            ConstrainedNode = ChSharedPtr<ChNodeFEAxyzD>(
                                TireMeshes[TireNo]->GetNode(ConstNodeInx).DynamicCastTo<ChNodeFEAxyzD>());
                            NodePosRim = ChSharedPtr<ChLinkPointFrame>(new ChLinkPointFrame);
                            NodePosRim->Initialize(ConstrainedNode, Body_4);
                            my_system.Add(NodePosRim);

                            NodeDirRim = ChSharedPtr<ChLinkDirFrame>(new ChLinkDirFrame);
                            NodeDirRim->Initialize(ConstrainedNode, Body_4);
                            NodeDirRim->SetDirectionInAbsoluteCoords(ConstrainedNode->D);
                            my_system.Add(NodeDirRim);

                            // Second node to be constrained (other side of the rim)
                            ConstrainedNode2 = ChSharedPtr<ChNodeFEAxyzD>(TireMeshes[TireNo]
                                ->GetNode(ConstNodeInx + N_Thread - 1)
                                .DynamicCastTo<ChNodeFEAxyzD>());

                            NodePosRim2 = ChSharedPtr<ChLinkPointFrame>(new ChLinkPointFrame);
                            NodePosRim2->Initialize(ConstrainedNode2, Body_4);
                            my_system.Add(NodePosRim2);

                            NodeDirRim2 = ChSharedPtr<ChLinkDirFrame>(new ChLinkDirFrame);
                            NodeDirRim2->Initialize(ConstrainedNode2, Body_4);
                            NodeDirRim2->SetDirectionInAbsoluteCoords(ConstrainedNode2->D);
                            my_system.Add(NodeDirRim2);
                            break;
                        case 3:
                            ConstrainedNode = ChSharedPtr<ChNodeFEAxyzD>(
                                TireMeshes[TireNo]->GetNode(ConstNodeInx).DynamicCastTo<ChNodeFEAxyzD>());
                            NodePosRim = ChSharedPtr<ChLinkPointFrame>(new ChLinkPointFrame);
                            NodePosRim->Initialize(ConstrainedNode, Body_5);
                            my_system.Add(NodePosRim);

                            NodeDirRim = ChSharedPtr<ChLinkDirFrame>(new ChLinkDirFrame);
                            NodeDirRim->Initialize(ConstrainedNode, Body_5);
                            NodeDirRim->SetDirectionInAbsoluteCoords(ConstrainedNode->D);
                            my_system.Add(NodeDirRim);

                            // Second node to be constrained (other side of the rim)
                            ConstrainedNode2 = ChSharedPtr<ChNodeFEAxyzD>(TireMeshes[TireNo]
                                ->GetNode(ConstNodeInx + N_Thread - 1)
                                .DynamicCastTo<ChNodeFEAxyzD>());

                            NodePosRim2 = ChSharedPtr<ChLinkPointFrame>(new ChLinkPointFrame);
                            NodePosRim2->Initialize(ConstrainedNode2, Body_5);
                            my_system.Add(NodePosRim2);

                            NodeDirRim2 = ChSharedPtr<ChLinkDirFrame>(new ChLinkDirFrame);
                            NodeDirRim2->Initialize(ConstrainedNode2, Body_5);
                            NodeDirRim2->SetDirectionInAbsoluteCoords(ConstrainedNode2->D);
                            my_system.Add(NodeDirRim2);
                            break;
                        }

                    }
                }
            }  // Loop over all 4 tires
        }

        if (addGroundForces) {
            for (int TireNo = 0; TireNo < 4; TireNo++) {
                // Select on which nodes we are going to apply a load
                std::vector<ChSharedPtr<ChLoadable>> NodeList;
                for (int iNode = 0; iNode < TotalNumNodes; iNode++) {
                    ChSharedPtr<ChNodeFEAxyzD> NodeLoad(
                        TireMeshes[TireNo]->GetNode(iNode).DynamicCastTo<ChNodeFEAxyzD>());
                    NodeList.push_back(NodeLoad);
                }
                // Instance load object. This require a list of ChLoadable objects
                // (these are our two nodes,pay attention to the sequence order), and add to
                // container.
                ChSharedPtr<MyLoadCustomMultiple> Mloadcustommultiple(new MyLoadCustomMultiple(NodeList));
                Mloadcontainer4[TireNo]->Add(Mloadcustommultiple);
                my_system.Add(Mloadcontainer4[TireNo]);

                // Switch off mesh class gravity: TireMeshes[TireNo] still does not implement this element's gravity
                // forces
                TireMeshes[TireNo]->SetAutomaticGravity(false);
                // Remember to add the mesh to the system
                my_system.Add(TireMeshes[TireNo]);
            }  // End loop over tires
        }
    }

    // This is mandatory
    my_system.SetupInitial();
    ChLcpMklSolver* mkl_solver_stab = new ChLcpMklSolver;  // MKL Solver option
    ChLcpMklSolver* mkl_solver_speed = new ChLcpMklSolver;
    my_system.ChangeLcpSolverStab(mkl_solver_stab);
    my_system.ChangeLcpSolverSpeed(mkl_solver_speed);
    mkl_solver_speed->SetSparsityPatternLock(true);
    mkl_solver_stab->SetSparsityPatternLock(true);

    // INT_HHT or INT_EULER_IMPLICIT
    my_system.SetIntegrationType(ChSystem::INT_HHT);

    ChSharedPtr<ChTimestepperHHT> mystepper = my_system.GetTimestepper().DynamicCastTo<ChTimestepperHHT>();
    mystepper->SetAlpha(-0.2);  // Important for convergence
    mystepper->SetMaxiters(5);
    mystepper->SetTolerance(5e-05);
    mystepper->SetMode(ChTimestepperHHT::POSITION);
    mystepper->SetScaling(true);  //
    mystepper->SetVerbose(true);

    ChMatrixNM<double, 3, 1> Cp;
    ChMatrixNM<double, 2, 1> Cd;  // Matrices for storing constraint violations

    // Visualization
    Body_2->SetPos_dt(ChVector<>(0, 0, ForVel));
    ChSharedPtr<ChObjShapeFile> mobjmesh(new ChObjShapeFile);
    mobjmesh->SetFilename(GetChronoDataFile("fea/tractor_wheel_rim.obj"));
    Body_2->AddAsset(mobjmesh);
    double start = std::clock();
    if (showVisual) {
        for (int TireNo = 0; TireNo < 4; TireNo++) {  // Add visualization to all meshes
            ChSharedPtr<ChVisualizationFEAmesh> mvisualizemeshC(
                new ChVisualizationFEAmesh(*(TireMeshes[TireNo].get_ptr())));
            mvisualizemeshC->SetFEMglyphType(ChVisualizationFEAmesh::E_GLYPH_NODE_DOT_POS);
            mvisualizemeshC->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NONE);
            mvisualizemeshC->SetSymbolsThickness(0.003);
            TireMeshes[TireNo]->AddAsset(mvisualizemeshC);

            ChSharedPtr<ChVisualizationFEAmesh> mvisualizemeshwire(
                new ChVisualizationFEAmesh(*(TireMeshes[TireNo].get_ptr())));
            mvisualizemeshwire->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NONE);
            mvisualizemeshwire->SetWireframe(true);
            TireMeshes[TireNo]->AddAsset(mvisualizemeshwire);

            ChSharedPtr<ChVisualizationFEAmesh> mvisualizemesh(
                new ChVisualizationFEAmesh(*(TireMeshes[TireNo].get_ptr())));
            mvisualizemesh->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NODE_SPEED_NORM);
            mvisualizemesh->SetColorscaleMinMax(0.0, 30);
            mvisualizemesh->SetSmoothFaces(true);
            TireMeshes[TireNo]->AddAsset(mvisualizemesh);
            application.AssetBindAll();
            application.AssetUpdateAll();

        }  // End of add visualization to all meshes
        video::ITexture* cubeMap = application.GetVideoDriver()->getTexture(GetChronoDataFile("concrete.jpg").c_str());
        video::ITexture* rockMap = application.GetVideoDriver()->getTexture(GetChronoDataFile("rock.jpg").c_str());
        // Create the a plane using body of 'box' type:
        ChBodySceneNode* mrigidBody;
        mrigidBody = (ChBodySceneNode*)addChBodySceneNode_easyBox(
            &my_system, application.GetSceneManager(), 100.0, ChVector<>(0, GroundLoc, 0), ChQuaternion<>(1, 0, 0, 0),
            ChVector<>(10, 0.000001, 10));
        mrigidBody->GetBody()->SetBodyFixed(true);
        mrigidBody->GetBody()->GetMaterialSurface()->SetFriction(0.5);
        mrigidBody->setMaterialTexture(0, cubeMap);

        my_system.Setup();
        my_system.Update();

        GetLog() << "\n\nREADME\n\n"
                 << " - Press SPACE to start dynamic simulation \n - Press F10 for nonlinear statics - Press F11 for "
                    "linear statics. \n";

        // at beginning, no analysis is running..
        application.SetPaused(true);
        int AccuNoIterations = 0;
        application.SetStepManage(true);
        application.SetTimestep(time_step);
        application.SetTryRealtime(true);
        double ChTime = 0.0;

        utils::CSV_writer out("\t");
        out.stream().setf(std::ios::scientific | std::ios::showpos);
        out.stream().precision(7);
        const double VerForce = 0;
        const double HorForce = 4000;
        const double tini = 0.1;
        const double tend = 0.2;
        const double interval = tend - tini;
        while (application.GetDevice()->run()) {
            Body_2->Empty_forces_accumulators();
            // Body_2->Set_Scr_force(const ChVector<>& mf) { Scr_force = mf; }
            if (my_system.GetChTime() >= tini && my_system.GetChTime() <= tend) {
                Body_2->Set_Scr_torque(ChVector<>(HorForce * (my_system.GetChTime() - tini) / interval, 0, 0));
            } else if (my_system.GetChTime() > tend) {
                Body_2->Set_Scr_torque(ChVector<>(HorForce, 0, 0));
            }
            application.BeginScene();
            application.DrawAll();
            application.DoStep();
            application.EndScene();
            if (!application.GetPaused()) {
                std::cout << "Time t = " << my_system.GetChTime() << "s \n";
                AccuNoIterations += mystepper->GetNumIterations();
                printf("Forward position of rim X:      %12.4e ", Body_2->coord.pos.x);
                out << my_system.GetChTime() << Body_2->GetPos().x << Body_2->GetPos().y << Body_2->GetPos().z
                    << Body_3->GetPos().x << Body_3->GetPos().y << Body_3->GetPos().z << Body_4->GetPos().x
                    << Body_4->GetPos().y << Body_4->GetPos().z << std::endl;
                out.write_to_file("../VertPosRim.txt");
            }
        }
        double duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
        chrono::GetLog() << "Computation Time: " << duration;
        system("pause");
    } else {
        for (unsigned int it = 0; it < num_steps; it++) {
            my_system.DoStepDynamics(time_step);
            std::cout << "Time t = " << my_system.GetChTime() << "s \n";
            // std::cout << "nodetip->pos.z = " << nodetip->pos.z << "\n";
            std::cout << "mystepper->GetNumIterations()= " << mystepper->GetNumIterations() << "\n";
            if (addConstRim) {
                Cp = NodePosRim->GetC();
                printf("Point constraint violations:      %12.4e  %12.4e  %12.4e\n", Cp.GetElement(0, 0),
                       Cp.GetElement(1, 0), Cp.GetElement(2, 0));
                Cd = NodeDirRim->GetC();
                printf("Direction constraint violations:  %12.4e  %12.4e\n", Cd.GetElement(0, 0), Cd.GetElement(1, 0));

                printf("Vertical position of the rim:  %12.4e m\n", Body_2->coord.pos.y);
            }
        }
        double duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
        chrono::GetLog() << "Computation Time: " << duration;
        system("pause");
    }

    return 0;
}
