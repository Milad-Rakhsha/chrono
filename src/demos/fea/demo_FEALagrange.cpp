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

bool addConstRim = true;  //
bool addBodies = true;
bool addGroundForces = true;
bool showVisual = true;
bool addSecondMesh = true;

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
ChSharedPtr<ChLinkPointFrame> NodePosRim3;
ChSharedPtr<ChNodeFEAxyzD> ConstrainedNode3;
ChSharedPtr<ChLinkPointFrame> NodePosRim4;
ChSharedPtr<ChNodeFEAxyzD> ConstrainedNode4;
ChSharedPtr<ChLinkPointFrame> NodePosRimA;
ChSharedPtr<ChLinkPointFrame> NodePosRimB;
ChSharedPtr<ChLinkPointFrame> NodePosRimC;
ChSharedPtr<ChLinkPointFrame> NodePosRimD;

ChSharedPtr<ChNodeFEAxyzD> ConstrainedNodeA;
ChSharedPtr<ChNodeFEAxyzD> ConstrainedNodeB;
ChSharedPtr<ChNodeFEAxyzD> ConstrainedNodeC;
ChSharedPtr<ChNodeFEAxyzD> ConstrainedNodeD;
// Some model parameters

double GroundLoc = -0.11;
const int num_steps = 1500;  // Number of time steps for unit test (range 1 to 4000) 550
double time_step = 0.0002;   // Time step: 0.001 It works like that, but rotates at the beginning
// 0.008
// Specification of the mesh
// ChLoadCustomMultiple to include basic node-Ground contact interaction
class MyLoadCustomMultiple : public ChLoadCustomMultiple {
  public:
    MyLoadCustomMultiple(std::vector<ChSharedPtr<ChLoadable>>& mloadables) : ChLoadCustomMultiple(mloadables){};
    virtual void ComputeQ(ChState* state_x,      ///< state position to evaluate Q
                          ChStateDelta* state_w  ///< state speed to evaluate Q
                          ) {
        std::vector<ChSharedPtr<ChLoadable>> NodeList;
        ChVector<> Node1_Pos;
        ChVector<> Node1_Vel;
        ChVector<> Node1_Grad;
        ChVector<> Node1_GradVel;
        // this->load_Q.FillElem(0);
        double KGround = 1e4;
        double CGround = 0.001 * KGround;
        double NormalForceNode = 0;
        double FrictionCoeff = 0.0;

        // Calculation of number of nodes in contact with the Ground
        KGround = 1e4;
        CGround = 0.001 * KGround;
        GroundLoc = -0.11;
        if (state_x && state_w) {
            for (int iiinode = 0; iiinode < 8; iiinode++){
                Node1_Pos = state_x->ClipVector(iiinode * 6, 0);
                Node1_Vel = state_w->ClipVector(iiinode * 6, 0);
                if (Node1_Pos.y < GroundLoc) {
                    double Penet = abs(Node1_Pos.y - GroundLoc); // +CGround*abs(Node1_Vel.y);
                    NormalForceNode = KGround * Penet;
                this->load_Q(iiinode * 6 + 1) = NormalForceNode; }             // +CGround * abs(Node1_Vel.y)*Penet;
            }
        }
    }
    virtual bool IsStiff() { return false; }
};

int main(int argc, char* argv[]) {
    // Definition of the model
    ChSystem my_system;
    std::vector<ChSharedPtr<ChMesh>> my_mesh(2);

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


    my_mesh[0] = ChSharedPtr<ChMesh>(new ChMesh);
    my_mesh[1] = ChSharedPtr<ChMesh>(new ChMesh);

    // Create one ANCF shell element nodes
    ChSharedPtr<ChNodeFEAxyzD> node1(new ChNodeFEAxyzD(ChVector<>(-0.5, -0.1, -0.5), ChVector<>(0, 1, 0)));
    node1->SetMass(0);
    my_mesh[0]->AddNode(node1);

    ChSharedPtr<ChNodeFEAxyzD> node2(new ChNodeFEAxyzD(ChVector<>(-0.5, -0.1, 0.5), ChVector<>(0, 1, 0)));
    node2->SetMass(0);
    my_mesh[0]->AddNode(node2);

    ChSharedPtr<ChNodeFEAxyzD> node3(new ChNodeFEAxyzD(ChVector<>(0.5, -0.1, 0.5), ChVector<>(0, 1, 0)));
    node3->SetMass(0);
    my_mesh[0]->AddNode(node3);

    ChSharedPtr<ChNodeFEAxyzD> node4(new ChNodeFEAxyzD(ChVector<>(0.5, -0.1, -0.5), ChVector<>(0, 1, 0)));
    node4->SetMass(0);
    my_mesh[0]->AddNode(node4);

    if (addSecondMesh){
        ChSharedPtr<ChNodeFEAxyzD> node1_1(new ChNodeFEAxyzD(ChVector<>(-0.5, -0.1, -0.5), ChVector<>(0, 1, 0)));
        node1_1->SetMass(0);
        my_mesh[1]->AddNode(node1_1);

        ChSharedPtr<ChNodeFEAxyzD> node2_1(new ChNodeFEAxyzD(ChVector<>(-0.5, -0.1, 0.5), ChVector<>(0, 1, 0)));
        node2_1->SetMass(0);
        my_mesh[1]->AddNode(node2_1);

        ChSharedPtr<ChNodeFEAxyzD> node3_1(new ChNodeFEAxyzD(ChVector<>(0.5, -0.1, 0.5), ChVector<>(0, 1, 0)));
        node3->SetMass(0);
        my_mesh[1]->AddNode(node3_1);

        ChSharedPtr<ChNodeFEAxyzD> node4_1(new ChNodeFEAxyzD(ChVector<>(0.5, -0.1, -0.5), ChVector<>(0, 1, 0)));
        node4_1->SetMass(0);
        my_mesh[1]->AddNode(node4_1);
    }


    ChSharedPtr<ChMaterialShellANCF> mat(new ChMaterialShellANCF(500, 5.0e5, 0.3));

    // Create the element and set its nodes.
    ChSharedPtr<ChElementShellANCF> element(new ChElementShellANCF);
    ChSharedPtr<ChElementShellANCF> element1(new ChElementShellANCF);
    element->SetNodes(
        my_mesh[0]->GetNode(0).DynamicCastTo<ChNodeFEAxyzD>(), my_mesh[0]->GetNode(1).DynamicCastTo<ChNodeFEAxyzD>(),
        my_mesh[0]->GetNode(3).DynamicCastTo<ChNodeFEAxyzD>(), my_mesh[0]->GetNode(2).DynamicCastTo<ChNodeFEAxyzD>());
    if (addSecondMesh){
        element1->SetNodes(
            my_mesh[1]->GetNode(0).DynamicCastTo<ChNodeFEAxyzD>(), my_mesh[1]->GetNode(1).DynamicCastTo<ChNodeFEAxyzD>(),
            my_mesh[1]->GetNode(3).DynamicCastTo<ChNodeFEAxyzD>(), my_mesh[1]->GetNode(2).DynamicCastTo<ChNodeFEAxyzD>());
        element1->SetDimensions(1.0, 1.0);

        // Add a single layers with a fiber angle of 0 degrees.
        element1->AddLayer(0.01, 0.0 * CH_C_DEG_TO_RAD, mat);

        // Set other element properties
        element1->SetAlphaDamp(0.00);  // Structural damping for this element
        element1->SetGravityOn(true);  // gravitational forces
        my_mesh[1]->AddElement(element);
        my_mesh[1]->SetAutomaticGravity(false);
        my_system.Add(my_mesh[1]);
    }
    // Set element dimensions
    element->SetDimensions(1.0, 1.0);

    // Add a single layers with a fiber angle of 0 degrees.
    element->AddLayer(0.01, 0.0 * CH_C_DEG_TO_RAD, mat);

    // Set other element properties
    element->SetAlphaDamp(0.00);  // Structural damping for this element
    element->SetGravityOn(true);  // gravitational forces

    // Add element to mesh
    my_mesh[0]->AddElement(element);

    // Switch off mesh class gravity
    my_mesh[0]->SetAutomaticGravity(false);


    // Add tmy_system.Addhe mesh to the system
    my_system.Add(my_mesh[0]);


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
        Body_2->SetMass(10);
        Body_2->SetInertiaXX(ChVector<>(0.3, 0.3, 0.3));
        Body_2->SetPos(ChVector<>(0, 0, 0));  // Y = -1m
        // Defining the Body 2: Rim for TireMeshes[0]
        Body_3 = ChSharedPtr<ChBody>(new ChBody);
        my_system.AddBody(Body_3);
        Body_3->SetIdentifier(3);
        Body_3->SetBodyFixed(false);
        Body_3->SetCollide(false);
        Body_3->SetMass(10);
        Body_3->SetInertiaXX(ChVector<>(0.3, 0.3, 0.3));
        Body_3->SetPos(ChVector<>(0, 0, 0));  // Y = -1m
        my_system.Set_G_acc(ChVector<>(0, -9.81, 0));
    }

    ChSharedPtr<ChLoadContainer> MloadcontainerGround(new ChLoadContainer);

    if (addConstRim) {
        ConstrainedNode =
            ChSharedPtr<ChNodeFEAxyzD>(my_mesh[0]->GetNode(0).DynamicCastTo<ChNodeFEAxyzD>());
        NodePosRim = ChSharedPtr<ChLinkPointFrame>(new ChLinkPointFrame);
        NodePosRim->Initialize(ConstrainedNode, Body_2);
        my_system.Add(NodePosRim);
        // GetLog() << "ConstrainedNode: " << ConstNodeInx << "\n";

        // Second node to be constrained (other side of the rim)
        ConstrainedNode2 = ChSharedPtr<ChNodeFEAxyzD>(
            my_mesh[0]->GetNode(1).DynamicCastTo<ChNodeFEAxyzD>());

        NodePosRim2 = ChSharedPtr<ChLinkPointFrame>(new ChLinkPointFrame);
        NodePosRim2->Initialize(ConstrainedNode2, Body_2);
        my_system.Add(NodePosRim2);

        // Second node to be constrained (other side of the rim)
        ConstrainedNode3 = ChSharedPtr<ChNodeFEAxyzD>(
            my_mesh[0]->GetNode(2).DynamicCastTo<ChNodeFEAxyzD>());

        NodePosRim3 = ChSharedPtr<ChLinkPointFrame>(new ChLinkPointFrame);
        NodePosRim3->Initialize(ConstrainedNode3, Body_2);
        my_system.Add(NodePosRim3);

        // Second node to be constrained (other side of the rim)
        ConstrainedNode4 = ChSharedPtr<ChNodeFEAxyzD>(
            my_mesh[0]->GetNode(3).DynamicCastTo<ChNodeFEAxyzD>());

        NodePosRim4 = ChSharedPtr<ChLinkPointFrame>(new ChLinkPointFrame);
        NodePosRim4->Initialize(ConstrainedNode4, Body_2);
        my_system.Add(NodePosRim4);

        if (addSecondMesh){
        ConstrainedNodeA =
            ChSharedPtr<ChNodeFEAxyzD>(my_mesh[1]->GetNode(0).DynamicCastTo<ChNodeFEAxyzD>());
        NodePosRimA = ChSharedPtr<ChLinkPointFrame>(new ChLinkPointFrame);
        NodePosRimA->Initialize(ConstrainedNodeA, Body_3);
        my_system.Add(NodePosRimA);
        // GetLog() << "ConstrainedNode: " << ConstNodeInx << "\n";

        // Second node to be constrained (other side of the rim)
        ConstrainedNodeB = ChSharedPtr<ChNodeFEAxyzD>(
            my_mesh[1]->GetNode(1).DynamicCastTo<ChNodeFEAxyzD>());

        NodePosRimB = ChSharedPtr<ChLinkPointFrame>(new ChLinkPointFrame);
        NodePosRimB->Initialize(ConstrainedNodeB, Body_3);
        my_system.Add(NodePosRimB);

        // Second node to be constrained (other side of the rim)
        ConstrainedNodeC = ChSharedPtr<ChNodeFEAxyzD>(
            my_mesh[1]->GetNode(2).DynamicCastTo<ChNodeFEAxyzD>());

        NodePosRimC = ChSharedPtr<ChLinkPointFrame>(new ChLinkPointFrame);
        NodePosRimC->Initialize(ConstrainedNodeC, Body_3);
        my_system.Add(NodePosRimC);

        // Second node to be constrained (other side of the rim)
        ConstrainedNodeD = ChSharedPtr<ChNodeFEAxyzD>(
            my_mesh[1]->GetNode(3).DynamicCastTo<ChNodeFEAxyzD>());

        NodePosRimD = ChSharedPtr<ChLinkPointFrame>(new ChLinkPointFrame);
        NodePosRimD->Initialize(ConstrainedNodeD, Body_3);
        my_system.Add(NodePosRimD);
        }

    }

    if (addGroundForces) {
        // Select on which nodes we are going to apply a load
        std::vector<ChSharedPtr<ChLoadable>> NodeList;
        for (int iNode = 0; iNode < 4; iNode++) {
                ChSharedPtr<ChNodeFEAxyzD> NodeLoad(my_mesh[0]->GetNode(iNode).DynamicCastTo<ChNodeFEAxyzD>());
                NodeList.push_back(NodeLoad);
        }
        for (int iNode = 0; iNode < 4; iNode++) {
            if (addSecondMesh){
                ChSharedPtr<ChNodeFEAxyzD> NodeLoad1(my_mesh[1]->GetNode(iNode).DynamicCastTo<ChNodeFEAxyzD>());
                NodeList.push_back(NodeLoad1);
            }
        }


        ChSharedPtr<MyLoadCustomMultiple> Mloadcustommultiple(new MyLoadCustomMultiple(NodeList));
        MloadcontainerGround->Add(Mloadcustommultiple);
    }  // End loop over tires

    ChSharedPtr<ChNodeFEAxyzD> NodeMesh0(
        my_mesh[0]->GetNode(0).DynamicCastTo<ChNodeFEAxyzD>());

    ChSharedPtr<ChNodeFEAxyzD> NodeMesh1(
        my_mesh[1]->GetNode(0).DynamicCastTo<ChNodeFEAxyzD>());

    my_system.Add(MloadcontainerGround);
    // Mark completion of system construction
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
    mystepper->SetAlpha(-0.3);  // Important for convergence
    mystepper->SetMaxiters(15);
    mystepper->SetTolerance(1e-6);
    mystepper->SetMode(ChTimestepperHHT::POSITION);
    mystepper->SetScaling(true);  //
    mystepper->SetVerbose(true);

    ChMatrixNM<double, 3, 1> Cp;
    ChMatrixNM<double, 2, 1> Cd;  // Matrices for storing constraint violations

    double start = std::clock();
    if (showVisual) {
            ChSharedPtr<ChVisualizationFEAmesh> mvisualizemeshC(
                new ChVisualizationFEAmesh(*(my_mesh[0].get_ptr())));
            mvisualizemeshC->SetFEMglyphType(ChVisualizationFEAmesh::E_GLYPH_NODE_DOT_POS);
            mvisualizemeshC->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NONE);
            mvisualizemeshC->SetSymbolsThickness(0.003);
            my_mesh[0]->AddAsset(mvisualizemeshC);

            ChSharedPtr<ChVisualizationFEAmesh> mvisualizemeshwire(
                new ChVisualizationFEAmesh(*(my_mesh[0].get_ptr())));
            mvisualizemeshwire->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NONE);
            mvisualizemeshwire->SetWireframe(true);
            my_mesh[0]->AddAsset(mvisualizemeshwire);

            ChSharedPtr<ChVisualizationFEAmesh> mvisualizemesh(
                new ChVisualizationFEAmesh(*(my_mesh[0].get_ptr())));
            mvisualizemesh->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NODE_SPEED_NORM);
            mvisualizemesh->SetColorscaleMinMax(0.0, 30);
            mvisualizemesh->SetSmoothFaces(true);
            my_mesh[0]->AddAsset(mvisualizemesh);

        application.AssetBindAll();
        application.AssetUpdateAll();

        video::ITexture* cubeMap = application.GetVideoDriver()->getTexture(GetChronoDataFile("concrete.jpg").c_str());
        video::ITexture* rockMap = application.GetVideoDriver()->getTexture(GetChronoDataFile("rock.jpg").c_str());
        // Create the a plane using body of 'box' type:
        ChBodySceneNode* mrigidBody;
        mrigidBody = (ChBodySceneNode*)addChBodySceneNode_easyBox(
            &my_system, application.GetSceneManager(), 100.0, ChVector<>(0, GroundLoc, 0), ChQuaternion<>(1, 0, 0, 0),
            ChVector<>(10, 0.000001, 10));
        mrigidBody->GetBody()->SetBodyFixed(true);
        mrigidBody->GetBody()->GetMaterialSurface()->SetFriction(0.0);
        mrigidBody->setMaterialTexture(0, cubeMap);
        mrigidBody->GetBody()->SetCollide(false);
        my_system.Setup();
        my_system.Update();

        chrono::GetLog()
            << "\n\nREADME\n\n"
            << " - Press SPACE to start dynamic simulation \n - Press F10 for nonlinear statics - Press F11 for "
               "linear statics. \n";

        // at beginning, no analysis is running..
        application.SetPaused(false);
        int AccuNoIterations = 0;
        application.SetStepManage(true);
        application.SetTimestep(time_step);
        application.SetTryRealtime(false);
        double ChTime = 0.0;

        utils::CSV_writer out("\t");
        out.stream().setf(std::ios::scientific | std::ios::showpos);
        out.stream().precision(7);
        const double VerForce = 0;
        const double HorForce = 4000;
        const double tini = 0.02;
        const double tend = 0.2;
        const double interval = tend - tini;
        while (application.GetDevice()->run()) {
            application.BeginScene();
            application.DrawAll();
            application.DoStep();
            application.EndScene();
            if (!application.GetPaused()) {
                std::cout << "Time t = " << my_system.GetChTime() << "s \n";
                printf(
                    "Vertical position of Nodes:      %12.4e    %12.4e \t Body  %12.4e %12.4e \t Lateral Pos Bodies %12.4e %12.4e",
                    NodeMesh0->pos.y, NodeMesh1->pos.y, Body_2->GetPos().y, Body_3->GetPos().y, Body_2->GetPos().z, Body_3->GetPos().z);
            }
        }
        double duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
        chrono::GetLog() << "Computation Time: " << duration;
        system("pause");
    } else {
        GetLog() << "error: Finish simulation \n";
        system("pause");
    }

    return 0;
}
