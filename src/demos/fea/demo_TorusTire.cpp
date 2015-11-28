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
// tire/ground interaction. The tire is made up of ANCF shell elements.
// The number of ANCF shell elements is user-defined.
// Gravity acc. vec. of the ANCF shell element must be modified for it to act along the Y axis.
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
#include "physics/ChLoadContainer.h"   // For Alessandro classes
  // For visualization
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
bool showVisual = true; // Also adds pushing force
bool addSingleLoad = false;

ChSharedPtr<ChBody> BGround;
ChSharedPtr<ChBody> Body_2;
ChSharedPtr<ChLinkPointFrame> NodePosRim;
ChSharedPtr<ChLinkDirFrame> NodeDirRim;
ChSharedPtr<ChNodeFEAxyzD> ConstrainedNode;
ChSharedPtr<ChLinkPointFrame> NodePosRim2;
ChSharedPtr<ChLinkDirFrame> NodeDirRim2;
ChSharedPtr<ChNodeFEAxyzD> ConstrainedNode2;

// Some model parameters

const int num_steps = 40000;     // Number of time steps for unit test (range 1 to 4000)
const double time_step = 0.0002;  // Time step
const double GroundLoc = -0.35; // Where the ground is 
const double ForVel = 2.0; // Initial velocity of the tire
// Geometry of the plate
// double plate_lenght_x = 1; // The length along X axis is defined parametrically in
// the initial coordinate section
// double plate_lenght_y = 1;
double plate_lenght_z = 0.005;  // Thickness of EACH layer (number of layers defined below)

// Specification of the mesh
const int numEl_Diameter = 24;  // Number of elements along diameter
const int numEl_Thread = 6;     // Number of elements along thread
const int numEl_z = 1;          // Single element along the thickness
const int N_Diameter = numEl_Diameter;
const int N_Thread = numEl_Thread + 1;
const int N_z = numEl_z + 1;
const double TorusSmallRadius = 0.050;
const double TorusRadius = 0.3;  // Radius of the cylinder
int TotalNumElements = numEl_Diameter * numEl_Thread;
int TotalNumNodes = (numEl_Diameter) * (numEl_Thread + 1);

double dz = plate_lenght_z / numEl_z;  // We are going to calculate distances through coordinates
double dx = 2 * CH_C_PI * (TorusRadius + TorusSmallRadius) / 2 / N_Diameter;
double dy = CH_C_PI * TorusSmallRadius / numEl_Thread;
// Ground Forces through ChLoadCustomMultiple
// Custom load n.1: a class that loads 4 ChNodeFEAxyzD
// with a stiff load - for testing.

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
        const double KGround = 1e6;
        double NormalForceNode = 0;
        double FrictionCoeff = 6.0;  // F1-like
        if (state_x && state_w) {
            for (unsigned int iii = 0; iii < TotalNumNodes; iii++) {
                Node1_Pos = state_x->ClipVector(iii * 6, 0);
                Node1_Grad = state_x->ClipVector(iii * 6 + 3, 0);
                Node1_Vel = state_w->ClipVector(iii * 6, 0);
                Node1_GradVel = state_w->ClipVector(iii * 6 + 3, 0);

                if (Node1_Pos.y < GroundLoc) {
                    NormalForceNode = KGround * abs(Node1_Pos.y - GroundLoc);

					 
                    this->load_Q(iii * 6 + 1) = NormalForceNode;  // Fy (Vertical)
                    this->load_Q(iii * 6 + 0) =
                        -NormalForceNode * FrictionCoeff *
                        (Node1_Vel.x / sqrt((pow(Node1_Vel.x, 2) + pow(Node1_Vel.z, 2))));  // Fx (Plane x)
                    this->load_Q(iii * 6 + 2) =
                        -NormalForceNode * FrictionCoeff *
                        (Node1_Vel.z / sqrt((pow(Node1_Vel.x, 2) + pow(Node1_Vel.z, 2))));  // Fz (Plane y)
					//printf("Node Number:      %4d  \n", iii);
					//printf("Interpenetration:      %12.4e  m\n", abs(Node1_Pos.y - GroundLoc));
					//printf("Long. Vel:      %12.4e  N\n", Node1_Vel.x);
                }
                Node1_Pos.SetNull();
            }

        } else {
            // explicit integrators might call ComputeQ(0,0), null pointers mean
            // that we assume current state, without passing state_x for efficiency
        }

        // store generalized forces as a contiguous vector in this->load_Q, with
        // same order of state_w
    }

    // OPTIONAL: if you want to provide an analytical jacobian, just implement the
    // following:
    //   virtual void ComputeJacobian(...)

    // Remember to set this as stiff, to enable the jacobians
    virtual bool IsStiff() { return true; }
};

int main(int argc, char* argv[]) {
    // Definition of the model
    ChSystem my_system;
    ChSharedPtr<ChMesh> my_mesh(new ChMesh);

    //
    // bind a simple user interface, etc. etc.)
    ChIrrApp application(&my_system, L"ANCF Rolling Tire", core::dimension2d<u32>(800, 600), false);
    // Easy shortcuts to add camera, lights, logo and sky in Irrlicht scene:
    application.AddTypicalLogo();
    application.AddTypicalSky();
    application.AddTypicalLights();
    application.AddTypicalCamera(core::vector3df(0.5f, 0.5f, 1.15f),   // camera location
                                 core::vector3df(0.65f, 0.0f, 0.0f));  // "look at" location
    //

    int numFlexBody = 1;

    utils::CSV_writer out("\t");
    out.stream().setf(std::ios::scientific | std::ios::showpos);
    out.stream().precision(7);

    int MaxMNUM = 0;
    int MTYPE = 0;
    int MaxLayNum = 0;
    ChMatrixDynamic<double> COORDFlex(TotalNumNodes, 6);
    ChMatrixDynamic<double> VELCYFlex(TotalNumNodes, 6);
    ChMatrixDynamic<int> NumNodes(TotalNumElements, 4);
    ChMatrixDynamic<int> LayNum(TotalNumElements, 1);
    ChMatrixDynamic<int> NDR(TotalNumNodes, 6);
    ChMatrixDynamic<double> ElemLengthXY(TotalNumElements, 2);
    ChMatrixNM<double, 10, 12> MPROP;
    ChMatrixNM<int, 10, 7> MNUM;
    ChMatrixNM<int, 10, 1> NumLayer;
    double LayPROP[10][7][2];

    //--------------- Element data--------------------
    for (int j = 0; j < N_Diameter; j++) {  // Start node numbering by zero
        for (int i = 0; i < numEl_Thread; i++) {
            LayNum(i + j * numEl_Thread, 0) = 1;  // Each element has one layer
            if (j == N_Diameter - 1) {
                NumNodes(i + j * numEl_Thread, 0) = i + j * (numEl_Thread + 1);
                NumNodes(i + j * numEl_Thread, 1) = i;
                NumNodes(i + j * numEl_Thread, 2) = i + 1;
                NumNodes(i + j * numEl_Thread, 3) = i + 1 + j * (numEl_Thread + 1);
            } else {
                NumNodes(i + j * numEl_Thread, 0) = i + j * (numEl_Thread + 1);
                NumNodes(i + j * numEl_Thread, 1) = i + (j + 1) * (numEl_Thread + 1);
                NumNodes(i + j * numEl_Thread, 2) = i + 1 + (j + 1) * (numEl_Thread + 1);
                NumNodes(i + j * numEl_Thread, 3) = i + 1 + j * (numEl_Thread + 1);
            }

            // Add shell dimensions based on coordinates
            ElemLengthXY(i + j * numEl_Thread, 0) = dx;  // dx
            ElemLengthXY(i + j * numEl_Thread, 1) = dy;  // dy

            if (MaxLayNum < LayNum(i, 0)) {
                MaxLayNum = LayNum(i, 0);
            }
        }
    }
    //  Fixing constraints, initial coordinates and velocities
    double thethaTorus = 0.0;
    double phiTorus = 0.0;

    for (int j = 0; j < N_Diameter; j++) {
        for (int i = 0; i < N_Thread; i++) {
            thethaTorus = CH_C_PI / 2 + CH_C_PI * i / (N_Thread - 1);
            thethaTorus = -CH_C_PI / 2 + CH_C_PI * i / (N_Thread - 1);
            phiTorus = 2 * CH_C_PI * j / N_Diameter;

            COORDFlex(i + j * N_Thread, 0) = (TorusRadius + TorusSmallRadius * cos(thethaTorus)) * cos(phiTorus);
            COORDFlex(i + j * N_Thread, 1) = (TorusRadius + TorusSmallRadius * cos(thethaTorus)) * sin(phiTorus);
            COORDFlex(i + j * N_Thread, 2) = TorusSmallRadius * sin(thethaTorus);

            COORDFlex(i + j * N_Thread, 3) = cos(thethaTorus) * cos(phiTorus);
            COORDFlex(i + j * N_Thread, 4) = cos(thethaTorus) * sin(phiTorus);
            COORDFlex(i + j * N_Thread, 5) = sin(thethaTorus);

            // Write position and gradient of nodes to a file and plot them in Matlab
            out << i + j * N_Thread << COORDFlex(i + j * N_Thread, 0) << COORDFlex(i + j * N_Thread, 1)
                << COORDFlex(i + j * N_Thread, 2) << COORDFlex(i + j * N_Thread, 3) << COORDFlex(i + j * N_Thread, 4)
                << COORDFlex(i + j * N_Thread, 5) << std::endl;
            out.write_to_file("../TorusTire.txt");

			VELCYFlex(i + j * N_Thread, 0) = ForVel;
            VELCYFlex(i + j * N_Thread, 1) = 0;
            VELCYFlex(i + j * N_Thread, 2) = 0;
            VELCYFlex(i + j * N_Thread, 3) = 0;
            VELCYFlex(i + j * N_Thread, 4) = 0;
            VELCYFlex(i + j * N_Thread, 5) = 0;
        }
    }
    //------------- Read Layer Data-------------------

    for (int i = 0; i < MaxLayNum; i++) {
        NumLayer(i, 0) = i + 1;

        for (int j = 0; j < NumLayer(i, 0); j++) {
            LayPROP[i][j][0] = dz;  // Height of each layer
            if (j == 0) {
                LayPROP[i][j][1] = 0;  // For first layer, fiber angle 20 degrees
            }                          // Fiber angle of each ply
            else {
                LayPROP[i][j][1] = 0;  // For second layer, fiber angle -20 degrees
            }
            MNUM[i][j] = 1;  // Material_ID
            // In this example one single material ID (same material properties for the two layers)
            if (MaxMNUM < MNUM(i, j))
                MaxMNUM = MNUM(i, j);
        }
    }
    //------------------------------------------------
    //------------ Input Material Data----------------
    //------------------------------------------------
    for (int i = 0; i < MaxMNUM; i++) {
        double nu_coef = 0.3;
        MTYPE = 2;  // The user must use orthotropic input (MTYPE=2)
        if (MTYPE == 2) {
            MPROP(i, 0) = 500;                                    // Density [kg/m3]
            MPROP(i, 1) = 6.0E+06;                                // Ex//
            MPROP(i, 2) = 6.0E+06;                                // Ey
            MPROP(i, 3) = 6.0E+06;                                // Ez
            MPROP(i, 4) = 0.3;                                    // nuxy
            MPROP(i, 5) = 0.3;                                    // nuxz
            MPROP(i, 6) = 0.3;                                    // nuyz
            MPROP(i, 7) = MPROP(i, 1) / 2.0 / (1 + MPROP(i, 6));  // Gxy
            MPROP(i, 8) = MPROP(i, 2) / 2.0 / (1 + MPROP(i, 6));  // Gxz
            MPROP(i, 9) = MPROP(i, 3) / 2.0 / (1 + MPROP(i, 6));  // Gyz
        }
    }
    ChSharedPtr<ChContinuumElastic> mmaterial(new ChContinuumElastic);
    // Adding the nodes to the mesh
    int i = 0;

    while (i < TotalNumNodes) {
        ChSharedPtr<ChNodeFEAxyzD> node(
            new ChNodeFEAxyzD(ChVector<>(COORDFlex(i, 0), COORDFlex(i, 1), COORDFlex(i, 2)),
                              ChVector<>(COORDFlex(i, 3), COORDFlex(i, 4), COORDFlex(i, 5))));
        node->SetMass(0.0);
        my_mesh->AddNode(node);
        // if (NDR(i, 0) == 1 && NDR(i, 1) == 1 && NDR(i, 2) == 1 && NDR(i, 3) == 1 && NDR(i, 4) == 1 && NDR(i, 5) == 1)
        // {
        //    node->SetFixed(true);
        //}
        i++;
    }

    int elemcount = 0;
    while (elemcount < TotalNumElements) {
        ChSharedPtr<ChElementShellANCF> element(new ChElementShellANCF);
        // Save material data into InertFlexVec(98x1) at each layer
        ChMatrixNM<double, 98, 1> InertFlexVec;
        InertFlexVec.Reset();
        double TotalThickness;  // Element thickness: Summation of all layers' thickness
        TotalThickness = 0.0;
        int i = elemcount;
        for (int j = 0; j < NumLayer(LayNum(i, 0) - 1, 0); j++) {  // For each element, define material properties
            int ij = 14 * j;
            InertFlexVec(ij) = MPROP[MNUM[LayNum(i, 0) - 1][j] - 1][0];  // Density
            InertFlexVec(ij + 1) = ElemLengthXY(i, 0);                   // EL
            InertFlexVec(ij + 2) = ElemLengthXY(i, 1);                   // EW
            InertFlexVec(ij + 3) = LayPROP[LayNum(i, 0) - 1][j][0];      // Thickness per layer
            TotalThickness += InertFlexVec(ij + 3);
            InertFlexVec(ij + 4) = LayPROP[LayNum(i, 0) - 1][j][1];           // Fiber angle
            InertFlexVec(ij + 5) = MPROP[MNUM[LayNum(i, 0) - 1][j] - 1][1];   // Ex
            InertFlexVec(ij + 6) = MPROP[MNUM[LayNum(i, 0) - 1][j] - 1][2];   // Ey
            InertFlexVec(ij + 7) = MPROP[MNUM[LayNum(i, 0) - 1][j] - 1][3];   // Ez
            InertFlexVec(ij + 8) = MPROP[MNUM[LayNum(i, 0) - 1][j] - 1][4];   // nuxy
            InertFlexVec(ij + 9) = MPROP[MNUM[LayNum(i, 0) - 1][j] - 1][5];   // nuxz
            InertFlexVec(ij + 10) = MPROP[MNUM[LayNum(i, 0) - 1][j] - 1][6];  // nuyz
            InertFlexVec(ij + 11) = MPROP[MNUM[LayNum(i, 0) - 1][j] - 1][7];  // Gxy
            InertFlexVec(ij + 12) = MPROP[MNUM[LayNum(i, 0) - 1][j] - 1][8];  // Gxz
            InertFlexVec(ij + 13) = MPROP[MNUM[LayNum(i, 0) - 1][j] - 1][9];  // Gyz
        }
        // Automatically selects the Gauss range to integrate inertia/internal forces over the thickness
        // Each layer's thickness will be integrated using 2 Gauss integration points
        // LimLow and LimHigh represent the thickness range for a given layer
        ChMatrixNM<double, 7, 2> GaussZRange;
        GaussZRange.Reset();
        double CurrentHeight = 0.0;
        for (int j = 0; j < NumLayer(LayNum(i, 0) - 1, 0); j++) {
            double LimLow = (CurrentHeight / TotalThickness - 0.5) * 2.0;
            CurrentHeight += LayPROP[LayNum(i, 0) - 1][j][0];
            double LimHigh = (CurrentHeight / TotalThickness - 0.5) * 2.0;
            GaussZRange(j, 0) = LimLow;
            GaussZRange(j, 1) = LimHigh;
        }
        // Now we give some parameters element by element
        element->SetInertFlexVec(InertFlexVec);
        element->SetGaussZRange(GaussZRange);
        element->SetNodes(my_mesh->GetNode(NumNodes[elemcount][0]).DynamicCastTo<ChNodeFEAxyzD>(),
                          my_mesh->GetNode(NumNodes[elemcount][1]).DynamicCastTo<ChNodeFEAxyzD>(),
                          my_mesh->GetNode(NumNodes[elemcount][3]).DynamicCastTo<ChNodeFEAxyzD>(),
                          my_mesh->GetNode(NumNodes[elemcount][2]).DynamicCastTo<ChNodeFEAxyzD>());
        element->SetMaterial(mmaterial);
        element->SetNumLayers(NumLayer(LayNum(i, 0) - 1, 0));
        element->SetThickness(TotalThickness);
        element->SetElemNum(elemcount);
        element->SetAlphaDamp(0.05);               // Structural damping for this
        element->Setdt(time_step);                 // dt to calculate DampingCoefficient
        element->SetGravityOn(true);               // turn gravity on/off
		element->SetAirPressureOn(false);
        ChMatrixNM<double, 35, 1> StockAlpha_EAS;  // StockAlpha(5*7,1): Max #Layer is 7
        StockAlpha_EAS.Reset();
        element->SetStockAlpha(StockAlpha_EAS);
        my_mesh->AddElement(element);
        elemcount++;
    }

	if (addSingleLoad)
	{
		ChSharedPtr<ChNodeFEAxyzD> nodetip(my_mesh->GetNode(TotalNumNodes - 1).DynamicCastTo<ChNodeFEAxyzD>());
		nodetip->SetForce(ChVector<>(0, 0, -10));
	}


    if (addBodies) {
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
        // Defining the Body 2
        Body_2 = ChSharedPtr<ChBody>(new ChBody);
        my_system.AddBody(Body_2);
        Body_2->SetIdentifier(2);
        Body_2->SetBodyFixed(false);
        Body_2->SetCollide(false);
        Body_2->SetMass(6);
        Body_2->SetInertiaXX(ChVector<>(1, 0.3, 0.3));
        Body_2->SetPos(ChVector<>(0, 0, 0));  // Y = -1m
        Body_2->SetRot(rot);
        my_system.Set_G_acc(ChVector<>(0, -9.81*4, 0)); // Hey! 4G!
    }

    // Add constraints to the rim
    if (addConstRim) {
        int NoConstNodes = 2 * numEl_Diameter;
        int ConstNodeInx = -1;

        for (int j = 0; j < N_Diameter; j++) {
            for (int i = 0; i < N_Thread; i++) {
                ConstNodeInx = i + j * N_Thread;  // General node index

                if ((i + j * N_Thread) % (N_Thread) == 0) {
                    // First node to be constrained (one side of the rim)
                    ConstrainedNode =
                        ChSharedPtr<ChNodeFEAxyzD>(my_mesh->GetNode(ConstNodeInx).DynamicCastTo<ChNodeFEAxyzD>());
                    NodePosRim = ChSharedPtr<ChLinkPointFrame>(new ChLinkPointFrame);
                    NodePosRim->Initialize(ConstrainedNode, Body_2);
                    my_system.Add(NodePosRim);

                    NodeDirRim = ChSharedPtr<ChLinkDirFrame>(new ChLinkDirFrame);
                    NodeDirRim->Initialize(ConstrainedNode, Body_2);
                    NodeDirRim->SetDirectionInAbsoluteCoords(ConstrainedNode->D);
                    my_system.Add(NodeDirRim);

                    // Second node to be constrained (other side of the rim)

                    ConstrainedNode2 = ChSharedPtr<ChNodeFEAxyzD>(
                        my_mesh->GetNode(ConstNodeInx + N_Thread - 1).DynamicCastTo<ChNodeFEAxyzD>());

                    NodePosRim2 = ChSharedPtr<ChLinkPointFrame>(new ChLinkPointFrame);
                    NodePosRim2->Initialize(ConstrainedNode2, Body_2);
                    my_system.Add(NodePosRim2);

                    NodeDirRim2 = ChSharedPtr<ChLinkDirFrame>(new ChLinkDirFrame);
                    NodeDirRim2->Initialize(ConstrainedNode2, Body_2);
                    NodeDirRim2->SetDirectionInAbsoluteCoords(ConstrainedNode2->D);
                    my_system.Add(NodeDirRim2);

                    // chrono::GetLog() << "We are constraining the nodes:  " << ConstNodeInx << "\t and \t"
                    //                 << ConstNodeInx + N_Thread - 1;
                    // system("pause");
                    // i = 2000; j = 2000;
                }
            }
        }

        if (addGroundForces) {  // APPLY SOME LOADS! (using ChLoad classes)
                                // Select on which nodes we are going to apply a load
            // ChSharedPtr<ChNodeFEAxyzD> NodeLoad6(my_mesh->GetNode(6).DynamicCastTo<ChNodeFEAxyzD>());
            // ChSharedPtr<ChNodeFEAxyzD> NodeLoad7(my_mesh->GetNode(7).DynamicCastTo<ChNodeFEAxyzD>());
            std::vector<ChSharedPtr<ChLoadable>> NodeList;
            for (int iNode = 0; iNode < TotalNumNodes; iNode++) {
                ChSharedPtr<ChNodeFEAxyzD> NodeLoad(my_mesh->GetNode(iNode).DynamicCastTo<ChNodeFEAxyzD>());
                NodeList.push_back(NodeLoad);
            }
            // First: loads must be added to "load containers",
            // and load containers must be added to your ChSystem
            ChSharedPtr<ChLoadContainer> Mloadcontainer(new ChLoadContainer);
            my_system.Add(Mloadcontainer);

            // Instance load object. This require a list of ChLoadable objects
            // (these are our two nodes,pay attention to the sequence order), and add to
            // container.
            ChSharedPtr<MyLoadCustomMultiple> Mloadcustommultiple(new MyLoadCustomMultiple(NodeList));
            Mloadcontainer->Add(Mloadcustommultiple);
        }
    }
    // Switch off mesh class gravity: my_mesh still does not implement this element's gravity forces
    my_mesh->SetAutomaticGravity(false);
    // This is mandatory
    my_mesh->SetupInitial();
    // Remember to add the mesh to the system
    my_system.Add(my_mesh);

    ChLcpMklSolver* mkl_solver_stab = new ChLcpMklSolver;  // MKL Solver option
    ChLcpMklSolver* mkl_solver_speed = new ChLcpMklSolver;
    my_system.ChangeLcpSolverStab(mkl_solver_stab);
    my_system.ChangeLcpSolverSpeed(mkl_solver_speed);
    mkl_solver_stab->SetSparsityPatternLock(true);
    mkl_solver_stab->SetSparsityPatternLock(false);
    my_system.Update();

    // INT_HHT or INT_EULER_IMPLICIT
    my_system.SetIntegrationType(ChSystem::INT_HHT);

    ChSharedPtr<ChTimestepperHHT> mystepper = my_system.GetTimestepper().DynamicCastTo<ChTimestepperHHT>();
    mystepper->SetAlpha(-0.2);  // Important for convergence
    mystepper->SetMaxiters(1000);
    mystepper->SetTolerance(8e-03);
    mystepper->SetMode(ChTimestepperHHT::POSITION);
    mystepper->SetScaling(true);  //

    ChMatrixNM<double, 3, 1> Cp;
    ChMatrixNM<double, 2, 1> Cd;  // Matrices for storing constraint violations

    // Visualization
	// Body_2->Accumulate_force(ChVector<>(60, 0, 0), Body_2->GetPos(), false);
    // Options for visualization in irrlicht
	Body_2->SetPos_dt(ChVector<>(ForVel,0,0));

    if (showVisual) {
        ChSharedPtr<ChVisualizationFEAmesh> mvisualizemesh(new ChVisualizationFEAmesh(*(my_mesh.get_ptr())));
        mvisualizemesh->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_SURFACE);
        mvisualizemesh->SetSmoothFaces(false);
        my_mesh->AddAsset(mvisualizemesh);

        ChSharedPtr<ChVisualizationFEAmesh> mvisualizemeshC(new ChVisualizationFEAmesh(*(my_mesh.get_ptr())));
        mvisualizemeshC->SetFEMglyphType(ChVisualizationFEAmesh::E_GLYPH_NODE_DOT_POS);
        mvisualizemeshC->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NONE);
        mvisualizemeshC->SetSymbolsThickness(0.004);
        my_mesh->AddAsset(mvisualizemeshC);

        ChSharedPtr<ChVisualizationFEAmesh> mvisualizemeshwire(new ChVisualizationFEAmesh(*(my_mesh.get_ptr())));
        mvisualizemeshwire->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_SURFACE);
        mvisualizemeshwire->SetWireframe(true);
        my_mesh->AddAsset(mvisualizemeshwire);
        application.AssetBindAll();
        application.AssetUpdateAll();

		video::ITexture* cubeMap = application.GetVideoDriver()->getTexture(GetChronoDataFile("concrete.jpg").c_str());
		video::ITexture* rockMap = application.GetVideoDriver()->getTexture(GetChronoDataFile("rock.jpg").c_str());
		// Create the a plane using body of 'box' type:
		///ChBodySceneNode* mrigidBody;
		///mrigidBody =
		///	(ChBodySceneNode*) addChBodySceneNode_easyBox(&my_system, application.GetSceneManager(), 100.0, ChVector<>(0, -1.0, 0),
		///	ChQuaternion<>(1, 0, 0, 0), ChVector<>(10, (1.0-abs(GroundLoc)), 10));
		///mrigidBody->GetBody()->SetBodyFixed(true);
		///mrigidBody->GetBody()->GetMaterialSurface()->SetFriction(0.5);
		///mrigidBody->setMaterialTexture(0, rockMap);

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

        while (application.GetDevice()->run()) {
            application.GetVideoDriver()->beginScene(true, true);
            application.DrawAll();
            application.DoStep();
            application.GetVideoDriver()->endScene();
            if (application.GetPaused() == false) {
                std::cout << "Time t = " << my_system.GetChTime() << "s \n";
                AccuNoIterations += mystepper->GetNumIterations();
				printf("Forward position of rim X:      %12.4e ", Body_2->coord.pos.x);
            }
        }

        application.EndScene();

        if (my_system.GetStepcount() >= num_steps)
            application.GetDevice()->closeDevice();
    } else {
        double start = std::clock();
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
        std::cout << "Unit test check succeeded \n";
        system("pause");
    }

    return 0;
}
