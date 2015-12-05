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
// Authors: Antonio Recuero Alessandro Tasora
// =============================================================================
//
// Unit test for continuum-based bilinear shear deformable shell element using ANCF
//
// Successful execution of this unit test may validate: this element's elastic, isotropic
// force formulation and numerical integration implementations.
//
// The reference file data was validated by matching the steady-state response of the
// flat shell tip (i.e. only steady-state was validated) in the paper by Yamashita, Valkeapaa,
// Jayakumar, and Sugiyama, "Continuum Mechanics Based Bilinear Shear Deformable Shell
// Element Using Absolute Nodal Coordinate Formulation", ASME Journal of Computational and
// Nonlinear Dynamics, 10(5), 051012 (Sep 01, 2015). See its Figure 4.
//
// Gravity must be disabled; only a constant load of -50 N at a corner is used. Only
// 10 time steps are checked by default. User can increase this value up to 4000.
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
#include <cmath> 
#include "physics/ChLoadContainer.h" // For Alessandro classes
#include "chrono_mkl/ChLcpMklSolver.h"
#include "chrono_irrlicht/ChIrrApp.h"
#include "chrono_fea/ChVisualizationFEAmesh.h"
#include "chrono_fea/ChBuilderBeam.h"

// Remember to use the namespace 'chrono' because all classes
// of Chrono::Engine belong to this namespace and its children...

using namespace chrono;
using namespace fea;
using namespace irr;

const double precision = 1e-5;  // Used to accept/reject implementation


// Custom load n.1: a class that loads 4 ChNodeFEAxyzD 
// with a stiff load - for testing.

class MyLoadCustomMultiple : public ChLoadCustomMultiple {
public:   
        MyLoadCustomMultiple(std::vector< ChSharedPtr<ChLoadable> >& mloadables) :  ChLoadCustomMultiple(mloadables) {};

        // Compute Q=Q(x,v)
        // This is the function that you have to implement. It should return the generalized Q load 
        // (i.e.the force in generalized lagrangian coordinates).
        // Since here we have multiple connected ChLoadable objects (the two nodes), the rule is that
        // all the vectors (load_Q, state_x, state_w) are split in the same order that the loadable objects
        // are added to MyLoadCustomMultiple; in this case for instance Q={Efx,Efy,Efz,Ffx,Ffy,Ffz}.
        // As this is a stiff force field, dependency from state_x and state_y must be considered.
        virtual void ComputeQ(ChState*      state_x, ///< state position to evaluate Q
                                ChStateDelta* state_w  ///< state speed to evaluate Q
                        ) { 
            ChVector<> Node6_pos;
			ChVector<> Node6_vel;
			ChVector<> Node7_pos;
			ChVector<> Node7_vel;
			ChVector<> Node11_pos;
			ChVector<> Node11_vel;
			ChVector<> Node12_pos;
			ChVector<> Node12_vel;
			ChVector<> Node6_grad;
			ChVector<> Node6_gradvel;
			ChVector<> Node7_grad;
			ChVector<> Node7_gradvel;
			ChVector<> Node11_grad;
			ChVector<> Node11_gradvel;
			ChVector<> Node12_grad;
			ChVector<> Node12_gradvel;
            if (state_x && state_w) {
				Node6_pos  = state_x->ClipVector(0, 0);
				Node6_grad = state_x->ClipVector(3, 0);
				Node6_vel  = state_w->ClipVector(0, 0);
				Node6_gradvel = state_w->ClipVector(3, 0);

				Node7_pos  = state_x->ClipVector(6, 0);
				Node7_grad = state_x->ClipVector(9, 0);
				Node7_vel  = state_w->ClipVector(6, 0);
				Node7_gradvel = state_w->ClipVector(9, 0);

				Node11_pos  = state_x->ClipVector(12, 0);
				Node11_grad = state_x->ClipVector(15, 0);
				Node11_vel  = state_w->ClipVector(12, 0);
				Node11_gradvel = state_w->ClipVector(15, 0);

				Node12_pos  = state_x->ClipVector(18, 0);
				Node12_grad = state_x->ClipVector(21, 0);
				Node12_vel  = state_w->ClipVector(18, 0);
				Node12_gradvel = state_w->ClipVector(21, 0);
            }
            else { 
                // explicit integrators might call ComputeQ(0,0), null pointers mean
                // that we assume current state, without passing state_x for efficiency
				Node6_pos = loadables[0].DynamicCastTo<ChNodeFEAxyzD>()->GetPos();
				Node6_grad = loadables[0].DynamicCastTo<ChNodeFEAxyzD>()->GetD();
				Node6_vel = loadables[0].DynamicCastTo<ChNodeFEAxyzD>()->GetPos_dt();
				Node6_gradvel = loadables[0].DynamicCastTo<ChNodeFEAxyzD>()->GetD_dt();

				Node7_pos = loadables[1].DynamicCastTo<ChNodeFEAxyzD>()->GetPos();
				Node7_grad = loadables[1].DynamicCastTo<ChNodeFEAxyzD>()->GetD();
				Node7_vel = loadables[1].DynamicCastTo<ChNodeFEAxyzD>()->GetPos_dt();
				Node7_gradvel = loadables[1].DynamicCastTo<ChNodeFEAxyzD>()->GetD_dt();

				Node11_pos = loadables[2].DynamicCastTo<ChNodeFEAxyzD>()->GetPos();
				Node11_grad = loadables[2].DynamicCastTo<ChNodeFEAxyzD>()->GetD();
				Node11_vel = loadables[2].DynamicCastTo<ChNodeFEAxyzD>()->GetPos_dt();
				Node11_gradvel = loadables[2].DynamicCastTo<ChNodeFEAxyzD>()->GetD_dt();

				Node12_pos = loadables[3].DynamicCastTo<ChNodeFEAxyzD>()->GetPos();
				Node12_grad = loadables[3].DynamicCastTo<ChNodeFEAxyzD>()->GetD();
				Node12_vel = loadables[3].DynamicCastTo<ChNodeFEAxyzD>()->GetPos_dt();
				Node12_gradvel = loadables[3].DynamicCastTo<ChNodeFEAxyzD>()->GetD_dt();

            }

            // store generalized forces as a contiguous vector in this->load_Q, with same order of state_w
            this->load_Q.FillElem(0);
            /*
			this->load_Q(2) = -150 * pow(Node7_pos.x, 9) + 1500 * pow(Node7_vel.z*Node12_pos.z, 11)*abs(Node7_pos.z); // Fz 1
			this->load_Q(5) = 0.0 ; // RFz 1
			this->load_Q(8) = -150*pow(Node6_pos.x, 9) ; // Fz 2 
			this->load_Q(11) = 0.0; // RFz 2
			this->load_Q(14) = -150*pow(Node12_pos.x, 8) ; // Fz 3 
			this->load_Q(17) = 0.0; // RFz 3 
			this->load_Q(20) = -150 * pow(Node7_pos.x, 3) -  1500 * pow(Node6_pos.x, 2)* pow(Node7_vel.z, 3);// Fz 4 
            */
            this->load_Q(2) = -3000 * Node6_pos.z; // Fz 1
			this->load_Q(5) = 0.0 ; // RFz 1
			this->load_Q(8) = -3000 * Node7_pos.z; // Fz 2 
			this->load_Q(11) = 0.0; // RFz 2
			this->load_Q(14) = -10  * Node11_pos.z + 5; // Fz 3 
			this->load_Q(17) = 0.0; // RFz 3 
			this->load_Q(20) = -30000 * Node12_pos.z + 1;// Fz 4 
			this->load_Q(23) = 0.0;// RFz 4 
        }

        // OPTIONAL: if you want to provide an analytical jacobian, just implement the following:
        //   virtual void ComputeJacobian(...)

        // Remember to set this as stiff, to enable the jacobians
        virtual bool IsStiff() {return true;}
};



int main(int argc, char* argv[]) {

        

    // Utils to open/read files: Load reference solution ("golden") file
    ChMatrixDynamic<> FileInputMat(4000, 2);
    std::string shell_validation_file = GetChronoDataPath() + "testing/" + "UT_ANCFShellIso.txt";
    std::ifstream fileMid(shell_validation_file);
    if (!fileMid.is_open()) {
        fileMid.open(shell_validation_file);
    }
    if (!fileMid) {
        std::cout << "Cannot open file.\n";
        exit(1);
    }
    for (int x = 0; x < 4000; x++) {
        fileMid >> FileInputMat[x][0] >> FileInputMat[x][1];
    }
    fileMid.close();

    // Definition of the model
    ChSystem my_system;


    // Create the Irrlicht visualization (open the Irrlicht device,
    // bind a simple user interface, etc. etc.)
    ChIrrApp application(&my_system, L"ANCF & stiff loads", core::dimension2d<u32>(800, 600), false, true);

    // Easy shortcuts to add camera, lights, logo and sky in Irrlicht scene:
    application.AddTypicalLogo();
    application.AddTypicalSky();
    application.AddTypicalLights();
    application.AddTypicalCamera(core::vector3df(1.5f, 0.5f, 1.5f),  // camera location
                                 core::vector3df(0.0f, 0.0f, 0.0f));  // "look at" location


    // The physical system: it contains all physical objects.
    // Create a mesh, that is a container for groups
    // of elements and their referenced nodes.
    ChSharedPtr<ChMesh> my_mesh(new ChMesh);

    bool do_shells = true;
    bool do_cable = true;
    bool apply_custom_stiff_loads = true;

    ChSharedPtr<ChNodeFEAxyzD> node_log1;
    ChSharedPtr<ChNodeFEAxyzD> node_log2;

    const int num_steps = 200;        // Number of time steps for unit test (range 1 to 4000)
    const double time_step = 0.01;  // Time step


    // 
    // Object n.1:  
    // a nxn-nodes  ANCF shell:
    // 
    if (do_shells) {
        
        // Geometry of the plate
        double plate_lenght_x = 1;
        double plate_lenght_y = 1;
        double plate_lenght_z = 0.01;  // small thickness
        // Specification of the mesh
        const int numDiv_x = 4;
        const int numDiv_y = 4;
        const int numDiv_z = 1;
        const int N_x = numDiv_x + 1;
        const int N_y = numDiv_y + 1;
        const int N_z = numDiv_z + 1;
        // Number of elements in the z direction is considered as 1
        int TotalNumElements = numDiv_x * numDiv_y;
        int TotalNumNodes = (numDiv_x + 1) * (numDiv_y + 1);
        // For uniform mesh
        double dx = plate_lenght_x / numDiv_x;
        double dy = plate_lenght_y / numDiv_y;
        double dz = plate_lenght_z / numDiv_z;
        int MaxMNUM = 0;
        int MTYPE = 0;
        int MaxLayNum = 0;
        ChMatrixDynamic<double> COORDFlex(TotalNumNodes, 6);
        ChMatrixDynamic<double> VELCYFlex(TotalNumNodes, 6);
        ChMatrixDynamic<int> NumNodes(TotalNumElements, 4);
        ChMatrixDynamic<int> LayNum(TotalNumElements, 1);  // Only one layer in this unit test
        ChMatrixDynamic<int> NDR(TotalNumNodes, 6);
        ChMatrixDynamic<double> ElemLengthXY(TotalNumElements, 2);
        ChMatrixNM<double, 10, 12> MPROP;
        ChMatrixNM<int, 10, 7> MNUM;
        ChMatrixNM<int, 10, 1> NumLayer;
        double LayPROP[10][7][2];
        //!------------------------------------------------!
        //!--------------- Element data--------------------!
        //!------------------------------------------------!
        for (int i = 0; i < TotalNumElements; i++) {
            // All the elements belong to the same layer, e.g layer number 1.
            LayNum(i, 0) = 1;
            // Node number of the 4 nodes which creates element i.
            // The nodes are distributed this way. First in the x direction for constant
            // y when max x is reached go to the next level for y by doing the same
            // distribution but for y+1 and keep doing until y limit is reached. Node
            // number start from 1.
            NumNodes(i, 0) = (i / (numDiv_x)) * (N_x) + i % numDiv_x;
            NumNodes(i, 1) = (i / (numDiv_x)) * (N_x) + i % numDiv_x + 1;
            NumNodes(i, 2) = (i / (numDiv_x)) * (N_x) + i % numDiv_x + N_x;
            NumNodes(i, 3) = (i / (numDiv_x)) * (N_x) + i % numDiv_x + 1 + N_x;
            // Let's keep the element length a fixed number in both direction. (uniform
            // distribution of nodes in both direction)
            // If the user interested in non-uniform mesh it can be manipulated here.
            ElemLengthXY(i, 0) = dx;
            ElemLengthXY(i, 1) = dy;
            if (MaxLayNum < LayNum(i, 0)) {
                MaxLayNum = LayNum(i, 0);
            }
        }
        //!----------------------------------------------!
        //!--------- NDR,COORDFlex,VELCYFlex-------------!
        //!----------------------------------------------!

        for (int i = 0; i < TotalNumNodes; i++) {
            // If the node is the first node from the left side fix the x,y,z degree of
            // freedom+ d/dx, d/dy,d/dz as well. 1for constrained 0 for ...
            //-The NDR array is used to define the degree of freedoms that are
            // constrained in the 6 variable explained above.
            NDR(i, 0) = (i % (numDiv_x + 1) == 0) ? 1 : 0;
            NDR(i, 1) = (i % (numDiv_x + 1) == 0) ? 1 : 0;
            NDR(i, 2) = (i % (numDiv_x + 1) == 0) ? 1 : 0;
            NDR(i, 3) = (i % (numDiv_x + 1) == 0) ? 1 : 0;
            NDR(i, 4) = (i % (numDiv_x + 1) == 0) ? 1 : 0;
            NDR(i, 5) = (i % (numDiv_x + 1) == 0) ? 1 : 0;

            //-COORDFlex are the initial coordinates for each node,
            // the first three are the position and the last three are the slope
            // coordinates for the tangent vector.
            // Note that i starts from 0 but nodes starts from 1
            COORDFlex(i, 0) = (i % (numDiv_x + 1)) * dx;
            COORDFlex(i, 1) = (i / (numDiv_x + 1)) % (numDiv_y + 1) * dy;
            COORDFlex(i, 2) = (i) / ((numDiv_x + 1) * (numDiv_y + 1)) * dz;
            COORDFlex(i, 3) = 0;
            COORDFlex(i, 4) = 0;
            COORDFlex(i, 5) = 1;
            //-VELCYFlex is essentially the same as COORDFlex, but for the initial
            // velocity instead of position.
            // let's assume zero initial velocity for nodes
            VELCYFlex(i, 0) = 0;
            VELCYFlex(i, 1) = 0;
            VELCYFlex(i, 2) = 0;
            VELCYFlex(i, 3) = 0;
            VELCYFlex(i, 4) = 0;
            VELCYFlex(i, 5) = 0;
        }
        //!------------------------------------------------!
        //!------------- Read Layer Data-------------------!
        //!------------------------------------------------!

        for (int i = 0; i < MaxLayNum; i++) {
            NumLayer(i, 0) = i + 1;
            // For each layer define some properties
            // For multilayer problem the user should modify the following loop as they
            // would like
            for (int j = 0; j < NumLayer(i, 0); j++) {
                LayPROP[i][j][0] = dz;  // Layerheight
                LayPROP[i][j][1] = 0;   // PlyAngle
                MNUM[i][j] = 1;         // Material_ID
                if (MaxMNUM < MNUM(i, j))
                    MaxMNUM = MNUM(i, j);
            }
        }
        //!------------------------------------------------!
        //!------------ Input Material Data----------------!
        //!------------------------------------------------!
        for (int i = 0; i < MaxMNUM; i++) {
            double nu_coef = 0.3;
            MTYPE = 2;  // The user must use orthotropic input (MTYPE=2) to introduce isotropic material
                        // properties for this unit test

            if (MTYPE == 2) {
                MPROP(i, 0) = 500;      // Density [kg/m3]
                MPROP(i, 1) = 2.1E+08;  // Ex//
                MPROP(i, 2) = 2.1E+08;  // Ey
                MPROP(i, 3) = 2.1E+08;  // Ez
                // Additional information for the Type 2 of Material.
                MPROP(i, 4) = 0.3;                                    // nuxy
                MPROP(i, 5) = 0.3;                                    // nuxz
                MPROP(i, 6) = 0.3;                                    // nuyz
                MPROP(i, 7) = MPROP(i, 1) / 2.0 / (1 + MPROP(i, 6));  // Gxy
                MPROP(i, 8) = MPROP(i, 1) / 2.0 / (1 + MPROP(i, 6));  // Gxz
                MPROP(i, 9) = MPROP(i, 1) / 2.0 / (1 + MPROP(i, 6));  // Gyz
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
            if (NDR(i, 0) == 1 && NDR(i, 1) == 1 && NDR(i, 2) == 1 && NDR(i, 3) == 1 && NDR(i, 4) == 1 && NDR(i, 5) == 1) {
                node->SetFixed(true);
            }
            i++;
        }

        ChSharedPtr<ChNodeFEAxyzD> nodetip(my_mesh->GetNode(TotalNumNodes - 1).DynamicCastTo<ChNodeFEAxyzD>());
        // A random node just to try
        ChSharedPtr<ChNodeFEAxyzD> noderand(my_mesh->GetNode(TotalNumNodes / 2).DynamicCastTo<ChNodeFEAxyzD>());
        // A constrained node
        ChSharedPtr<ChNodeFEAxyzD> noderclamped(my_mesh->GetNode(0).DynamicCastTo<ChNodeFEAxyzD>());

        int elemcount = 0;
        while (elemcount < TotalNumElements) {
            ChSharedPtr<ChElementShellANCF> element(new ChElementShellANCF);
            // Save material data into InertFlexVec(98x1) at each layer
            ChMatrixNM<double, 98, 1> InertFlexVec;
            InertFlexVec.Reset();
            double TotalThickness;  // element thickness
            TotalThickness = 0.0;
            int i = elemcount;
            for (int j = 0; j < NumLayer(LayNum(i, 0) - 1, 0); j++) {
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
            ChMatrixNM<double, 7, 2> GaussZRange;
            GaussZRange.Reset();
            double CurrentHeight = 0.0;
            for (int j = 0; j < NumLayer(LayNum(i, 0) - 1, 0); j++) {
                double AA = (CurrentHeight / TotalThickness - 0.5) * 2.0;
                CurrentHeight += LayPROP[LayNum(i, 0) - 1][j][0];
                double AAA = (CurrentHeight / TotalThickness - 0.5) * 2.0;
                GaussZRange(j, 0) = AA;
                GaussZRange(j, 1) = AAA;
            }
            element->SetInertFlexVec(InertFlexVec);
            element->SetGaussZRange(GaussZRange);
            element->SetNodes(my_mesh->GetNode(NumNodes[elemcount][0]).DynamicCastTo<ChNodeFEAxyzD>(),
                              my_mesh->GetNode(NumNodes[elemcount][1]).DynamicCastTo<ChNodeFEAxyzD>(),
                              my_mesh->GetNode(NumNodes[elemcount][2]).DynamicCastTo<ChNodeFEAxyzD>(),
                              my_mesh->GetNode(NumNodes[elemcount][3]).DynamicCastTo<ChNodeFEAxyzD>());
            element->SetMaterial(mmaterial);
            element->SetNumLayers(NumLayer(LayNum(i, 0) - 1, 0));
            element->SetThickness(TotalThickness);
            element->SetElemNum(elemcount);
            element->SetAlphaDamp(0.08);
            element->Setdt(0.001);                     // dt to calculate DampingCoefficient
            element->SetGravityOn(false);              // turn gravity on/off
            ChMatrixNM<double, 35, 1> StockAlpha_EAS;  // StockAlpha(5*7,1): Max #Layer is 7
            StockAlpha_EAS.Reset();
            element->SetStockAlpha(StockAlpha_EAS);
            my_mesh->AddElement(element);
            elemcount++;
        }

        // APPLY SOME CONSTANT LOAD (the easy way, as a force in a single node)
        nodetip->SetForce(ChVector<>(0,0,-10));

	    // APPLY SOME LOADS! (using ChLoad classes)

	    ChSharedPtr<ChNodeFEAxyzD> NodeLoad6(my_mesh->GetNode(6).DynamicCastTo<ChNodeFEAxyzD>());
	    // A random node just to try
	    ChSharedPtr<ChNodeFEAxyzD> NodeLoad7(my_mesh->GetNode(7).DynamicCastTo<ChNodeFEAxyzD>());
	    // A constrained node
	    ChSharedPtr<ChNodeFEAxyzD> NodeLoad11(my_mesh->GetNode(11).DynamicCastTo<ChNodeFEAxyzD>());
	    // A constrained node
	    ChSharedPtr<ChNodeFEAxyzD> NodeLoad12(my_mesh->GetNode(12).DynamicCastTo<ChNodeFEAxyzD>());

	    // First: loads must be added to "load containers", 
	    // and load containers must be added to your ChSystem 
	    ChSharedPtr< ChLoadContainer > mloadcontainer(new ChLoadContainer);
 	    my_system.Add(mloadcontainer);

  
       // Instance load object. This require a list of ChLoadable objects
       // (these are our two nodes,pay attention to the sequence order), and add to container.
       std::vector< ChSharedPtr< ChLoadable > > mnodelist;
       mnodelist.push_back(NodeLoad6);
       mnodelist.push_back(NodeLoad7);
       mnodelist.push_back(NodeLoad11);
       mnodelist.push_back(NodeLoad12);
       ChSharedPtr< MyLoadCustomMultiple > mloadcustommultiple(new MyLoadCustomMultiple(mnodelist));
       if (apply_custom_stiff_loads)
            mloadcontainer->Add(mloadcustommultiple);

       node_log1 = nodetip;
       node_log2 = NodeLoad11;
    }


    // 
    // Object n.2:  
    // a 4-nodes  ANCF cable:
    //

    if (do_cable) {

	    ChSharedPtr<ChBeamSectionCable> msection_cable2(new ChBeamSectionCable);
	    msection_cable2->SetDiameter(0.05);
	    msection_cable2->SetYoungModulus (0.01e9);
	    msection_cable2->SetBeamRaleyghDamping(0.05);

	    ChBuilderBeamANCF builder;

	    builder.BuildBeam(	my_mesh,		// the mesh where to put the created nodes and elements 
						    msection_cable2,// the ChBeamSectionCable to use for the ChElementBeamANCF elements
						    4,				// the number of ChElementBeamANCF to create
						    ChVector<>(0, 1.2, 0),	// the 'A' point in space (beginning of beam)
						    ChVector<>(1, 1.2, 0));	// the 'B' point in space (end of beam)

        builder.GetLastBeamNodes().front()->SetFixed(true);

        // APPLY SOME CONSTANT LOAD! (the easy way, as a force in a single node)
	    builder.GetLastBeamNodes().back()->SetForce( ChVector<>(0,0,-1));

        // APPLY SOME STIFF LOADS!  (using ChLoad classes)

	    // First: loads must be added to "load containers", 
	    // and load containers must be added to your ChSystem 
	    ChSharedPtr< ChLoadContainer > mloadcontainer(new ChLoadContainer);
 	    my_system.Add(mloadcontainer);

  
        // Instance load object. This require a list of ChLoadable objects
        // (these are our two nodes,pay attention to the sequence order), and add to container.
        std::vector< ChSharedPtr< ChLoadable > > mnodelist;
        mnodelist.push_back(builder.GetLastBeamNodes()[0]);
        mnodelist.push_back(builder.GetLastBeamNodes()[1]);
        mnodelist.push_back(builder.GetLastBeamNodes()[2]);
        mnodelist.push_back(builder.GetLastBeamNodes()[3]);
        ChSharedPtr< MyLoadCustomMultiple > mloadcustommultiple_for_beam(new MyLoadCustomMultiple(mnodelist));
        if (apply_custom_stiff_loads)
            mloadcontainer->Add(mloadcustommultiple_for_beam);

    }


    // Switch off mesh class gravity
    my_mesh->SetAutomaticGravity(false);
    // This is mandatory
    my_mesh->SetupInitial();
    // Remember to add the mesh to the system!
    my_system.Add(my_mesh);
    
    /*
    my_system.SetLcpSolverType(ChSystem::LCP_ITERATIVE_MINRES);  // <- NEEDED because other solvers can't
                                                                 // handle stiffness matrices
    chrono::ChLcpIterativeMINRES* msolver = (chrono::ChLcpIterativeMINRES*)my_system.GetLcpSolverSpeed();
    msolver->SetDiagonalPreconditioning(true);
    my_system.SetIterLCPwarmStarting(true);  // this helps a lot to speedup convergence in this class of problems
    my_system.SetIterLCPmaxItersSpeed(100);
    my_system.SetIterLCPmaxItersStab(100);
    my_system.SetTolForce(1e-10);
    */

    ChLcpMklSolver * mkl_solver_stab = new ChLcpMklSolver; // MKL Solver option
    ChLcpMklSolver * mkl_solver_speed = new ChLcpMklSolver;
    my_system.ChangeLcpSolverStab(mkl_solver_stab);
    my_system.ChangeLcpSolverSpeed(mkl_solver_speed);
    mkl_solver_stab->SetSparsityPatternLock(true);
	mkl_solver_speed->SetSparsityPatternLock(true);

    my_system.Setup();
    my_system.Update();

    // INT_HHT or INT_EULER_IMPLICIT
    my_system.SetIntegrationType(ChSystem::INT_HHT);

    ChSharedPtr<ChTimestepperHHT> mystepper = my_system.GetTimestepper().DynamicCastTo<ChTimestepperHHT>();
    //mystepper->SetAlpha(0.0);
    mystepper->SetMaxiters(20);
    mystepper->SetTolerance(1e-06);
    mystepper->SetMode(ChTimestepperHHT::POSITION);
    mystepper->SetScaling(true);

    utils::Data m_data;
    m_data.resize(2);
    for (size_t col = 0; col < 2; col++)
        m_data[col].resize(num_steps);
    utils::CSV_writer csv(" ");
    std::ifstream file2("UT_ANCFShellWith.txt");



    // Options for visualization in irrlicht

    ChSharedPtr<ChVisualizationFEAmesh> mvisualizemesh(new ChVisualizationFEAmesh(*(my_mesh.get_ptr())));
    mvisualizemesh->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_SURFACE);
    mvisualizemesh->SetSmoothFaces(false);
    my_mesh->AddAsset(mvisualizemesh);

    ChSharedPtr<ChVisualizationFEAmesh> mvisualizemeshwire(new ChVisualizationFEAmesh(*(my_mesh.get_ptr())));
    mvisualizemeshwire->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_SURFACE);
    mvisualizemeshwire->SetWireframe(true);
    my_mesh->AddAsset(mvisualizemeshwire);

    ChSharedPtr<ChVisualizationFEAmesh> mvisualizemeshC(new ChVisualizationFEAmesh(*(my_mesh.get_ptr())));
    mvisualizemeshC->SetFEMglyphType(ChVisualizationFEAmesh::E_GLYPH_NODE_DOT_POS);
    mvisualizemeshC->SetFEMdataType(ChVisualizationFEAmesh::E_PLOT_NONE);
    mvisualizemeshC->SetSymbolsThickness(0.01);
    my_mesh->AddAsset(mvisualizemeshC);


    application.AssetBindAll();
    application.AssetUpdateAll();

    my_system.Setup();
    my_system.Update();

    GetLog() << "\n\nREADME\n\n" << " - Press SPACE to start dynamic simulation \n - Press F10 for nonlinear statics - Press F11 for linear statics. \n";
    
    // at beginning, no analysis is running..
    application.SetPaused(true);

    int AccuNoIterations = 0;

    while (application.GetDevice()->run()) {
        application.BeginScene();
        application.DrawAll();
        application.DoStep();
        
        if (application.GetPaused() == false) {
            std::cout << "Time t = " << my_system.GetChTime() << "s \n";
            if (node_log1 && &node_log2) {
                std::cout << "node1->pos.z = " << node_log1->pos.z << "\n";
                std::cout << "node2->pos.z = " << node_log2->pos.z << "\n";
                AccuNoIterations += mystepper->GetNumIterations();
	            m_data[0][my_system.GetStepcount()] = my_system.GetChTime();
	            m_data[1][my_system.GetStepcount()] = node_log2->pos.z;
	            csv << m_data[0][my_system.GetStepcount()] << m_data[1][my_system.GetStepcount()] << std::endl;
	            csv.write_to_file("UT_ANCFShellWith.txt");
            }
        }

        application.EndScene();

        if(my_system.GetStepcount() >= num_steps)
            application.GetDevice()->closeDevice();
    }


    // Code snippet to generate golden file

    return 0;
}
