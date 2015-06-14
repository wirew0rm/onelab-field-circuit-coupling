Include "transformer_data.geo";
Include "BH.pro"; // nonlinear BH caracteristic of magnetic material

Dir="res/";
ExtGmsh     = ".pos";
ExtGnuplot  = ".dat";
po = "Output/";

Function {
	// Geometry Parameters
  DefineConstant[
		Flag_InitFromPrevious = {0, Choices{0,1}, Name "Input/42InitSolutionFromPrevious"},
		Flag_BC_Type = {1, Choices{0="Neumann",1="Dirichlet"}, ReadOnly (Flag_Infinity==1),
			Name "Input/20Boundary condition at infinity",
			Highlight "Blue"},
		Lz = 1,
		SymmetryFactor = 1
  ];
}

// Newton Raphson Parameters
Function {
  DefineConstant[
		Flag_NL = { 1, Choices{0,1},
			Name "Input/40Nonlinear BH-curve"}
		Nb_max_iter = 30,
		relaxation_factor = 1,
		stop_criterion = 1e-5,
		reltol = 1e-7,
		abstol = 1e-5,
		Flag_NL_Newton_Raphson = {1, Choices{0,1}, Name "Input/41Newton-Raphson iteration",
			Visible Flag_NL}
  ];
}

Group {
  Core = Region[ {ECORE, ICORE} ];
  Ind_1      = Region[{COIL}] ;
  Ind_2      = Region[{(COIL+1)}] ;
	Ind_1_ = Region[{(COIL+2)}] ;
	Ind_2_ = Region[{(COIL+3)}] ;
	Inds1   = Region[ {Ind_1, Ind_1_} ];
	Inds2   = Region[ {Ind_2, Ind_2_} ];

  AirGap = Region[ AIRGAP ];

  Air  = Region[ AIR ];
  AirInf = Region[ AIRINF ];

  If(Flag_OpenCore)
    Air  += Region[ {AirGap} ];
  EndIf
  If(!Flag_OpenCore)
    Core += Region[ {AirGap} ];
  EndIf

	// Surfaces for imposing boundary conditions
	DefineGroup[Surf_bn0];
	If(Flag_BC_Type==1)
		Surf_Inf = Region[ SURF_AIROUT ];
	EndIf
	If(Flag_Symmetry)
		Surf_bn0 = Region[ {AXIS_Y} ];
	EndIf

	DomainB   = Region[ {Inds1, Inds2} ];
	DomainB1  = Region[ {Inds1} ];
	DomainB2  = Region[ {Inds2} ];

	If(Flag_Infinity)
		DomainInf = Region[ {AirInf} ];
	EndIf

	DomainCC = Region[ {Air, AirInf, Inds1, Inds2, Core} ];
	DomainC  = Region[ { } ];
	Domain  = Region[ {DomainCC, DomainC} ];

	If(Flag_NL)
		DomainNL = Region[ {Core} ];
		DomainL  = Region[ {Domain,-DomainNL} ];
	EndIf
	DomainDummy = #12345 ; // Dummy region number for postpro with functions
}

// Source Current Distribution
Function {
  DefineConstant[
    II1 = { IA1,
      Name "Input/4Coil Parameters/0Current1 (rms) [A]", Highlight "AliceBlue" },
    NbWires1 = { Nw1,
      Name "Input/4Coil Parameters/1Number of turns1", Highlight "AliceBlue"},
    II2 = { IA2,
      Name "Input/4Coil Parameters/0Current2 (rms) [A]", Highlight "AliceBlue" },
    NbWires2 = { Nw2,
      Name "Input/4Coil Parameters/1Number of turns2", Highlight "AliceBlue"}
  ];
  NbWires1[]  = NbWires1 ;
  NbWires2[]  = NbWires2 ;
	SurfCoil1[] = SurfaceArea[]{COIL} ;
	SurfCoil2[] = SurfaceArea[]{(COIL+1)} ;
	Idir1[#{COIL}]     =  1. ;
	Idir1[#{(COIL+2)}] = -1. ;
	vDir1[]   = Vector[ 0, 0, Idir1[]] ;
	Idir2[#{(COIL+1)}] =  1. ;
	Idir2[#{(COIL+3)}] = -1. ;
	vDir2[]   = Vector[ 0, 0, Idir2[]] ;
}

// Material properties
Function {
  DefineConstant[
    mur_fe = { 2000., Min 100, Max 2000, Step 100,
      Name "Input/42Core relative permeability", Highlight "AliceBlue",
      Visible (!Flag_NL)}
  ];
  mu0 = 4.e-7 * Pi ;
  nu [#{Air, AirInf, Inds1, Inds2}]  = TensorDiag[1.0, 1.0, 1.0] * 1./mu0 ;
  If(!Flag_NL)
    nu [#{Core}]  = 1/(mur_fe*mu0) ;
  EndIf
  If(Flag_NL)
    nu [ DomainNL ] = TensorDiag[1.0, 1.0, 1.0] * nu_EIcore[$1] ;
    dhdb_NL [ DomainNL ] = dhdb_EIcore_NL[$1];
    dhdb [ DomainNL ] = dhdb_EIcore[$1];
  EndIf
}


//-------------------------------------------------------------------------------------

Jacobian {
  { Name Vol;
    Case {
      { Region DomainInf ; Jacobian VolSphShell{Val_Rint, Val_Rext} ; }
      { Region All ; Jacobian Vol; }
    }
  }
}

Integration {
  { Name I1 ; Case {
      { Type Gauss ;
        Case {
          { GeoElement Triangle   ; NumberOfPoints  6 ; }
				  { GeoElement Quadrangle ; NumberOfPoints  4 ; }
	  			{ GeoElement Line       ; NumberOfPoints  13 ; }
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------

Constraint {
  { Name MVP_2D ;
    Case {
			// outer boundary condition
      { Region Surf_Inf ; Type Assign ; Value 0. ; }
			// Symmetry Boundary condition
      { Region Surf_bn0 ; Type Assign ; Value 0. ; }
    }
  }

	// Initialize a with the value from the previous run.
	// getdp has to be called with -gmshread res/a.pos flag
	{ Name InitA ;
		Case {
			{ Region Domain ; Type Init ; Value Field[XYZ[]] ; }
		}
	}
}

//-----------------------------------------------------------------------------------------------

FunctionSpace {
  // Magnetic Vector Potential
  { Name Hcurl_a_2D ; Type Form1P ;
    BasisFunction {
      { Name se1 ; NameOfCoef ae1 ; Function BF_PerpendicularEdge ;
        Support Region[{Domain}] ; Entity NodesOf [ All ] ; }
   }
    Constraint {
			// Boundary conditions and symmetry
      { NameOfCoef ae1 ; EntityType NodesOf ; NameOfConstraint MVP_2D ; }
			// Apply initial Values to vector Potential
			If (Flag_InitFromPrevious)
				{ NameOfCoef ae1 ; EntityType NodesOf ; NameOfConstraint InitA ; }
			EndIf
    }
  }
	// Magnetic Vector Potential for postprocessing
  { Name Hcurl_a_2D_L ; Type Form1P ;
    BasisFunction {
      { Name se1 ; NameOfCoef ae1 ; Function BF_PerpendicularEdge ;
        Support Region[{Domain}] ; Entity NodesOf [ All ] ; }
   }
    Constraint {
			// Boundary conditions and symmetry
      { NameOfCoef ae1 ; EntityType NodesOf ; NameOfConstraint MVP_2D ; }
    }
  }
}

//-----------------------------------------------------------------------------------------------

Formulation {
  { Name MagStaDyn_a_2D ; Type FemEquation ;
    Quantity {
      { Name a  ; Type Local  ; NameOfSpace Hcurl_a_2D ; }
    }
    Equation {
      Galerkin { [ nu[{d a}] * Dof{d a}  , {d a} ] ;
        In Domain ; Jacobian Vol ; Integration I1 ; }

			If(Flag_NL)
      If(Flag_NL_Newton_Raphson)
      Galerkin { JacNL [ dhdb_NL[{d a}] * Dof{d a} , {d a} ] ;
        In DomainNL ; Jacobian Vol ; Integration I1 ; }
      EndIf
      EndIf

			Galerkin { [ -NbWires1[]/SurfCoil1[] * vDir1[] * II1, {a} ] ;
			  In DomainB1 ; Jacobian Vol ; Integration I1 ; }
			Galerkin { [ -NbWires2[]/SurfCoil2[] * vDir2[] * II2, {a} ] ;
			  In DomainB2 ; Jacobian Vol ; Integration I1 ; }
    }
  }

  { Name MagStaDyn_a_2D_L1; Type FemEquation ;
    Quantity {
      { Name a  ; Type Local  ; NameOfSpace Hcurl_a_2D_L ; }
    }

    Equation {
      Galerkin { [ Field[XYZ[]]{55} * Dof{d a}  , {d a} ] ;
        In Domain ; Jacobian Vol ; Integration I1 ; }

      Galerkin { [ -NbWires1[]/SurfCoil1[] * vDir1[] , {a} ] ;
        In DomainB1 ; Jacobian Vol ; Integration I1 ; }
    }
	}
  { Name MagStaDyn_a_2D_L2; Type FemEquation ;
    Quantity {
      { Name a  ; Type Local  ; NameOfSpace Hcurl_a_2D_L ; }
    }

    Equation {
      Galerkin { [ Field[XYZ[]]{55} * Dof{d a}  , {d a} ] ;
        In Domain ; Jacobian Vol ; Integration I1 ; }

      Galerkin { [ -NbWires2[]/SurfCoil2[] * vDir2[] , {a} ] ;
        In DomainB2 ; Jacobian Vol ; Integration I1 ; }
    }
	}

  { Name MagStaDyn_a_2D_Ld1; Type FemEquation ;
    Quantity {
      { Name a  ; Type Local  ; NameOfSpace Hcurl_a_2D_L ; }
    }

    Equation {
      Galerkin { [ Field[XYZ[]]{66} * Dof{d a}  , {d a} ] ;
        In Domain ; Jacobian Vol ; Integration I1 ; }

      Galerkin { [ -NbWires1[]/SurfCoil1[] * vDir1[] , {a} ] ;
        In DomainB1 ; Jacobian Vol ; Integration I1 ; }
    }
	}
  { Name MagStaDyn_a_2D_Ld2; Type FemEquation ;
    Quantity {
      { Name a  ; Type Local  ; NameOfSpace Hcurl_a_2D_L ; }
    }

    Equation {
      Galerkin { [ Field[XYZ[]]{66} * Dof{d a}  , {d a} ] ;
        In Domain ; Jacobian Vol ; Integration I1 ; }

      Galerkin { [ -NbWires2[]/SurfCoil2[] * vDir2[] , {a} ] ;
        In DomainB2 ; Jacobian Vol ; Integration I1 ; }
    }
	}
}

//-----------------------------------------------------------------------------------------------

Resolution {
  { Name Analysis ;
    System {
		{ Name A ; NameOfFormulation MagStaDyn_a_2D ; }
		{ Name A_L1 ; NameOfFormulation MagStaDyn_a_2D_L1 ; }
		{ Name A_L2 ; NameOfFormulation MagStaDyn_a_2D_L2 ; }
		{ Name A_Ld1 ; NameOfFormulation MagStaDyn_a_2D_Ld1 ; }
		{ Name A_Ld2 ; NameOfFormulation MagStaDyn_a_2D_Ld2 ; }
    }
    Operation {
      CreateDir["res/"];
      InitSolution[A] ;
			If(!Flag_NL)
				Generate[A] ; Solve[A] ;
			EndIf
			If(Flag_NL)
				IterativeLoop[Nb_max_iter, stop_criterion, relaxation_factor]{
					GenerateJac[A] ; SolveJac[A] ; }
			EndIf
			SaveSolution[A] ;
			PostOperation[Get_LocalFields] ;
			PostOperation[Get_GlobalQuantities] ;

			// get Inductances
      InitSolution[A_L1] ;
      Generate[A_L1] ; Solve[A_L1]; SaveSolution[A_L1] ;
      PostOperation[Get_GlobalQuantities_L1] ;  
      InitSolution[A_L2] ;
      Generate[A_L2] ; Solve[A_L2]; SaveSolution[A_L2] ;
      PostOperation[Get_GlobalQuantities_L2] ;  

			// get differential Inductances
      InitSolution[A_Ld1] ;
      Generate[A_Ld1] ; Solve[A_Ld1]; SaveSolution[A_Ld1] ;
      PostOperation[Get_GlobalQuantities_Ld1] ;  
      InitSolution[A_Ld2] ;
      Generate[A_Ld2] ; Solve[A_Ld2]; SaveSolution[A_Ld2] ;
      PostOperation[Get_GlobalQuantities_Ld2] ;  
    }
  }
}

//-----------------------------------------------------------------------------------------------

PostProcessing {

  { Name MagStaDyn_a_2D ; NameOfFormulation MagStaDyn_a_2D ;
    PostQuantity {
      { Name az ; Value { Term { [ CompZ[{a}] ] ; In Domain ; Jacobian Vol ; } } }

      { Name dhdb ; Value {
				If(!Flag_NL)
          Term { [ nu[{d a}]]; In Domain; Jacobian Vol; }
				EndIf
				If(Flag_NL)
          Term { [ nu[{d a}]]; In DomainL; Jacobian Vol; }
          Term { [ dhdb[{d a}] ] ; In DomainNL ; Jacobian Vol ; }
				EndIf
        }
      }
      { Name nu ; Value {
          Term { [ nu[{d a}]]; In Domain; Jacobian Vol; }
        }
      }

      { Name MagEnergy ; Value {
          Integral { [ SymmetryFactor*Lz* 1/2 *nu[{d a}]*{d a}*{d a} ] ;
            In Domain ; Jacobian Vol ; Integration I1 ; } 
				} 
			}

      { Name Flux1 ; Value {
          Integral { [ SymmetryFactor*Lz*NbWires1[]/SurfCoil1[] * Idir1[]* CompZ[{a}] ] ;
            In Inds1  ; Jacobian Vol ; Integration I1 ; }
        }
      }

      { Name Flux2 ; Value {
          Integral { [ SymmetryFactor*Lz*NbWires2[]/SurfCoil2[] * Idir2[]* CompZ[{a}] ] ;
            In Inds2  ; Jacobian Vol ; Integration I1 ; }
        }
      }
    }//PostQuantity
  }// PostProcessing:MagStaDyn_a_2D

  { Name MagStaDyn_a_2D_L1 ; NameOfFormulation MagStaDyn_a_2D_L1 ;
    PostQuantity {
      { Name L1; Value {
          Integral { [ SymmetryFactor*Lz*NbWires1[]/SurfCoil1[] * Idir1[]* CompZ[{a}] ] ;
            In Inds1  ; Jacobian Vol ; Integration I1 ; }
        }
      }
      { Name LC; Value {
          Integral { [ SymmetryFactor*Lz*NbWires2[]/SurfCoil2[] * Idir2[]* CompZ[{a}] ] ;
            In Inds2  ; Jacobian Vol ; Integration I1 ; }
        }
      }
    }
  }
  { Name MagStaDyn_a_2D_L2 ; NameOfFormulation MagStaDyn_a_2D_L2 ;
    PostQuantity {
      { Name L2; Value {
          Integral { [ SymmetryFactor*Lz*NbWires2[]/SurfCoil2[] * Idir2[]* CompZ[{a}] ] ;
            In Inds2  ; Jacobian Vol ; Integration I1 ; }
        }
      }
    }
  }

  { Name MagStaDyn_a_2D_Ld1 ; NameOfFormulation MagStaDyn_a_2D_Ld1 ;
    PostQuantity {
      { Name Ld1; Value {
          Integral { [ SymmetryFactor*Lz*NbWires1[]/SurfCoil1[] * Idir1[]* CompZ[{a}] ] ;
            In Inds1  ; Jacobian Vol ; Integration I1 ; }
        }
      }
      { Name LdC; Value {
          Integral { [ SymmetryFactor*Lz*NbWires2[]/SurfCoil2[] * Idir2[]* CompZ[{a}] ] ;
            In Inds2  ; Jacobian Vol ; Integration I1 ; }
        }
      }
    }
  }
  { Name MagStaDyn_a_2D_Ld2 ; NameOfFormulation MagStaDyn_a_2D_Ld2 ;
    PostQuantity {
      { Name Ld2; Value {
          Integral { [ SymmetryFactor*Lz*NbWires2[]/SurfCoil2[] * Idir2[]* CompZ[{a}] ] ;
            In Inds2  ; Jacobian Vol ; Integration I1 ; }
        }
      }
    }
  }
}// PostProcessing

//-----------------------------------------------------------------------------------------------

PostOperation Get_LocalFields UsingPost MagStaDyn_a_2D {
  Print[ az, OnElementsOf Domain, File StrCat[Dir,StrCat["a",ExtGmsh]], LastTimeStepOnly ];
  Print[ nu, OnElementsOf Domain, LastTimeStepOnly, StoreInField 55 ];
  Print[ dhdb, OnElementsOf Domain, LastTimeStepOnly, StoreInField 66 ];
}

PostOperation Get_GlobalQuantities UsingPost MagStaDyn_a_2D {
  Print[ Flux1[Inds1], OnGlobal, Format TimeTable,
    File > StrCat[Dir,StrCat["Flux1",ExtGnuplot]],Store 11,
    SendToServer StrCat[po,"40Flux [Wb]1"],  Color "LightYellow" ];
  Print[ Flux2[Inds2], OnGlobal, Format TimeTable,
    File > StrCat[Dir,StrCat["Flux2",ExtGnuplot]],Store 12,
    SendToServer StrCat[po,"40Flux [Wb]2"],  Color "LightYellow" ];

  Print[ MagEnergy[Domain], OnGlobal, Format TimeTable,
    File > StrCat[Dir,StrCat["ME",ExtGnuplot]], LastTimeStepOnly, Store 22,
    SendToServer StrCat[po,"99Magnetic Energy [W]"],  Color "LightYellow" ];
}

PostOperation Get_GlobalQuantities_L1 UsingPost MagStaDyn_a_2D_L1 {
  Print[ L1[Inds1], OnGlobal, Format TimeTable,
         File StrCat[Dir,StrCat["Inductance_L1",ExtGnuplot]],
         SendToServer StrCat[po,"54Inductance from Flux[H]1"], Color "LightYellow" ];
  Print[ LC[Inds2], OnGlobal, Format TimeTable,
         File StrCat[Dir,StrCat["Inductance_LC",ExtGnuplot]],
         SendToServer StrCat[po,"54Inductance from Flux[H]C"], Color "LightYellow" ];
}  
PostOperation Get_GlobalQuantities_L2 UsingPost MagStaDyn_a_2D_L2 {
  Print[ L2[Inds2], OnGlobal, Format TimeTable,
         File StrCat[Dir,StrCat["Inductance_L2",ExtGnuplot]],
         SendToServer StrCat[po,"54Inductance from Flux[H]2"], Color "LightYellow" ];
}  
PostOperation Get_GlobalQuantities_Ld1 UsingPost MagStaDyn_a_2D_Ld1 {
  Print[ Ld1[Inds1], OnGlobal, Format TimeTable,
         File StrCat[Dir,StrCat["Inductance_Ld1",ExtGnuplot]],
         SendToServer StrCat[po,"55Inductance from Flux Diff[H]1"], Color "LightYellow" ];
  Print[ LdC[Inds2], OnGlobal, Format TimeTable,
         File StrCat[Dir,StrCat["Inductance_LdC",ExtGnuplot]],
         SendToServer StrCat[po,"55Inductance from Flux Diff[H]C"], Color "LightYellow" ];
}  
PostOperation Get_GlobalQuantities_Ld2 UsingPost MagStaDyn_a_2D_Ld2 {
  Print[ Ld2[Inds2], OnGlobal, Format TimeTable,
         File StrCat[Dir,StrCat["Inductance_Ld2",ExtGnuplot]],
         SendToServer StrCat[po,"55Inductance from Flux Diff[H]2"], Color "LightYellow" ];
}  
 
DefineConstant[
  R_ = {"Analysis", Name "GetDP/1ResolutionChoices", Visible 0},
  C_ = {"-solve -v2 -gmshread res/a.pos", Name "GetDP/9ComputeCommand", Visible 0},
  P_ = {"", Name "GetDP/2PostOperationChoices", Visible 0}
];
