Include "inductor_data.geo";
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

	Ind_1_ = Region[{(COIL+1)}] ;
	Inds   = Region[ {Ind_1, Ind_1_} ];

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

	DomainB   = Region[ {Inds} ];

	If(Flag_Infinity)
		DomainInf = Region[ {AirInf} ];
	EndIf

	DomainCC = Region[ {Air, AirInf, Inds, Core} ];
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
    II = { IA,
//      Name "Input/4Coil Parameters/0Current (rms) [A]", Highlight "AliceBlue" , Graph "2000" },
      Name "Input/4Coil Parameters/0Current (rms) [A]", Highlight "AliceBlue" },
    NbWires = { Nw,
      Name "Input/4Coil Parameters/1Number of turns", Highlight "AliceBlue"}
  ];
  NbWires[]  = NbWires ;
	SurfCoil[] = SurfaceArea[]{COIL} ;
	Idir[#{COIL}]     =  1. ;
	Idir[#{(COIL+1)}] = -1. ;
	vDir[]   = Vector[ 0, 0, Idir[]] ;
  js0[] = NbWires[]/SurfCoil[] * vDir[] ;
}

// Material properties
Function {
  DefineConstant[
    mur_fe = { 2000., Min 100, Max 2000, Step 100,
      Name "Input/42Core relative permeability", Highlight "AliceBlue",
      Visible (!Flag_NL)}
  ];
  mu0 = 4.e-7 * Pi ;
  nu [#{Air, AirInf, Inds}]  = TensorDiag[1.0, 1.0, 1.0] * 1./mu0 ;
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
  { Name Hcurl_a_2D_Ld ; Type Form1P ;
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

			Galerkin { [ -NbWires[]/SurfCoil[] * vDir[] * II , {a} ] ;
			  In DomainB ; Jacobian Vol ; Integration I1 ; }
    }
  }

  { Name MagStaDyn_a_2D_Ld; Type FemEquation ;
    Quantity {
      { Name a  ; Type Local  ; NameOfSpace Hcurl_a_2D_Ld ; }
    }

    Equation {
      Galerkin { [ Field[XYZ[]]{66} * Dof{d a}  , {d a} ] ;
        In Domain ; Jacobian Vol ; Integration I1 ; }

      Galerkin { [ -NbWires[]/SurfCoil[] * vDir[] , {a} ] ;
        In DomainB ; Jacobian Vol ; Integration I1 ; }
    }
  }
}

//-----------------------------------------------------------------------------------------------

Resolution {
  { Name Analysis ;
    System {
		{ Name A ; NameOfFormulation MagStaDyn_a_2D ; }
		{ Name A_Ld ; NameOfFormulation MagStaDyn_a_2D_Ld ; }
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

			// get differential Inductance
      InitSolution[A_Ld] ;
      Generate[A_Ld] ; Solve[A_Ld]; SaveSolution[A_Ld] ;
      PostOperation[Get_GlobalQuantities_Ld] ;  
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

      { Name MagEnergy ; Value {
          Integral { [ SymmetryFactor*Lz* 1/2 *nu[{d a}]*{d a}*{d a} ] ;
            In Domain ; Jacobian Vol ; Integration I1 ; } } }

      { Name Flux ; Value {
          Integral { [ SymmetryFactor*Lz*Idir[]*NbWires[]/SurfCoil[]* CompZ[{a}] ] ;
            In Inds  ; Jacobian Vol ; Integration I1 ; }
        }
      }

      { Name Inductance_from_Flux ; Value { 
				Term { Type Global; [ #11/II ] ; // Flux stored in register #11
					In DomainDummy ; }
				} 
			}

      { Name Inductance_from_MagEnergy ; Value {
				Term { Type Global; [ 2*#22/(II*II) ] ; // MagEnergy stored in register #22
					In DomainDummy ; }
				}
			}
    }//PostQuantity
  }// PostProcessing:MagStaDyn_a_2D

  { Name MagStaDyn_a_2D_Ld ; NameOfFormulation MagStaDyn_a_2D_Ld ;
    PostQuantity {
      { Name Ld; Value {
          Integral { [ SymmetryFactor*Lz*Idir[]*NbWires[]/SurfCoil[]* CompZ[{a}] ] ;
            In Inds  ; Jacobian Vol ; Integration I1 ; }
        }
      }
    }
  } //Postprocessing: MagStaDyn_a_2D_Ld
}// PostProcessing

//-----------------------------------------------------------------------------------------------

PostOperation Get_LocalFields UsingPost MagStaDyn_a_2D {
  Print[ az, OnElementsOf Domain, File StrCat[Dir,StrCat["a",ExtGmsh]], LastTimeStepOnly ];
  Print[ dhdb, OnElementsOf Domain, LastTimeStepOnly, StoreInField 66 ];
}

PostOperation Get_GlobalQuantities UsingPost MagStaDyn_a_2D {
  Print[ Flux[Inds], OnGlobal, Format TimeTable,
    File > StrCat[Dir,StrCat["Flux",ExtGnuplot]], LastTimeStepOnly, Store 11,
    SendToServer StrCat[po,"40Flux [Wb]"],  Color "LightYellow" ];

  Print[ MagEnergy[Domain], OnGlobal, Format TimeTable,
    File > StrCat[Dir,StrCat["ME",ExtGnuplot]], LastTimeStepOnly, Store 22,
    SendToServer StrCat[po,"41Magnetic Energy [W]"],  Color "LightYellow" ];

  Print[ Inductance_from_Flux, OnRegion DomainDummy, Format Table, LastTimeStepOnly,
    File StrCat[Dir,StrCat["Inductance",ExtGnuplot]],
    SendToServer StrCat[po,"50Inductance from Flux [H]"], Color "LightYellow" ];

  Print[ Inductance_from_MagEnergy, OnRegion DomainDummy, Format Table, LastTimeStepOnly,
    File StrCat[Dir,StrCat["Inductance",ExtGnuplot]],
    SendToServer StrCat[po,"51Inductance from Magnetic Energy [H]"], Color "LightYellow" ];
}
PostOperation Get_GlobalQuantities_Ld UsingPost MagStaDyn_a_2D_Ld {
  Print[ Ld[Inds], OnGlobal, Format TimeTable,
         File StrCat[Dir,StrCat["Inductance_Ld",ExtGnuplot]],
         SendToServer StrCat[po,"55Inductance from Flux Diff[H]"], Color "LightYellow" ];
}  
 
DefineConstant[
//	Flux_ = {0, Name "Output/50Inductance from Flux [H]",ReadOnly 1, Graph "0200"},
//	Ld_ = {0, Name "Output/55Inductance from Flux Diff[H]", ReadOnly 1, Graph "0002"},
  R_ = {"Analysis", Name "GetDP/1ResolutionChoices", Visible 0},
  C_ = {"-solve -v2 -gmshread res/a.pos", Name "GetDP/9ComputeCommand", Visible 0},
  P_ = {"", Name "GetDP/2PostOperationChoices", Visible 0}
];
