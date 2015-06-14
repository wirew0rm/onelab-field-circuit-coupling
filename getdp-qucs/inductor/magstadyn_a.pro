Group {
   DefineGroup[
     DomainM, DomainB, DomainS, DomainInf,
     DomainL, DomainNL,
     Surf_bn0, Surf_Inf
  ];
}

Function{
 DefineConstant[
    Val_Rint, Val_Rext,
    Lz = 1,
    SymmetryFactor = 1,
    Nb_max_iter = 30,
    relaxation_factor = 1,
    stop_criterion = 1e-5,
    reltol = 1e-7,
    abstol = 1e-5,
    T = 1/Freq, // Fundamental period in s
    time0 = 0,
    NbT = 1,
    timemax = NbT*T,
    NbSteps = 100,
    delta_time = T/NbSteps,
    FillFactor_Winding = {1, Name "Input/4Coil Parameters/3Fill factor",
      Highlight "AliceBlue", Visible 1},
    Factor_R_3DEffects = {1, Name "Input/4Coil Parameters/9fact", Label "3D factor",
      Highlight "AliceBlue", Visible 1},
    // Increasing the resistance by 50% == 1.5
    II, VV,
    Flag_NL = 0,
    Flag_NL_Newton_Raphson = {1, Choices{0,1}, Name "Input/41Newton-Raphson iteration",
      Visible Flag_NL},
    po = "Output/"
		Flag_InitFromPrevious = {0, Choices{0,1}, Name "Input/42InitSolutionFromPrevious"}
  ];

  DefineFunction[
    dhdb_NL, dhdb, br, js0
  ];
  //=====================================================================================
  //===================================================================================== 
  dhdb_Diff[] = TensorField[XYZ[]]{0} ;
  //=====================================================================================
  //===================================================================================== 
}

Include "BH.pro"; // nonlinear BH caracteristic of magnetic material

Group {
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

Function {
  nu [#{Air, AirInf, Inds}]  = TensorDiag[1.0, 1.0, 1.0] * 1./mu0 ;

  If(!Flag_NL)
    nu [#{Core}]  = 1/(mur_fe*mu0) ;
  EndIf
  If(Flag_NL)
    nu [ DomainNL ] = TensorDiag[1.0, 1.0, 1.0] * nu_EIcore[$1] ;
    dhdb_NL [ DomainNL ] = dhdb_EIcore_NL[$1];
    dhdb [ DomainNL ] = dhdb_EIcore[$1];
  EndIf

  sigma[#{Inds}] = sigma_coil ;
  rho[] = 1/sigma[] ;

  Rb[] = Factor_R_3DEffects*FillFactor_Winding*Lz*NbWires[]^2/SurfCoil[]/sigma[] ;
  Resistance[#{Inds}] = Rb[] ;
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
      { Region Surf_Inf ; Type Assign ; Value 0. ; }
      { Region Surf_bn0 ; Type Assign ; Value 0. ; }
    }
  }

	// ====================================================
	// Initialize a with the value from the previous run.
	// getdp has to be called with -gmshread res/a.pos flag
	// ====================================================
	{ Name InitA ;
		Case {
			{ Region Domain ; Type Init ; Value Field[XYZ[]] ; }
		}
	}
	// ====================================================

  { Name Current_2D ;
    Case {
      { Region Inds ; Value II*Idir[] ; TimeFunction IA[]; }
    }
  }

  { Name Voltage_2D ;
    Case {
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
      { NameOfCoef ae1 ; EntityType NodesOf ; NameOfConstraint MVP_2D ; }
			// == Apply initial Values to vector Potential ==
			If (Flag_InitFromPrevious)
				{ NameOfCoef ae1 ; EntityType NodesOf ; NameOfConstraint InitA ; }
			EndIf
    }
  }

  // Gradient of Electric scalar potential (2D)
  { Name Hregion_u_Mag_2D ; Type Form1P ;
    BasisFunction {
      { Name sr ; NameOfCoef ur ; Function BF_RegionZ ;
        Support DomainC ; Entity DomainC ; }
    }
    GlobalQuantity {
      { Name U ; Type AliasOf        ; NameOfCoef ur ; }
      { Name I ; Type AssociatedWith ; NameOfCoef ur ; }
    }
    Constraint {
      { NameOfCoef U ; EntityType GroupsOfNodesOf ; NameOfConstraint Voltage_2D ; }
      { NameOfCoef I ; EntityType GroupsOfNodesOf ; NameOfConstraint Current_2D ; }
    }
  }

  { Name Hregion_i_Mag_2D ; Type Vector ;
    BasisFunction {
      { Name sr ; NameOfCoef ir ; Function BF_RegionZ ;
        Support DomainB ; Entity DomainB ; }
    }
    GlobalQuantity {
      { Name Ib ; Type AliasOf        ; NameOfCoef ir ; }
      { Name Ub ; Type AssociatedWith ; NameOfCoef ir ; }
    }
    Constraint {
      { NameOfCoef Ub ; EntityType Region ; NameOfConstraint Voltage_2D ; }
      { NameOfCoef Ib ; EntityType Region ; NameOfConstraint Current_2D ; }
    }
  }

}

//-----------------------------------------------------------------------------------------------

Formulation {

  { Name MagStaDyn_a_2D ; Type FemEquation ;
    Quantity {
      { Name a  ; Type Local  ; NameOfSpace Hcurl_a_2D ; }
      { Name ur ; Type Local  ; NameOfSpace Hregion_u_Mag_2D ; }
      { Name I  ; Type Global ; NameOfSpace Hregion_u_Mag_2D [I] ; }
      { Name U  ; Type Global ; NameOfSpace Hregion_u_Mag_2D [U] ; }
      
      { Name ir ; Type Local  ; NameOfSpace Hregion_i_Mag_2D ; }
      { Name Ub ; Type Global ; NameOfSpace Hregion_i_Mag_2D [Ub] ; }
      { Name Ib ; Type Global ; NameOfSpace Hregion_i_Mag_2D [Ib] ; }
    }

    Equation {
      Galerkin { [ nu[{d a}] * Dof{d a}  , {d a} ] ;
        In Domain ; Jacobian Vol ; Integration I1 ; }
      If(Flag_NL_Newton_Raphson)
      Galerkin { JacNL [ dhdb_NL[{d a}] * Dof{d a} , {d a} ] ;
        In DomainNL ; Jacobian Vol ; Integration I1 ; }
      EndIf

      Galerkin { [ -nu[] * br[] , {d a} ] ;
        In DomainM ; Jacobian Vol ; Integration I1 ; }

      Galerkin { DtDof[ sigma[] * Dof{a} , {a} ] ;
        In DomainC ; Jacobian Vol ; Integration I1 ; }
      Galerkin { [ sigma[] * Dof{ur}, {a} ] ;
        In DomainC ; Jacobian Vol ; Integration I1 ; }

      Galerkin { [ -js0[] , {a} ] ;
        In DomainS ; Jacobian Vol ; Integration I1 ; }

      Galerkin { DtDof[ sigma[] * Dof{a} , {ur} ] ;
        In DomainC ; Jacobian Vol ; Integration I1 ; }
      Galerkin { [ sigma[] * Dof{ur} , {ur} ] ;
        In DomainC ; Jacobian Vol ; Integration I1 ; }
      GlobalTerm { [ Dof{I} , {U} ] ; In DomainC ; }
        
      Galerkin { [ -NbWires[]/SurfCoil[] * Dof{ir} , {a} ] ;
        In DomainB ; Jacobian Vol ; Integration I1 ; }
      Galerkin { DtDof [ Lz * NbWires[]/SurfCoil[] * Dof{a} , {ir} ] ;
        In DomainB ; Jacobian Vol ; Integration I1 ; }
      GlobalTerm { [ Dof{Ub}/SymmetryFactor , {Ib} ] ; In DomainB ; }
      Galerkin { [ Rb[]/SurfCoil[]* Dof{ir} , {ir} ] ;
        In DomainB ; Jacobian Vol ; Integration I1 ; } // Resistance term
    }
  }

  //=====================================================================================
  //===================================================================================== 
  { Name MagStaDyn_a_2D_Diff ; Type FemEquation ;
    Quantity {
      { Name a  ; Type Local  ; NameOfSpace Hcurl_a_2D ; }

      { Name ir ; Type Local  ; NameOfSpace Hregion_i_Mag_2D ; }
      { Name Ub ; Type Global ; NameOfSpace Hregion_i_Mag_2D [Ub] ; }
      { Name Ib ; Type Global ; NameOfSpace Hregion_i_Mag_2D [Ib] ; }
    }

    Equation {
      Galerkin { [ dhdb_Diff[] * Dof{d a}  , {d a} ] ;
        In Domain ; Jacobian Vol ; Integration I1 ; }

      Galerkin { [ -NbWires[]/SurfCoil[] * Dof{ir} , {a} ] ;
        In DomainB ; Jacobian Vol ; Integration I1 ; }
      Galerkin { DtDof [ Lz * NbWires[]/SurfCoil[] * Dof{a} , {ir} ] ;
        In DomainB ; Jacobian Vol ; Integration I1 ; }
      GlobalTerm { [ Dof{Ub}/SymmetryFactor , {Ib} ] ; In DomainB ; }
      Galerkin { [ Rb[]/SurfCoil[]* Dof{ir} , {ir} ] ;
        In DomainB ; Jacobian Vol ; Integration I1 ; } // Resistance term
    }
  }
  //=====================================================================================
  //===================================================================================== 

}

//-----------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------

Resolution {

  { Name Analysis ;
    System {
      If(Flag_AnalysisType==2)
        { Name A ; NameOfFormulation MagStaDyn_a_2D ; Type ComplexValue ; Frequency Freq ; }
      EndIf
      If(Flag_AnalysisType<2)
      { Name A ; NameOfFormulation MagStaDyn_a_2D ; }
      //=====================================================================================
      //===================================================================================== 
      { Name A_Diff ; NameOfFormulation MagStaDyn_a_2D_Diff ; }
      //=====================================================================================
      //===================================================================================== 
      EndIf
    }
    Operation {
      CreateDir["res/"];
      InitSolution[A] ;
      PostOperation[Get_Analytical] ; // Values from magnetic circuit
      If(Flag_AnalysisType==0 || Flag_AnalysisType==2)
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
      EndIf // End Flag_AnalysisType==0 (Static) Flag_AnalysisType==2 (Frequency)

        //=====================================================================================
        //=====================================================================================
        //=====================================================================================
        //=====================================================================================
        //=====================================================================================
        //===================================================================================== 
      CreateDir["res_diff/"] ;
      GmshRead["res/dhdb.pos"] ;
      InitSolution[A_Diff] ;
      //If(Flag_NL)
      //  IterativeLoop[Nb_max_iter, stop_criterion, relaxation_factor]{
      //GenerateJac[A_Diff] ; SolveJac[A_Diff] ; }
      //EndIf
      Generate[A_Diff] ; Solve[A_Diff]; SaveSolution[A_Diff] ;
      PostOperation[Get_LocalFields_Diff] ;
      PostOperation[Get_GlobalQuantities_Diff] ;  
      //=====================================================================================
      //=====================================================================================
      //===================================================================================== 
          
      If(Flag_AnalysisType==1)
        TimeLoopTheta[time0, timemax, delta_time, 1.]{ // Euler implicit (1) -- Crank-Nicolson (0.5)
          If(!Flag_NL)
            Generate[A]; Solve[A];
          EndIf
          If(Flag_NL)
            IterativeLoop[Nb_max_iter, stop_criterion, relaxation_factor] {
              GenerateJac[A] ; SolveJac[A] ; }
          EndIf
          SaveSolution[A];

          PostOperation[Get_LocalFields] ;
          Test[ $TimeStep > 1 ]{
            PostOperation[Get_GlobalQuantities];
          }
        }
      EndIf // End Flag_AnalysisType==1 (Time domain)
    }
  }
}

//-----------------------------------------------------------------------------------------------

PostProcessing {

  { Name MagStaDyn_a_2D ; NameOfFormulation MagStaDyn_a_2D ;
    PostQuantity {
      { Name a  ; Value { Term { [ {a} ] ; In Domain ; Jacobian Vol ; } } }
      { Name az ; Value { Term { [ CompZ[{a}] ] ; In Domain ; Jacobian Vol ; } } }

      { Name b  ; Value { Term { [ {d a} ] ; In Domain ; Jacobian Vol ; } } }
      { Name nb  ; Value { Term { [ Norm[{d a}] ] ; In Domain ; Jacobian Vol ; } } }
      { Name br ; Value { Term { [ br[] ] ; In DomainM ; Jacobian Vol ; } } }

      { Name h ; Value { Term { [ nu[{d a}] * {d a} ] ; In Domain ; Jacobian Vol ; } } }

      { Name js0 ; Value { Term { [ js0[] ] ; In DomainS ; Jacobian Vol ; } } }

      { Name j  ; Value {
          Term { [ -sigma[]*(Dt[{a}]+{ur}) ]        ; In DomainC ; Jacobian Vol ; }
        }
      }

      { Name ir ; Value { Term { [ {ir} ] ; In Inds ; Jacobian Vol ; } } }

      { Name jz ; Value {
          Term { [ CompZ[-sigma[]*(Dt[{a}]+{ur})] ]       ; In DomainC ; Jacobian Vol ; }
        }
      }

      { Name dhdb ; Value {
          Term { [ nu[{d a}]]; In DomainL; Jacobian Vol; }
          Term { [ dhdb[{d a}] ] ; In DomainNL ; Jacobian Vol ; }
        }
      }

      { Name rhoj2 ;
        Value {
          Term { [ sigma[]*SquNorm[ Dt[{a}]+{ur}] ] ; In Region[{DomainC}] ; Jacobian Vol ; }
          Term { [ 1./sigma[]*SquNorm[ IA[]*{ir} ] ] ; In Inds  ; Jacobian Vol ; }
        }
      }

      { Name JouleLosses ;
        Value {
          Integral { [ SymmetryFactor*Lz*sigma[] * SquNorm[ Dt[{a}]+{ur} ] ];
            In Region[{DomainC}] ; Jacobian Vol ; Integration I1 ; }
          Integral { [ SymmetryFactor*Lz*1./sigma[]*SquNorm[ IA[]*{ir} ] ];
            In Inds  ; Jacobian Vol ; Integration I1 ; }
        }
      }

      { Name MagEnergy ; Value {
          Integral { [ SymmetryFactor*Lz* 1/2 *nu[{d a}]*{d a}*{d a} ] ;
            In Domain ; Jacobian Vol ; Integration I1 ; } } }

      { Name diffMagEnergy ; Value {
          Integral { [ SymmetryFactor*Lz* 1/2 *nu[{d a}]*{d a}*{d a} ] ;
            In Domain ; Jacobian Vol ; Integration I1 ; }
          Integral { [ SymmetryFactor*Lz* 1/2 *dhdb_NL[{d a}]*{d a}*{d a} ] ;
            In DomainNL ; Jacobian Vol ; Integration I1 ; }
						}}
			
      { Name Flux ; Value {
          Integral { [ SymmetryFactor*Lz*Idir[]*NbWires[]/SurfCoil[]* CompZ[{a}] ] ;
            In Inds  ; Jacobian Vol ; Integration I1 ; }
        }
      }

      { Name ComplexPower ; // S = P + i*Q
        Value {
          Integral { [ Complex[ sigma[]*SquNorm[Dt[{a}]+{ur}], nu[]*SquNorm[{d a}] ] ] ;
            In Region[{DomainC}] ; Jacobian Vol ; Integration I1 ; }
        }
      }

      { Name U ; Value {
          Term { [ {U} ]   ; In DomainC ; }
          Term { [ {Ub} ]  ; In DomainB ; }
        }
      }

      { Name I ; Value {
          Term { [ {I} ]   ; In DomainC ; }
          Term { [ {Ib} ]  ; In DomainB ; }
        }
      }

      { Name S ; Value {
          Term { [ {U}*Conj[{I}] ]    ; In DomainC ; }
          Term { [ {Ub}*Conj[{Ib}] ]  ; In DomainB ; }
        }
      }

      // Getting the value of some functions
     For k In {0:NbAvailableMagCircuits-1}
       { Name Reluctance~{k} ; Value { Term { Type Global; [ Reluctance~{k}[] ] ; In DomainDummy ; } } }
       { Name Inductance~{k} ; Value { Term { Type Global; [ Inductance~{k}[] ] ; In DomainDummy ; } } }
     EndFor

      { Name Inductance_from_Flux ; Value { Term { Type Global; [ #11/II ] ; In DomainDummy ; } } } // Flux stored in register #11
      { Name Inductance_from_MagEnergy ; Value { Term { Type Global; [ 2*#22/(II*II) ] ; In DomainDummy ; } } } // MagEnergy stored in register #22
      { Name Ld ; Value { Term { Type Global; [ 2*#23/(II*II) ] ; In DomainDummy ; } } } // MagEnergy stored in register #22

    }//PostQuantity
  }// MagStaDyn_a_2D

  { Name MagStaDyn_a_2D_Diff ; NameOfFormulation MagStaDyn_a_2D_Diff ;
    PostQuantity {
      { Name dhdb_Diff ; Value {
          Term { [ dhdb_Diff[{d a}] ] ; In Domain ; Jacobian Vol ; }
        }
      } 
      { Name Flux_Diff ; Value {
          Integral { [ SymmetryFactor*Lz*Idir[]*NbWires[]/SurfCoil[]* CompZ[{a}] ] ;
            In Inds  ; Jacobian Vol ; Integration I1 ; }
        }
      }
      { Name Inductance_from_Flux_Diff ; Value { Term { Type Global; [ #55/II ] ; In DomainDummy ; } } }

    }
  }
}// PostProcessing


//-----------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------

PostOperation Get_LocalFields UsingPost MagStaDyn_a_2D {
  Print[ ir, OnElementsOf Inds,   File StrCat[Dir,StrCat["ir",ExtGmsh]], LastTimeStepOnly ] ;
  Print[ b,  OnElementsOf Domain, File StrCat[Dir,StrCat["b",ExtGmsh]], LastTimeStepOnly ] ;
  Print[ nb,  OnElementsOf Domain, File StrCat[Dir,StrCat["nb",ExtGmsh]], LastTimeStepOnly ] ;
  Print[ az, OnElementsOf Domain, File StrCat[Dir,StrCat["a",ExtGmsh]], LastTimeStepOnly ];
  Print[ dhdb, OnElementsOf Domain, File StrCat[Dir,StrCat["dhdb",ExtGmsh]], LastTimeStepOnly ];
}

PostOperation Get_Analytical UsingPost MagStaDyn_a_2D {
  For k In {0:NbAvailableMagCircuits-1}
    Print[ Reluctance~{k}, OnRegion DomainDummy, Format Table, LastTimeStepOnly,
      File StrCat[Dir,StrCat[Sprintf("Reluctance%g",k),ExtGnuplot]] ];
    Print[ Inductance~{k}, OnRegion DomainDummy, Format Table, LastTimeStepOnly,
      File StrCat[Dir, StrCat[Sprintf("Inductance%g",k),ExtGnuplot]],
      SendToServer StrCat[po,Sprintf("6%gInductance Magnetic Circuit %g [mH]", k, k)], Color "LightYellow" ];
  EndFor
}

PostOperation Get_GlobalQuantities UsingPost MagStaDyn_a_2D {
  Print[ I, OnRegion Ind_1, Format Table,
    File > StrCat[Dir,StrCat["I",ExtGnuplot]], LastTimeStepOnly,
    SendToServer StrCat[po,"20I [A]"], Color "LightYellow" ];
  Print[ U, OnRegion Ind_1, Format Table,
    File > StrCat[Dir,StrCat["U",ExtGnuplot]], LastTimeStepOnly,
    SendToServer StrCat[po,"30U [V]"], Color "LightYellow" ];

  Print[ Flux[Inds], OnGlobal, Format TimeTable,
    File > StrCat[Dir,StrCat["Flux",ExtGnuplot]], LastTimeStepOnly, Store 11,
    SendToServer StrCat[po,"40Flux [Wb]"],  Color "LightYellow" ];

  Print[ MagEnergy[Domain], OnGlobal, Format TimeTable,
    File > StrCat[Dir,StrCat["ME",ExtGnuplot]], LastTimeStepOnly, Store 22,
    SendToServer StrCat[po,"41Magnetic Energy [W]"],  Color "LightYellow" ];

  Print[ diffMagEnergy[Domain], OnGlobal, Format TimeTable,
    File > StrCat[Dir,StrCat["dME",ExtGnuplot]], LastTimeStepOnly, Store 23,
    SendToServer StrCat[po,"42diffMagnetic Energy [W]"],  Color "LightYellow" ];

  Print[ Ld, OnRegion DomainDummy, Format Table, LastTimeStepOnly,
    File StrCat[Dir,StrCat["Ld",ExtGnuplot]],
    SendToServer StrCat[po,"53Ld [H]"], Color "LightYellow" ]; // export Ld

  Print[ Inductance_from_Flux, OnRegion DomainDummy, Format Table, LastTimeStepOnly,
    File StrCat[Dir,StrCat["Inductance",ExtGnuplot]],
    SendToServer StrCat[po,"50Inductance from Flux [H]"], Color "LightYellow" ];
  Print[ Inductance_from_MagEnergy, OnRegion DomainDummy, Format Table, LastTimeStepOnly,
    File StrCat[Dir,StrCat["Inductance",ExtGnuplot]],
    SendToServer StrCat[po,"51Inductance from Magnetic Energy [H]"], Color "LightYellow" ];
}


 PostOperation Get_LocalFields_Diff UsingPost MagStaDyn_a_2D_Diff {
   Print[ dhdb_Diff, OnElementsOf Domain, File StrCat[Dir,StrCat["dhdb_Diff",ExtGmsh]], LastTimeStepOnly ]; 
 }
PostOperation Get_GlobalQuantities_Diff UsingPost MagStaDyn_a_2D_Diff {
  Print[ Flux_Diff[Inds], OnGlobal, Format TimeTable,
         File > StrCat[Dir,StrCat["Flux_Diff",ExtGnuplot]], LastTimeStepOnly, Store 55,
         SendToServer StrCat[po,"44Flux Diff[Wb]"],  Color "LightYellow" ];
  Print[ Inductance_from_Flux_Diff, OnRegion DomainDummy, Format Table, LastTimeStepOnly,
         File StrCat[Dir,StrCat["Inductance_Diff",ExtGnuplot]],
         SendToServer StrCat[po,"55Inductance from Flux Diff[H]"], Color "LightYellow" ];
}  
 
DefineConstant[
  R_ = {"Analysis", Name "GetDP/1ResolutionChoices", Visible 0},
  C_ = {"-solve -v2", Name "GetDP/9ComputeCommand", Visible 0},
  P_ = {"", Name "GetDP/2PostOperationChoices", Visible 0}
];
