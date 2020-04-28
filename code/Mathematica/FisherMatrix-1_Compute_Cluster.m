(* ::Package:: *)

(* ::Chapter:: *)
(*Structural complexity of functional response models - Part 1*)


(* ::Subsection:: *)
(*Define function to compute fisher matrix of a function given the vector of its parameters*)


(* ::Input::Initialization:: *)
fisher[fun_,pars_]:=-D[Apply[fun,pars],{pars,2}]


(* ::Subsection:: *)
(* Functional responses - Number of prey eaten*)


(* ::Input::Initialization:: *)
H1 = a \[CapitalNu] P T;
R = (a \[CapitalNu] P T)/P;
H2 = (a \[CapitalNu] P T)/(1+a h \[CapitalNu]);
HV = (a \[CapitalNu] P T)/P^m;
AG = (a \[CapitalNu] P T)/(P+a h \[CapitalNu]);
BD = (a \[CapitalNu] P T)/(1+a h \[CapitalNu] + \[Gamma] P);
CM = (a \[CapitalNu] P T)/((1+a h \[CapitalNu])(1+\[Gamma] P));
AA=(a \[CapitalNu] P T)/(P^m+a h \[CapitalNu]);


(* ::Subsection:: *)
(*Experimental Designs*)


(* ::Input::Initialization:: *)
subs={T-> 1};
\[CapitalNu]vals={10,15,25,45};
Pvals={1,2,3};

\[CapitalNu]probs=ConstantArray[1/Length[\[CapitalNu]vals],Length[\[CapitalNu]vals]];
Pprobs=ConstantArray[1/Length[Pvals],Length[Pvals]];
\[CapitalNu]DistE=EmpiricalDistribution[\[CapitalNu]probs->\[CapitalNu]vals];
PDistE=EmpiricalDistribution[Pprobs->Pvals];


(* ::Subsection:: *)
(*Define function for Expected Fisher Matrix given experimental design*)


(* ::Input::Initialization:: *)
EFisher[FisherMatrix_,func_]:=Expectation[Expectation[FisherMatrix,x\[Distributed]PoissonDistribution[func]],{P\[Distributed]PDistE,\[CapitalNu]\[Distributed]\[CapitalNu]DistE}]


(* ::Subsection:: *)
(*Likelihood and substitutions*)


(* ::Input::Initialization:: *)
PoislL[func_]:=-n func+Log[func] x/.{n->1}


(* ::Chapter:: *)
(*Fisher information matrices and determinants*)


(* ::Section:: *)
(*One-parameter models*)


(* ::Subsubsection:: *)
(*Holling Type I model - Poisson*)


(* ::Input::Initialization:: *)
lLH1Poiss[a_]:=PoislL[H1]
FisherH1Poiss=fisher[lLH1Poiss,{a}];
EFisherH1Poiss=EFisher[FisherH1Poiss/.subs,H1/.subs];
DetH1=Det[EFisherH1Poiss];


(* ::Subsubsection:: *)
(*Ratio - Poisson*)


(* ::Input::Initialization:: *)
lLRPoiss[a_]:=PoislL[R]
FisherRPoiss=fisher[lLRPoiss,{a}]//FullSimplify;
EFisherRPoiss=EFisher[FisherRPoiss/.subs,R/.subs];
DetR=Det[EFisherRPoiss]//FullSimplify;


(* ::Section:: *)
(*Two-parameter models*)


(* ::Subsubsection:: *)
(*Holling Type II model - Poisson*)


(* ::Input::Initialization:: *)
lLH2Poiss[a_,h_]:=PoislL[H2]
FisherH2Poiss=fisher[lLH2Poiss,{a, h}];
EFisherH2Poiss=EFisher[FisherH2Poiss/.subs,H2/.subs];
DetH2=Det[EFisherH2Poiss]//FullSimplify;


(* ::Subsubsection:: *)
(*Hassell-Varley - Poisson*)


(* ::Input::Initialization:: *)
lLHVPoiss[a_,m_]:=PoislL[HV]
FisherHVPoiss=fisher[lLHVPoiss,{a,m}]//FullSimplify;
EFisherHVPoiss=EFisher[FisherHVPoiss/.subs,HV/.subs];
DetHV=Det[EFisherHVPoiss];


(* ::Subsubsection:: *)
(*Arditi-Ginzburg - Poisson*)


(* ::Input::Initialization:: *)
lLAGPoiss[a_,h_]:=PoislL[AG]
FisherAGPoiss=fisher[lLAGPoiss,{a, h}]//FullSimplify;
EFisherAGPoiss=EFisher[FisherAGPoiss/.subs,AG/.subs];
DetAG=Det[EFisherAGPoiss];


(* ::Section:: *)
(*Three-parameter models*)


(* ::Subsubsection:: *)
(*Crowley-Martin model - Poisson*)


(* ::Input::Initialization:: *)
lLCMPoiss[a_,h_,\[Gamma]_]:=PoislL[CM]
FisherCMPoiss=fisher[lLCMPoiss,{a, h,\[Gamma]}]//FullSimplify;
EFisherCMPoiss=EFisher[FisherCMPoiss/.subs,CM/.subs];
DetCM=Det[EFisherCMPoiss];


(* ::Subsubsection:: *)
(*Beddington-DeAngelis model - Poisson*)


(* ::Input::Initialization:: *)
lLBDPoiss[a_,h_,\[Gamma]_]:=PoislL[BD]
FisherBDPoiss=fisher[lLBDPoiss,{a, h,\[Gamma]}]//FullSimplify;
EFisherBDPoiss=EFisher[FisherBDPoiss/.subs,BD/.subs];
DetBD=Det[EFisherBDPoiss];


(* ::Subsubsection:: *)
(*Arditi-Akcakaya model - Poisson*)


(* ::Input::Initialization:: *)
lLAAPoiss[a_,h_,m_]:=PoislL[AA]
FisherAAPoiss=fisher[lLAAPoiss,{a, h,m}]//FullSimplify;
EFisherAAPoiss=EFisher[FisherAAPoiss/.subs,AA/.subs];
DetAA=Det[EFisherAAPoiss];


(* ::Text:: *)
(*Override save all results (for all models)*)


(* ::Input::Initialization:: *)
Save["FisherMatrix-1_Output_k1k2k3","`*"];


(* ::Chapter:: *)
(*Integrate square root of Expected Fisher Information Matrix over parameters, and take natural log*)


(* ::Subsection:: *)
(*Numeric Integration*)


(* ::Input::Initialization:: *)
P\[CapitalNu]rats=Sort[DeleteDuplicates[Join@@Outer[Divide,Pvals,\[CapitalNu]vals]]];

e=10;
ParmRangeA={
{a,0,1}, (* Note maximum!!!! *)
Flatten[{h,0,P\[CapitalNu]rats,10^e}],
Flatten[{\[Gamma],0,1/Pvals,10^e}],
Flatten[{m,0,P\[CapitalNu]rats,10^e}]
};
precgoal=2;


(* ::Subsubsection:: *)
(*One-parameter models*)


(* ::Input::Initialization:: *)
NIntH1=Log[NIntegrate[Sqrt[DetH1],
ParmRangeA[[1]],
PrecisionGoal->precgoal]];


(* ::Input::Initialization:: *)
NIntR=Log[NIntegrate[Sqrt[DetR],
ParmRangeA[[1]],
PrecisionGoal->precgoal]];


(* ::Subsubsection:: *)
(*Two-parameter models*)


(* ::Input::Initialization:: *)
NIntHV=Log[NIntegrate[Sqrt[DetHV],
ParmRangeA[[1]],
ParmRangeA[[4]],
PrecisionGoal->precgoal]];


(* ::Input::Initialization:: *)
NIntH2=Log[NIntegrate[Sqrt[DetH2],
ParmRangeA[[1]],
ParmRangeA[[2]],
PrecisionGoal->precgoal]];


(* ::Input::Initialization:: *)
NIntAG=Log[NIntegrate[Sqrt[DetAG],
ParmRangeA[[1]],
ParmRangeA[[2]],
PrecisionGoal->precgoal,Method->{"GlobalAdaptive","SingularityHandler"->"DuffyCoordinates"},
MaxRecursion->500]];


(* ::Subsubsection:: *)
(*Three-parameter models*)


(* ::Input::Initialization:: *)
NIntBD=Log[NIntegrate[Sqrt[DetBD],
ParmRangeA[[1]],
ParmRangeA[[2]],
ParmRangeA[[3]],
PrecisionGoal->precgoal,Method->{"GlobalAdaptive","SingularityHandler"->"DuffyCoordinates"},
MaxRecursion->500]];


(* ::Input::Initialization:: *)
NIntCM=Log[NIntegrate[Sqrt[DetCM],
ParmRangeA[[1]],
ParmRangeA[[2]],
ParmRangeA[[3]],
PrecisionGoal->3,Method->{"GlobalAdaptive","SingularityHandler"->"DuffyCoordinates"},
MaxRecursion->500]];


(* ::Input::Initialization:: *)
NIntAA=Log[NIntegrate[Sqrt[DetAA],
ParmRangeA[[1]],
ParmRangeA[[2]],
ParmRangeA[[4]],
PrecisionGoal->precgoal,Method->{"GlobalAdaptive","SingularityHandler"->"DuffyCoordinates"},
MaxRecursion->500]];


(* ::Text:: *)
(*Save all results (for all models)*)


(* ::Input::Initialization:: *)
Save["FisherMatrix-1_Output_k1k2k3_NInt","`*"];


(* ::Subsection:: *)
(*Partial numeric Integration*)


(* ::Input::Initialization:: *)
ParmRangeB=Drop[ParmRangeA,1];
precgoal=2;


(* ::Input::Initialization:: *)
aNIntH1=Sqrt[DetH1];
aNIntR=Sqrt[DetR];
aNIntHV=NIntegrate[Sqrt[DetHV],
ParmRangeB[[3]],
PrecisionGoal->precgoal];
aNIntH2[\[Alpha]_?NumericQ]:=NIntegrate[Sqrt[DetH2/.{a->\[Alpha]}],
ParmRangeB[[1]],
PrecisionGoal->precgoal]
aNIntAG[\[Alpha]_?NumericQ]:=NIntegrate[Sqrt[DetAG/.{a->\[Alpha]}],
ParmRangeB[[1]],
PrecisionGoal->precgoal];
aNIntBD[\[Alpha]_?NumericQ]:=NIntegrate[Sqrt[DetBD/.{a->\[Alpha]}],
ParmRangeB[[1]],ParmRangeB[[2]],
PrecisionGoal->precgoal];
aNIntCM[\[Alpha]_?NumericQ]:=NIntegrate[Sqrt[DetCM/.{a->\[Alpha]}],
ParmRangeB[[1]],ParmRangeB[[2]],
PrecisionGoal->precgoal];
aNIntAA[\[Alpha]_?NumericQ]:=NIntegrate[Sqrt[DetAA/.{a->\[Alpha]}],
ParmRangeB[[1]],ParmRangeB[[3]],
PrecisionGoal->precgoal];


(* ::Text:: *)
(*Override save all results (for all models)*)


(* ::Input::Initialization:: *)
Save["FisherMatrix-1_Output_k1k2k3_aNInt","`*"];


(* ::Subsection:: *)
(*Investigate limits*)


(* ::Input::Initialization:: *)
LimsH1a={a,Limit[Sqrt[DetH1],{a->0}], Limit[Sqrt[DetH1],{a->\[Infinity]}]}
LimsRa={a,Limit[Sqrt[DetR],{a->0}], Limit[Sqrt[DetR],{a->\[Infinity]}]}
LimsHVa={a,Limit[Sqrt[DetHV],{a->0}], Limit[Sqrt[DetHV],{a->\[Infinity]}]}
LimsH2a={a,Limit[Sqrt[DetH2],{a->0}], Limit[Sqrt[DetH2],{a->\[Infinity]}]}
LimsAGa={a,Limit[Sqrt[DetAG],{a->0}],Limit[Sqrt[DetAG],{a->\[Infinity]}]}
LimsBDa={a,Limit[Sqrt[DetBD],{a->0}], Limit[Sqrt[DetBD],{a->\[Infinity]}]}
LimsCMa={a,Limit[Sqrt[DetCM],{a->0}],Limit[Sqrt[DetCM],{a->\[Infinity]}]}
LimsAAa={a,Limit[Sqrt[DetAA],{a->0}], Limit[Sqrt[DetAA],{a->\[Infinity]}]}


(* ::Input::Initialization:: *)
LimsH2h={h,Limit[Sqrt[DetH2],{h->0}], Limit[Sqrt[DetH2],{h->\[Infinity]}]}
LimsAGh={h,Limit[Sqrt[DetAG],{h->0}],Limit[Sqrt[DetAG],{h->\[Infinity]}]}
LimsBDh={h,Limit[Sqrt[DetBD],{h->0}], Limit[Sqrt[DetBD],{h->\[Infinity]}]}
LimsCMh={h,Limit[Sqrt[DetCM],{h->0}],Limit[Sqrt[DetCM],{h->\[Infinity]}]}
LimsAAh={h,Limit[Sqrt[DetAA],{h-> 0}],Limit[Sqrt[DetAA],{h->\[Infinity]}]}


(* ::Input::Initialization:: *)

LimsBD\[Gamma]={\[Gamma],Limit[Sqrt[DetBD],{\[Gamma]->0}], Limit[Sqrt[DetBD],{\[Gamma]->\[Infinity]}]}
LimsCM\[Gamma]={\[Gamma],Limit[Sqrt[DetCM],{\[Gamma]->0}],Limit[Sqrt[DetCM],{\[Gamma]->\[Infinity]}]}


(* ::Input::Initialization:: *)
LimsHVm={m,Limit[Sqrt[DetHV],{m->0}], Limit[Sqrt[DetHV],{m->\[Infinity]}]}
LimsAAm={m,Limit[Sqrt[DetAA],{m->0}],Limit[Sqrt[DetAA],{m->\[Infinity]}]}


(* ::Text:: *)
(*Override save all results (for all models)*)


(* ::Input::Initialization:: *)
Save["FisherMatrix-1_Output_k1k2k3_limits","`*"];