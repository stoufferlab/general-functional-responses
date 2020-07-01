(* ::Package:: *)

(* ::Input::Initialization:: *)
fisher[fun_,pars_]:=-D[Apply[fun,pars],{pars,2}]

H1 = a \[CapitalNu] P T;
R = (a \[CapitalNu] P T)/P;
H2 = (a \[CapitalNu] P T)/(1+a h \[CapitalNu]);
HV = (a \[CapitalNu] P T)/P^m;
AG = (a \[CapitalNu] P T)/(P+a h \[CapitalNu]);
BD = (a \[CapitalNu] P T)/(1+a h \[CapitalNu] + \[Gamma] (P-1));
CM = (a \[CapitalNu] P T)/((1+a h \[CapitalNu])(1+\[Gamma] (P-1)));
AA=(a \[CapitalNu] P T)/(P^m+a h \[CapitalNu]);
\[CapitalNu]vals={2,4,8,16,32,64,128,256};
Pvals={1,2,5};

\[CapitalNu]probs=ConstantArray[1/Length[\[CapitalNu]vals],Length[\[CapitalNu]vals]];
Pprobs=ConstantArray[1/Length[Pvals],Length[Pvals]];
\[CapitalNu]DistE=EmpiricalDistribution[\[CapitalNu]probs->\[CapitalNu]vals];
PDistE=EmpiricalDistribution[Pprobs->Pvals];

subs={T-> 1};

EFisher[FisherMatrix_,func_]:=
Expectation[
Expectation[FisherMatrix,x\[Distributed] PoissonDistribution[func]],
{P\[Distributed]PDistE, \[CapitalNu]\[Distributed]\[CapitalNu]DistE}] 

DetEFisherTrunc[FisherMatrix_,func_,\[Lambda]min_,\[Lambda]max_,subs_]:=
Det[
Expectation[
Expectation[Boole[\[Lambda]min <= func<=\[Lambda]max]FisherMatrix,x\[Distributed] PoissonDistribution[func]]//.subs,
{P\[Distributed]PDistE, \[CapitalNu]\[Distributed]\[CapitalNu]DistE}] 
]

PoislL[func_]:=-n func+Log[func] x/.{n->1}


(* ::Input::Initialization:: *)
lLH1Poiss[a_]:=PoislL[H1]
FisherH1Poiss=fisher[lLH1Poiss,{a}];

EFisherH1Poiss=EFisher[FisherH1Poiss,H1]/.subs;

DetH1=Det[EFisherH1Poiss];
DetH1Trunc=DetEFisherTrunc[FisherH1Poiss,H1,0,Max[\[CapitalNu]vals],subs];

lLRPoiss[a_]:=PoislL[R]
FisherRPoiss=fisher[lLRPoiss,{a}]//FullSimplify;

EFisherRPoiss=EFisher[FisherRPoiss,R]//.subs;

DetR=Det[EFisherRPoiss]//FullSimplify;
DetRTrunc=DetEFisherTrunc[FisherRPoiss,R,0,Max[\[CapitalNu]vals],subs];

lLHVPoiss[a_,m_]:=PoislL[HV]
FisherHVPoiss=fisher[lLHVPoiss,{a,m}]//FullSimplify;

EFisherHVPoiss=EFisher[FisherHVPoiss,HV]//.subs;

DetHV=Det[EFisherHVPoiss];
DetHVTrunc=DetEFisherTrunc[FisherHVPoiss,HV,0,Max[\[CapitalNu]vals],subs];

lLH2Poiss[a_,h_]:=PoislL[H2]
FisherH2Poiss=fisher[lLH2Poiss,{a, h}];

EFisherH2Poiss=EFisher[FisherH2Poiss,H2]//.subs;

DetH2=Det[EFisherH2Poiss]//FullSimplify;
DetH2Trunc=DetEFisherTrunc[FisherH2Poiss,H2,0,Max[\[CapitalNu]vals],subs];

lLAGPoiss[a_,h_]:=PoislL[AG]
FisherAGPoiss=fisher[lLAGPoiss,{a, h}]//FullSimplify;

EFisherAGPoiss=EFisher[FisherAGPoiss,AG]//.subs;

DetAG=Det[EFisherAGPoiss];
DetAGTrunc=DetEFisherTrunc[FisherAGPoiss,AG,0,Max[\[CapitalNu]vals],subs];

lLCMPoiss[a_,h_,\[Gamma]_]:=PoislL[CM]
FisherCMPoiss=fisher[lLCMPoiss,{a, h,\[Gamma]}]//FullSimplify;

EFisherCMPoiss=EFisher[FisherCMPoiss,CM]//.subs;

DetCM=Det[EFisherCMPoiss];
DetCMTrunc=DetEFisherTrunc[FisherCMPoiss,CM,0,Max[\[CapitalNu]vals],subs];

lLBDPoiss[a_,h_,\[Gamma]_]:=PoislL[BD]
FisherBDPoiss=fisher[lLBDPoiss,{a, h,\[Gamma]}]//FullSimplify;

EFisherBDPoiss=EFisher[FisherBDPoiss,BD]//.subs;

DetBD=Det[EFisherBDPoiss];
DetBDTrunc=DetEFisherTrunc[FisherBDPoiss,BD,0,Max[\[CapitalNu]vals],subs];

lLAAPoiss[a_,h_,m_]:=PoislL[AA]
FisherAAPoiss=fisher[lLAAPoiss,{a, h,m}]//FullSimplify;

EFisherAAPoiss=EFisher[FisherAAPoiss,AA]//.subs;

DetAA=Det[EFisherAAPoiss];
DetAATrunc=DetEFisherTrunc[FisherAAPoiss,AA,0,Max[\[CapitalNu]vals],subs];

Save["FisherMatrix-1_Output_k1k2k3","`*"];


(* ::Input::Initialization:: *)
accgoal=3;
P\[CapitalNu]rats=Sort[DeleteDuplicates[Join@@Outer[Divide,Pvals,\[CapitalNu]vals]]];
e=10;
ParmRangeA={
{a,0,1000},
DeleteDuplicates[Flatten[{{h},Sort[{0,P\[CapitalNu]rats,1,10^e},Less]}]],
DeleteDuplicates[Flatten[{{\[Gamma]},Sort[{0,1/Pvals,1,10^e},Less]}]],
DeleteDuplicates[Flatten[{{m},Sort[{0,P\[CapitalNu]rats,1,5},Less]}]]
};

NIntH1=Log[NIntegrate[Sqrt[DetH1],
ParmRangeA[[1]],
AccuracyGoal->accgoal]]

NIntH1Trunc=Log[NIntegrate[Sqrt[DetH1Trunc],
ParmRangeA[[1]],
AccuracyGoal->accgoal]]

NIntR=Log[NIntegrate[Sqrt[DetR],
ParmRangeA[[1]],
AccuracyGoal->accgoal]]

NIntRTrunc=Log[NIntegrate[Sqrt[DetRTrunc],
ParmRangeA[[1]],
AccuracyGoal->accgoal]]

NIntHV=Log[NIntegrate[Sqrt[DetHV],
ParmRangeA[[1]],
ParmRangeA[[4]],
AccuracyGoal->accgoal]]

NIntHVTrunc=Log[NIntegrate[Sqrt[DetHVTrunc],
ParmRangeA[[1]],
ParmRangeA[[4]],
AccuracyGoal->accgoal,
Method->{"LocalAdaptive"},
MaxRecursion->500]]

NIntH2=Log[NIntegrate[Sqrt[DetH2],
ParmRangeA[[1]],
ParmRangeA[[2]],
AccuracyGoal->accgoal,
Method->{"LocalAdaptive"},
MaxRecursion->500]]

NIntH2Trunc=Log[NIntegrate[Sqrt[DetH2Trunc],
ParmRangeA[[1]],
ParmRangeA[[2]],
AccuracyGoal->accgoal,
Method->{"LocalAdaptive"},
MaxRecursion->500]]

NIntAG=Log[NIntegrate[Sqrt[DetAG],
ParmRangeA[[1]],
ParmRangeA[[2]],
AccuracyGoal->accgoal,
Method->{"LocalAdaptive"},
MaxRecursion->500]]

NIntAGTrunc=Log[NIntegrate[Sqrt[DetAGTrunc],
ParmRangeA[[1]],
ParmRangeA[[2]],
AccuracyGoal->accgoal,
Method->{"LocalAdaptive"},
MaxRecursion->500]]

NIntBD=Log[NIntegrate[Sqrt[DetBD],
ParmRangeA[[1]],
ParmRangeA[[2]],
ParmRangeA[[3]],
AccuracyGoal->accgoal,
Method->{"LocalAdaptive"},
MaxRecursion->500]]

NIntBDTrunc=Log[NIntegrate[Sqrt[DetBDTrunc],
ParmRangeA[[1]],
ParmRangeA[[2]],
ParmRangeA[[3]],
AccuracyGoal->accgoal,
Method->{"LocalAdaptive"},
MaxRecursion->500]]

NIntCM=Log[NIntegrate[Sqrt[DetCM],
ParmRangeA[[1]],
ParmRangeA[[2]],
ParmRangeA[[3]],
AccuracyGoal->accgoal,
Method->{"LocalAdaptive"},
MaxRecursion->500]]

NIntCMTrunc=Log[NIntegrate[Sqrt[DetCMTrunc],
ParmRangeA[[1]],
ParmRangeA[[2]],
ParmRangeA[[3]],
AccuracyGoal->accgoal,
Method->{"LocalAdaptive"},
MaxRecursion->500]]

NIntAA=Log[NIntegrate[Sqrt[DetAA],
ParmRangeA[[1]],
ParmRangeA[[2]],
ParmRangeA[[4]],
AccuracyGoal->accgoal,
MaxRecursion->500]]

NIntAATrunc=Log[NIntegrate[Sqrt[DetAATrunc],
ParmRangeA[[1]],
ParmRangeA[[2]],
ParmRangeA[[4]],
AccuracyGoal->accgoal,
Method->{"LocalAdaptive"},
MaxRecursion->500]]

Save["FisherMatrix-1_Output_k1k2k3_NInt","`*"];


(* ::Input::Initialization:: *)
ParmRangeB=Drop[ParmRangeA,1];

aNIntH1=Sqrt[DetH1];
aNIntH1Trunc=Sqrt[DetH1Trunc];

aNIntR=Sqrt[DetR];
aNIntRTrunc=Sqrt[DetRTrunc];

aNIntHV=NIntegrate[Sqrt[DetHV],
ParmRangeB[[3]],
AccuracyGoal->accgoal];
aNIntHVTrunc[\[Alpha]_?NumericQ]:=NIntegrate[Sqrt[DetHVTrunc/.{a->\[Alpha]}],
ParmRangeB[[3]],
AccuracyGoal->accgoal];

aNIntH2[\[Alpha]_?NumericQ]:=NIntegrate[Sqrt[DetH2/.{a->\[Alpha]}],
ParmRangeB[[1]],
AccuracyGoal->accgoal];
aNIntH2Trunc[\[Alpha]_?NumericQ]:=NIntegrate[Sqrt[DetH2Trunc/.{a->\[Alpha]}],
ParmRangeB[[1]],
AccuracyGoal->accgoal];

aNIntAG[\[Alpha]_?NumericQ]:=NIntegrate[Sqrt[DetAG/.{a->\[Alpha]}],
ParmRangeB[[1]],
AccuracyGoal->accgoal];
aNIntAGTrunc[\[Alpha]_?NumericQ]:=NIntegrate[Sqrt[DetAGTrunc/.{a->\[Alpha]}],
ParmRangeB[[1]],
AccuracyGoal->accgoal];

aNIntBD[\[Alpha]_?NumericQ]:=NIntegrate[Sqrt[DetBD/.{a->\[Alpha]}],
ParmRangeB[[1]],
ParmRangeB[[2]],
AccuracyGoal->accgoal];
aNIntBDTrunc[\[Alpha]_?NumericQ]:=NIntegrate[Sqrt[DetBDTrunc/.{a->\[Alpha]}],
ParmRangeB[[1]],
ParmRangeB[[2]],
AccuracyGoal->accgoal];

aNIntCM[\[Alpha]_?NumericQ]:=NIntegrate[Sqrt[DetCM/.{a->\[Alpha]}],
ParmRangeB[[1]],
ParmRangeB[[2]],
AccuracyGoal->accgoal];
aNIntCMTrunc[\[Alpha]_?NumericQ]:=NIntegrate[Sqrt[DetCMTrunc/.{a->\[Alpha]}],
ParmRangeB[[1]],
ParmRangeB[[2]],
AccuracyGoal->accgoal];

aNIntAA[\[Alpha]_?NumericQ]:=NIntegrate[Sqrt[DetAA/.{a->\[Alpha]}],
ParmRangeB[[1]],
ParmRangeB[[3]],
AccuracyGoal->accgoal];
aNIntAATrunc[\[Alpha]_?NumericQ]:=NIntegrate[Sqrt[DetAATrunc/.{a->\[Alpha]}],
ParmRangeB[[1]],
ParmRangeB[[3]],
AccuracyGoal->accgoal];

Save["FisherMatrix-1_Output_k1k2k3_aNInt","`*"];
