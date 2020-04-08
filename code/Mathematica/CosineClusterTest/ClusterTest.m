(* ::Package:: *)

(* ::Subsection:: *)
(*This is a test of using Mathematica on the HPC.*)


(* ::Input::Initialization:: *)
a=3
b=4


(* ::Input::Initialization:: *)
a+b


(* ::Input::Initialization:: *)
Export["ClusterTest.txt",{a,b,a+b}]


(* ::Input::Initialization:: *)
Save["ClusterTestResults","`*"];
