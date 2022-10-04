(* ::Package:: *)

(* ::Title:: *)
(*Analytic Evaluation of Three-Electron Integrals*)


(* ::Author:: *)
(*Hao Feng*)
(**)
(*Xihua University*)


(* ::Date:: *)
(*2016. 05. 15 - 2016. 07. 21*)


(* ::Text:: *)
(*This is a Mathematica notebook file. If you want to use the definitions of the functions in other notebook, you should save it as "Mathematica package file" first and get it.*)


(* ::Text:: *)
(*The key to save this file as package is to store the defintions of functions as "code" instead of "input".*)


(* ::Section:: *)
(*Introduction*)


(* ::Subsection:: *)
(*Closed Form Formulae*)


(* ::Text:: *)
(*This paper is to evaluate the following integrals*)


(* ::EquationNumbered:: *)
(*J(Subscript[n, 1],Subscript[n, 2],Subscript[n, 3],Subscript[n, 12],Subscript[n, 23],Subscript[n, 31];Subscript[\[Alpha], 1],Subscript[\[Alpha], 2],Subscript[\[Alpha], 3],Subscript[\[Alpha], 12],Subscript[\[Alpha], 23],Subscript[\[Alpha], 31])=\[Integral]Subscript[r, 1]^(Subscript[n, 1]-1) Subscript[r, 2]^(Subscript[n, 2]-1) Subscript[r, 3]^(Subscript[n, 3]-1) Subscript[r, 12]^(Subscript[n, 12]-1) Subscript[r, 23]^(Subscript[n, 23]-1) Subscript[r, 31]^(Subscript[n, 31]-1) E^(-Subscript[\[Alpha], 1] Subscript[r, 1]-Subscript[\[Alpha], 2] Subscript[r, 2]-Subscript[\[Alpha], 3] Subscript[r, 3]-Subscript[\[Alpha], 12] Subscript[r, 12]-Subscript[\[Alpha], 23] Subscript[r, 23]-Subscript[\[Alpha], 31] Subscript[r, 31]) \[DifferentialD]^3Subscript[r, 1]\[DifferentialD]^3Subscript[r, 2]\[DifferentialD]^3Subscript[r, 3]*)


(* ::Text:: *)
(*for all non-negative integers Subscript[n, 1], Subscript[n, 2], Subscript[n, 3], Subscript[n, 12], Subscript[n, 23], Subscript[n, 31], where Subscript[r, ij]=|Subscript[*)
(*\!\(\*OverscriptBox[\(r\), \(\[RightVector]\)]\), i]- Subscript[*)
(*\!\(\*OverscriptBox[\(r\), \(\[RightVector]\)]\), j]|. This is achieved by deriving an analytic formula for the generating integral*)


(* ::EquationNumbered:: *)
(*I(Subscript[\[Alpha], 1],Subscript[\[Alpha], 2],Subscript[\[Alpha], 3],Subscript[\[Alpha], 12],Subscript[\[Alpha], 23],Subscript[\[Alpha], 31])=\[Integral](Subscript[r, 1] Subscript[r, 2] Subscript[r, 3] Subscript[r, 12] Subscript[r, 23] Subscript[r, 31])^-1 E^(-Subscript[\[Alpha], 1] Subscript[r, 1]-Subscript[\[Alpha], 2] Subscript[r, 2]-Subscript[\[Alpha], 3] Subscript[r, 3]-Subscript[\[Alpha], 12] Subscript[r, 12]-Subscript[\[Alpha], 23] Subscript[r, 23]-Subscript[\[Alpha], 31] Subscript[r, 31]) \[DifferentialD]^3Subscript[r, 1]\[DifferentialD]^3Subscript[r, 2]\[DifferentialD]^3Subscript[r, 3]*)


(* ::Text:: *)
(*The integral J is obtained from I by differentiation,*)


(* ::EquationNumbered:: *)
(*J(Subscript[n, 1],Subscript[n, 2],Subscript[n, 3],Subscript[n, 12],Subscript[n, 23],Subscript[n, 31];Subscript[\[Alpha], 1],Subscript[\[Alpha], 2],Subscript[\[Alpha], 3],Subscript[\[Alpha], 12],Subscript[\[Alpha], 23],Subscript[\[Alpha], 31])=(-(\[PartialD]/\[PartialD]Subscript[\[Alpha], 1]))^Subscript[n, 1] (-(\[PartialD]/\[PartialD]Subscript[\[Alpha], 2]))^Subscript[n, 2] (-(\[PartialD]/\[PartialD]Subscript[\[Alpha], 3]))^Subscript[n, 3] (-(\[PartialD]/\[PartialD]Subscript[\[Alpha], 12]))^Subscript[n, 12] (-(\[PartialD]/\[PartialD]Subscript[\[Alpha], 23]))^Subscript[n, 23] (-(\[PartialD]/\[PartialD]Subscript[\[Alpha], 31]))^Subscript[n, 31] I(Subscript[\[Alpha], 1],Subscript[\[Alpha], 2],Subscript[\[Alpha], 3],Subscript[\[Alpha], 12],Subscript[\[Alpha], 23],Subscript[\[Alpha], 31])*)


(* ::Subsection:: *)
(*Power Series Formulae*)


(* ::Text:: *)
(*The power series expansion method is also used to calculate the three-electron integrals. *)


(* ::EquationNumbered:: *)
(*I=\[Integral]Subscript[r, 1]^K Subscript[r, 2]^L Subscript[r, 3]^M Subscript[r, 23]^l Subscript[r, 31]^m Subscript[r, 12]^n e^(-Subscript[a, 1] Subscript[r, 1]-Subscript[a, 2] Subscript[r, 2]-Subscript[a, 3] Subscript[r, 3]) \[DifferentialD]^3Subscript[r, 1]\[DifferentialD]^3Subscript[r, 2]\[DifferentialD]^3Subscript[r, 3]*)


(* ::Text:: *)
(*Notice that the definition of I and J is different:*)


(* ::Item1:: *)
(*I does not have Subscript[a, 12], Subscript[a, 23], and Subscript[a, 31];*)


(* ::Item1:: *)
(*the order of l, m, n is different from Subscript[n, 12], Subscript[n, 23], Subscript[n, 31]*)


(* ::Text:: *)
(*If any of l, m, n is even, the series is finite. Otherwise the series is infinite and one could use some accelerator approach to improve the accuracy and computing speed. *)


(* ::Text:: *)
(*e.g., Levin u-transformation (J. S. Sims and S. A. Hagstrom, J. Phys. B, 37, 1519 (2004))*)


(* ::Section:: *)
(*Basic Integrals*)


(* ::Subsection:: *)
(*auxilary functions*)


(* F. E. Harris, Phys. Rev. A, 55, 1820 (1997) *)

Clear[sign,abs,hu,hv,hs,hmu,\[Sigma]2];

sign[x_]:=If[x<0,-1,1];
abs[x_]:=If[x<0,-x,x];

(* u and v functions *)
hu[z_]:=PolyLog[2,z]-PolyLog[2,1/z]; 

hv[z_]:=1/2 PolyLog[2,(1-z)/2]-1/4 (Log[(1-z)/2])^2-1/2 PolyLog[2,(1+z)/2]+1/4 ((Log[(1+z)/2]))^2;

(* \[Sigma] function; In current cases, a12, a23 and a31 will always be set to zero so \[Sigma] function is single-valued. 
   It means \[Sigma] is always real ! 
*)
hs[a1_,a2_,a3_,a12_,a23_,a31_]:= Sqrt[a1^2 a2^2 a3^2+a1^2 a23^2 (a1^2-a2^2-a3^2-a12^2+a23^2-a31^2)
                                  +a1^2 a12^2 a31^2+a2^2 a31^2 (-a1^2+a2^2-a3^2-a12^2-a23^2+a31^2)
                                  +a2^2 a23^2 a12^2+a3^2 a12^2 (-a1^2-a2^2+a3^2+a12^2-a23^2-a31^2)
                                  +a3^2 a31^2 a23^2];
\[Sigma]2[a1_,a2_,a3_,a12_,a23_,a31_]:= a1^2 a2^2 a3^2+a1^2 a23^2 (a1^2-a2^2-a3^2-a12^2+a23^2-a31^2)+
                                 a1^2 a12^2 a31^2+a2^2 a31^2 (-a1^2+a2^2-a3^2-a12^2-a23^2+a31^2)+
                                 a2^2 a23^2 a12^2+a3^2 a12^2 (-a1^2-a2^2+a3^2+a12^2-a23^2-a31^2)+
                                 a3^2 a31^2 a23^2;

(* \[Mu] functions; subscript is first; then superscript;*)
hmu[0,0,a1_,a2_,a3_,a12_,a23_,a31_]:=2a12 a23 a31;
hmu[0,1,a1_,a2_,a3_,a12_,a23_,a31_]:=a23(a2^2+a3^2-a1^2);
hmu[0,2,a1_,a2_,a3_,a12_,a23_,a31_]:=a31(a1^2+a3^2-a2^2);
hmu[0,3,a1_,a2_,a3_,a12_,a23_,a31_]:=a12(a1^2+a2^2-a3^2);

hmu[1,1,a1_,a2_,a3_,a12_,a23_,a31_]:=2a2 a23 a3;
hmu[1,0,a1_,a2_,a3_,a12_,a23_,a31_]:=a23(a12^2+a31^2-a1^2);
hmu[1,2,a1_,a2_,a3_,a12_,a23_,a31_]:=a3(a1^2+a31^2-a12^2);
hmu[1,3,a1_,a2_,a3_,a12_,a23_,a31_]:=a2(a1^2+a12^2-a31^2);

hmu[2,2,a1_,a2_,a3_,a12_,a23_,a31_]:=2a1 a31 a3;
hmu[2,0,a1_,a2_,a3_,a12_,a23_,a31_]:=a31(a12^2+a23^2-a2^2);
hmu[2,1,a1_,a2_,a3_,a12_,a23_,a31_]:=a3(a2^2+a23^2-a12^2);
hmu[2,3,a1_,a2_,a3_,a12_,a23_,a31_]:=a1(a2^2+a12^2-a23^2);

hmu[3,3,a1_,a2_,a3_,a12_,a23_,a31_]:=2a1 a12 a2;
hmu[3,0,a1_,a2_,a3_,a12_,a23_,a31_]:=a12(a31^2+a23^2-a3^2);
hmu[3,1,a1_,a2_,a3_,a12_,a23_,a31_]:=a2(a3^2+a23^2-a31^2);
hmu[3,2,a1_,a2_,a3_,a12_,a23_,a31_]:=a1(a3^2+a31^2-a23^2);

Clear[xfactor,hpt,hp];

(* p functions; To judge if \[Sigma] is equal to \[CapitalGamma] or \[Gamma] *)
xfactor = 0.1; (* small factor to judge how close z is from \[PlusMinus]1 *)

hpt[0,0,a1_,a2_,a3_,a12_,a23_,a31_]:=a1+a2+a3;
hpt[0,1,a1_,a2_,a3_,a12_,a23_,a31_]:=-a1+a2+a3;
hpt[0,2,a1_,a2_,a3_,a12_,a23_,a31_]:=a1-a2+a3;
hpt[0,3,a1_,a2_,a3_,a12_,a23_,a31_]:=a1+a2-a3;

hpt[1,1,a1_,a2_,a3_,a12_,a23_,a31_]:=a1+a12+a31;
hpt[1,0,a1_,a2_,a3_,a12_,a23_,a31_]:=-a1+a12+a31;
hpt[1,2,a1_,a2_,a3_,a12_,a23_,a31_]:=a1-a12+a31;
hpt[1,3,a1_,a2_,a3_,a12_,a23_,a31_]:=a1+a12-a31;

hpt[2,2,a1_,a2_,a3_,a12_,a23_,a31_]:=a12+a2+a23;
hpt[2,0,a1_,a2_,a3_,a12_,a23_,a31_]:=a12-a2+a23;
hpt[2,1,a1_,a2_,a3_,a12_,a23_,a31_]:=-a12+a2+a23;
hpt[2,3,a1_,a2_,a3_,a12_,a23_,a31_]:=a12+a2-a23;

hpt[3,3,a1_,a2_,a3_,a12_,a23_,a31_]:=a31+a23+a3;
hpt[3,0,a1_,a2_,a3_,a12_,a23_,a31_]:=a31+a23-a3;
hpt[3,1,a1_,a2_,a3_,a12_,a23_,a31_]:=-a31+a23+a3;
hpt[3,2,a1_,a2_,a3_,a12_,a23_,a31_]:=a31-a23+a3;

hp[0,0,a1_,a2_,a3_,a12_,a23_,a31_]:=a1+a2+a3;
hp[0,1,a1_,a2_,a3_,a12_,a23_,a31_]:=Module[{t},t=-a1+a2+a3;If[abs[t]<xfactor,1,-a1+a2+a3]];
hp[0,2,a1_,a2_,a3_,a12_,a23_,a31_]:=Module[{t},t=a1-a2+a3;If[abs[t]<xfactor,1,a1-a2+a3]];
hp[0,3,a1_,a2_,a3_,a12_,a23_,a31_]:=Module[{t},t=a1+a2-a3;If[abs[t]<xfactor,1,a1+a2-a3]];

hp[1,1,a1_,a2_,a3_,a12_,a23_,a31_]:=a1+a12+a31;
hp[1,0,a1_,a2_,a3_,a12_,a23_,a31_]:=Module[{t},t=-a1+a12+a31;If[abs[t]<xfactor,1,-a1+a12+a31]];
hp[1,2,a1_,a2_,a3_,a12_,a23_,a31_]:=Module[{t},t=a1-a12+a31;If[abs[t]<xfactor,1,a1-a12+a31]];
hp[1,3,a1_,a2_,a3_,a12_,a23_,a31_]:=Module[{t},t=a1+a12-a31;If[abs[t]<xfactor,1,a1+a12-a31]];

hp[2,2,a1_,a2_,a3_,a12_,a23_,a31_]:=a12+a2+a23;
hp[2,0,a1_,a2_,a3_,a12_,a23_,a31_]:=Module[{t},t=a12-a2+a23;If[abs[t]<xfactor,1,a12-a2+a23]];
hp[2,1,a1_,a2_,a3_,a12_,a23_,a31_]:=Module[{t},t=-a12+a2+a23;If[abs[t]<xfactor,1,-a12+a2+a23]];
hp[2,3,a1_,a2_,a3_,a12_,a23_,a31_]:=Module[{t},t=a12+a2-a23;If[abs[t]<xfactor,1,a12+a2-a23]];

hp[3,3,a1_,a2_,a3_,a12_,a23_,a31_]:=a31+a23+a3;
hp[3,0,a1_,a2_,a3_,a12_,a23_,a31_]:=Module[{t},t=a31+a23-a3;If[abs[t]<xfactor,1,a31+a23-a3]];
hp[3,1,a1_,a2_,a3_,a12_,a23_,a31_]:=Module[{t},t=-a31+a23+a3;If[abs[t]<xfactor,1,-a31+a23+a3]];
hp[3,2,a1_,a2_,a3_,a12_,a23_,a31_]:=Module[{t},t=a31-a23+a3;If[abs[t]<xfactor,1,a31-a23+a3]];

Clear[hg,hb,hq,hG];

(* \[Gamma] functions; subscript is first; then superscript;*)
hg[k_,j_,a1_,a2_,a3_,a12_,a23_,a31_]:=If[k==j,Sum[hmu[i,j,a1,a2,a3,a12,a23,a31],{i,0,3}], 
                                             hg[j,j,a1,a2,a3,a12,a23,a31]-2(hmu[j,j,a1,a2,a3,a12,a23,a31]+
                                                                            hmu[k,j,a1,a2,a3,a12,a23,a31])];

(* \[Beta] functions; subscript is first; then superscript;*)
hb[k_,j_,a1_,a2_,a3_,a12_,a23_,a31_]:=(hs[a1,a2,a3,a12,a23,a31]-hg[k,j,a1,a2,a3,a12,a23,a31])/(hs[a1,a2,a3,a12,a23,a31]+hg[k,j,a1,a2,a3,a12,a23,a31]);

(* \[CapitalGamma] functions;*)
hq[1,a1_,a2_,a3_,a12_,a23_,a31_]:=a2+a3+a12+a31;
hq[2,a1_,a2_,a3_,a12_,a23_,a31_]:=a1+a3+a12+a23;
hq[3,a1_,a2_,a3_,a12_,a23_,a31_]:=a1+a2+a31+a23;

hG[1,a1_,a2_,a3_,a12_,a23_,a31_]:=((a1^2+(a2+a3)(a12+a31))(a23^2+(a2+a12)(a3+a31)))/hq[1,a1,a2,a3,a12,a23,a31]-hq[1,a1,a2,a3,a12,a23,a31](a2 a31+a3 a12);
hG[2,a1_,a2_,a3_,a12_,a23_,a31_]:=((a2^2+(a1+a3)(a12+a23))(a31^2+(a1+a12)(a3+a23)))/hq[2,a1,a2,a3,a12,a23,a31]-hq[2,a1,a2,a3,a12,a23,a31](a1 a23+a3 a12);
hG[3,a1_,a2_,a3_,a12_,a23_,a31_]:=((a3^2+(a1+a2)(a31+a23))(a12^2+(a1+a31)(a2+a23)))/hq[3,a1,a2,a3,a12,a23,a31]-hq[3,a1,a2,a3,a12,a23,a31](a1 a23+a2 a31);

(* Abbreviations of Eq.[44] *)
Clear[hld];

hld[0,0,a1_,a2_,a3_,a12_,a23_,a31_]:= hp[1,1,a1,a2,a3,a12,a23,a31]hp[2,2,a1,a2,a3,a12,a23,a31]hp[3,3,a1,a2,a3,a12,a23,a31]*
                                      hp[1,0,a1,a2,a3,a12,a23,a31]hp[2,0,a1,a2,a3,a12,a23,a31]hp[3,0,a1,a2,a3,a12,a23,a31]*
                                      1/(abs[hg[0,0,a1,a2,a3,a12,a23,a31]]+hs[a1,a2,a3,a12,a23,a31])^2;

hld[0,1,a1_,a2_,a3_,a12_,a23_,a31_]:= hp[1,1,a1,a2,a3,a12,a23,a31]hp[1,0,a1,a2,a3,a12,a23,a31]hp[2,1,a1,a2,a3,a12,a23,a31]*
                                      hp[3,1,a1,a2,a3,a12,a23,a31]hp[2,3,a1,a2,a3,a12,a23,a31]hp[3,2,a1,a2,a3,a12,a23,a31]*
                                      1/(abs[hg[1,0,a1,a2,a3,a12,a23,a31]]+hs[a1,a2,a3,a12,a23,a31])^2;
hld[0,2,a1_,a2_,a3_,a12_,a23_,a31_]:= hp[2,2,a1,a2,a3,a12,a23,a31]hp[2,0,a1,a2,a3,a12,a23,a31]hp[1,2,a1,a2,a3,a12,a23,a31]*
                                      hp[3,2,a1,a2,a3,a12,a23,a31]hp[1,3,a1,a2,a3,a12,a23,a31]hp[3,1,a1,a2,a3,a12,a23,a31]*
                                      1/(abs[hg[2,0,a1,a2,a3,a12,a23,a31]]+hs[a1,a2,a3,a12,a23,a31])^2;
hld[0,3,a1_,a2_,a3_,a12_,a23_,a31_]:= hp[3,3,a1,a2,a3,a12,a23,a31]hp[3,0,a1,a2,a3,a12,a23,a31]hp[1,3,a1,a2,a3,a12,a23,a31]*
                                      hp[2,3,a1,a2,a3,a12,a23,a31]hp[1,2,a1,a2,a3,a12,a23,a31]hp[2,1,a1,a2,a3,a12,a23,a31]*
                                      1/(abs[hg[3,0,a1,a2,a3,a12,a23,a31]]+hs[a1,a2,a3,a12,a23,a31])^2;

hld[1,1,a1_,a2_,a3_,a12_,a23_,a31_]:= hp[0,0,a1,a2,a3,a12,a23,a31]hp[2,2,a1,a2,a3,a12,a23,a31]hp[3,3,a1,a2,a3,a12,a23,a31]*
                                      hp[0,1,a1,a2,a3,a12,a23,a31]hp[2,1,a1,a2,a3,a12,a23,a31]hp[3,1,a1,a2,a3,a12,a23,a31]*
                                      1/(abs[hg[1,1,a1,a2,a3,a12,a23,a31]]+hs[a1,a2,a3,a12,a23,a31])^2;
hld[1,0,a1_,a2_,a3_,a12_,a23_,a31_]:= hp[0,0,a1,a2,a3,a12,a23,a31]hp[0,1,a1,a2,a3,a12,a23,a31]hp[2,0,a1,a2,a3,a12,a23,a31]*
                                      hp[3,0,a1,a2,a3,a12,a23,a31]hp[2,3,a1,a2,a3,a12,a23,a31]hp[3,2,a1,a2,a3,a12,a23,a31]*
                                      1/(abs[hg[0,1,a1,a2,a3,a12,a23,a31]]+hs[a1,a2,a3,a12,a23,a31])^2;
hld[1,2,a1_,a2_,a3_,a12_,a23_,a31_]:= hp[2,2,a1,a2,a3,a12,a23,a31]hp[2,1,a1,a2,a3,a12,a23,a31]hp[0,2,a1,a2,a3,a12,a23,a31]*
                                      hp[3,2,a1,a2,a3,a12,a23,a31]hp[0,3,a1,a2,a3,a12,a23,a31]hp[3,0,a1,a2,a3,a12,a23,a31]*
                                      1/(abs[hg[2,1,a1,a2,a3,a12,a23,a31]]+hs[a1,a2,a3,a12,a23,a31])^2;
hld[1,3,a1_,a2_,a3_,a12_,a23_,a31_]:= hp[3,3,a1,a2,a3,a12,a23,a31]hp[3,1,a1,a2,a3,a12,a23,a31]hp[0,3,a1,a2,a3,a12,a23,a31]*
                                      hp[2,3,a1,a2,a3,a12,a23,a31]hp[0,2,a1,a2,a3,a12,a23,a31]hp[2,0,a1,a2,a3,a12,a23,a31]*
                                      1/(abs[hg[3,1,a1,a2,a3,a12,a23,a31]]+hs[a1,a2,a3,a12,a23,a31])^2;

hld[2,2,a1_,a2_,a3_,a12_,a23_,a31_]:= hp[1,1,a1,a2,a3,a12,a23,a31]hp[0,0,a1,a2,a3,a12,a23,a31]hp[3,3,a1,a2,a3,a12,a23,a31]*
                                      hp[1,2,a1,a2,a3,a12,a23,a31]hp[0,2,a1,a2,a3,a12,a23,a31]hp[3,2,a1,a2,a3,a12,a23,a31]*
                                      1/(abs[hg[2,2,a1,a2,a3,a12,a23,a31]]+hs[a1,a2,a3,a12,a23,a31])^2;
hld[2,1,a1_,a2_,a3_,a12_,a23_,a31_]:= hp[1,1,a1,a2,a3,a12,a23,a31]hp[1,2,a1,a2,a3,a12,a23,a31]hp[0,1,a1,a2,a3,a12,a23,a31]*
                                      hp[3,1,a1,a2,a3,a12,a23,a31]hp[0,3,a1,a2,a3,a12,a23,a31]hp[3,0,a1,a2,a3,a12,a23,a31]*
                                      1/(abs[hg[1,2,a1,a2,a3,a12,a23,a31]]+hs[a1,a2,a3,a12,a23,a31])^2;
hld[2,0,a1_,a2_,a3_,a12_,a23_,a31_]:= hp[0,0,a1,a2,a3,a12,a23,a31]hp[0,2,a1,a2,a3,a12,a23,a31]hp[1,0,a1,a2,a3,a12,a23,a31]*
                                      hp[3,0,a1,a2,a3,a12,a23,a31]hp[1,3,a1,a2,a3,a12,a23,a31]hp[3,1,a1,a2,a3,a12,a23,a31]*
                                      1/(abs[hg[0,2,a1,a2,a3,a12,a23,a31]]+hs[a1,a2,a3,a12,a23,a31])^2;
hld[2,3,a1_,a2_,a3_,a12_,a23_,a31_]:= hp[3,3,a1,a2,a3,a12,a23,a31]hp[3,2,a1,a2,a3,a12,a23,a31]hp[1,3,a1,a2,a3,a12,a23,a31]*
                                      hp[0,3,a1,a2,a3,a12,a23,a31]hp[1,0,a1,a2,a3,a12,a23,a31]hp[0,1,a1,a2,a3,a12,a23,a31]*
                                      1/(abs[hg[3,2,a1,a2,a3,a12,a23,a31]]+hs[a1,a2,a3,a12,a23,a31])^2;

hld[3,3,a1_,a2_,a3_,a12_,a23_,a31_]:= hp[1,1,a1,a2,a3,a12,a23,a31]hp[2,2,a1,a2,a3,a12,a23,a31]hp[0,0,a1,a2,a3,a12,a23,a31]*
                                      hp[1,3,a1,a2,a3,a12,a23,a31]hp[2,3,a1,a2,a3,a12,a23,a31]hp[0,3,a1,a2,a3,a12,a23,a31]*
                                      1/(abs[hg[3,3,a1,a2,a3,a12,a23,a31]]+hs[a1,a2,a3,a12,a23,a31])^2;
hld[3,1,a1_,a2_,a3_,a12_,a23_,a31_]:= hp[1,1,a1,a2,a3,a12,a23,a31]hp[1,3,a1,a2,a3,a12,a23,a31]hp[2,1,a1,a2,a3,a12,a23,a31]*
                                      hp[0,1,a1,a2,a3,a12,a23,a31]hp[2,0,a1,a2,a3,a12,a23,a31]hp[0,2,a1,a2,a3,a12,a23,a31]*
                                      1/(abs[hg[1,3,a1,a2,a3,a12,a23,a31]]+hs[a1,a2,a3,a12,a23,a31])^2;
hld[3,2,a1_,a2_,a3_,a12_,a23_,a31_]:= hp[2,2,a1,a2,a3,a12,a23,a31]hp[2,3,a1,a2,a3,a12,a23,a31]hp[1,2,a1,a2,a3,a12,a23,a31]*
                                      hp[0,2,a1,a2,a3,a12,a23,a31]hp[1,0,a1,a2,a3,a12,a23,a31]hp[0,1,a1,a2,a3,a12,a23,a31]*
                                      1/(abs[hg[2,3,a1,a2,a3,a12,a23,a31]]+hs[a1,a2,a3,a12,a23,a31])^2;
hld[3,0,a1_,a2_,a3_,a12_,a23_,a31_]:= hp[0,0,a1,a2,a3,a12,a23,a31]hp[0,3,a1,a2,a3,a12,a23,a31]hp[1,0,a1,a2,a3,a12,a23,a31]*
                                      hp[2,0,a1,a2,a3,a12,a23,a31]hp[1,2,a1,a2,a3,a12,a23,a31]hp[2,1,a1,a2,a3,a12,a23,a31]*
                                      1/(abs[hg[0,3,a1,a2,a3,a12,a23,a31]]+hs[a1,a2,a3,a12,a23,a31])^2;

(* Abbreviations of Eq.[46] *)
Clear[hdd,hcd];

hdd[a1_,a2_,a3_,a12_,a23_,a31_]:=Product[hp[j,j,a1,a2,a3,a12,a23,a31],{j,0,3}];

hcd[1,a1_,a2_,a3_,a12_,a23_,a31_]:=hdd[a1,a2,a3,a12,a23,a31]hp[0,1,a1,a2,a3,a12,a23,a31]hp[1,0,a1,a2,a3,a12,a23,a31]*
                                                            hp[3,2,a1,a2,a3,a12,a23,a31]hp[2,3,a1,a2,a3,a12,a23,a31]*
                                                            1/(hq[1,a1,a2,a3,a12,a23,a31](abs[hG[1,a1,a2,a3,a12,a23,a31]]+hs[a1,a2,a3,a12,a23,a31]))^2;
hcd[2,a1_,a2_,a3_,a12_,a23_,a31_]:=hdd[a1,a2,a3,a12,a23,a31]hp[0,2,a1,a2,a3,a12,a23,a31]hp[2,0,a1,a2,a3,a12,a23,a31]*
                                                            hp[1,3,a1,a2,a3,a12,a23,a31]hp[3,1,a1,a2,a3,a12,a23,a31]*
                                                            1/(hq[2,a1,a2,a3,a12,a23,a31](abs[hG[2,a1,a2,a3,a12,a23,a31]]+hs[a1,a2,a3,a12,a23,a31]))^2;
hcd[3,a1_,a2_,a3_,a12_,a23_,a31_]:=hdd[a1,a2,a3,a12,a23,a31]hp[0,3,a1,a2,a3,a12,a23,a31]hp[3,0,a1,a2,a3,a12,a23,a31]*
                                                            hp[1,2,a1,a2,a3,a12,a23,a31]hp[2,1,a1,a2,a3,a12,a23,a31]*
                                                            1/(hq[3,a1,a2,a3,a12,a23,a31](abs[hG[3,a1,a2,a3,a12,a23,a31]]+hs[a1,a2,a3,a12,a23,a31]))^2;

Clear[hV,hVN,hVN1,hDL,hdlseries,hAD];

(* V function *)                                       
hV[z_]:=If[Abs[z]<1,hv[z],
                    hv[1/z]+PolyLog[2,1-z]+PolyLog[2,-z]+Log[z]*Log[1+z]+(\[Pi]^2/12)];

hVN[d_,g_]:=sign[g](-Log[abs[d]]^2/4-\[Pi]^2/12+PolyLog[2,(1-abs[g])/2]+Log[(1+abs[g])/2]^2/2);

hVN1[z_]:=sign[z](-Log[abs[(1-abs[z])/(1+abs[z])]]^2/4-\[Pi]^2/12+PolyLog[2,(1-abs[z])/2]+Log[(1+abs[z])/2]^2/2);


hDL[x_]:=Module[{w},Which[abs[x]<=1/2,hdlseries[x],
                          x<-1,     w=1/(1-x);hdlseries[w]+Log[w] Log[w x^2]/2-\[Pi]^2/6,
                          x<-(1/2),    w=x/(x-1);-hdlseries[w]-Log[-w/x]^2/2,
                          x==1.,    \[Pi]^2/6,
                          x<1,      w=1-x;-hdlseries[w]-Log[x]Log[w]+\[Pi]^2/6,
                          x<=2,      w=1-1/x;hdlseries[w]-Log[x]Log[w]-Log[x]^2/2+\[Pi]^2/6,
                          True,     -hdlseries[1/x]-Log[x]^2/2+\[Pi]^2/3]];

hdlseries[x_]:=Sum[x^n/n^2,{n,1,90}]; 

hAD[j_,i1_,i2_,i3_,a1_,a2_,a3_,a12_,a23_,a31_]:=Module[{cj},If[abs[hp[j,i1,a1,a2,a3,a12,a23,a31]]>10^-16,
    cj=Min[i1+j,6-i1-j];Sign[hG[cj,a1,a2,a3,a12,a23,a31]]Log[abs[hp[j,i1,a1,a2,a3,a12,a23,a31]]]*
                        Log[hq[cj,a1,a2,a3,a12,a23,a31] (abs[hG[cj,a1,a2,a3,a12,a23,a31]]+hs[a1,a2,a3,a12,a23,a31])/(hq[cj,a1,a2,a3,a12,a23,a31]-hp[j,i1,a1,a2,a3,a12,a23,a31])]^2*
                        (abs[hg[i2,i1,a1,a2,a3,a12,a23,a31]]+hs[a1,a2,a3,a12,a23,a31])*
                        (abs[hg[i3,i1,a1,a2,a3,a12,a23,a31]]+hs[a1,a2,a3,a12,a23,a31])*
                         1/(abs[hg[i1,j,a1,a2,a3,a12,a23,a31]]+hs[a1,a2,a3,a12,a23,a31])*
                         1/(abs[hg[i1,i1,a1,a2,a3,a12,a23,a31]]+hs[a1,a2,a3,a12,a23,a31])*
                         1/(abs[hg[i3,i2,a1,a2,a3,a12,a23,a31]]+hs[a1,a2,a3,a12,a23,a31])*
                         1/(abs[hg[i2,i3,a1,a2,a3,a12,a23,a31]]+hs[a1,a2,a3,a12,a23,a31]),
   0]
];


Clear[hc,hD,hF];

(* \[Chi] function *)
hc[z_]:=1/2 PolyLog[2,z]-1/2 PolyLog[2,-z];

(* D function *)
hD[a1_,a2_,a3_,a12_,a23_,a31_]:= (8a12 a23 a31)/(hg[0,0,a1,a2,a3,a12,a23,a31]hg[1,0,a1,a2,a3,a12,a23,a31]hg[2,0,a1,a2,a3,a12,a23,a31]hg[3,0,a1,a2,a3,a12,a23,a31])+
                                 (8a2 a3 a23)/(hg[0,1,a1,a2,a3,a12,a23,a31]hg[1,1,a1,a2,a3,a12,a23,a31]hg[2,1,a1,a2,a3,a12,a23,a31]hg[3,1,a1,a2,a3,a12,a23,a31])-
  (1/(hg[3,3,a1,a2,a3,a12,a23,a31]hg[3,2,a1,a2,a3,a12,a23,a31])+1/(hg[2,2,a1,a2,a3,a12,a23,a31]hg[2,3,a1,a2,a3,a12,a23,a31])) 1/hG[1,a1,a2,a3,a12,a23,a31]-

  (1/(hg[0,0,a1,a2,a3,a12,a23,a31]hg[0,2,a1,a2,a3,a12,a23,a31])+1/(hg[1,1,a1,a2,a3,a12,a23,a31]hg[1,3,a1,a2,a3,a12,a23,a31])) 1/hG[2,a1,a2,a3,a12,a23,a31]-
  (1/(hg[1,1,a1,a2,a3,a12,a23,a31]hg[1,2,a1,a2,a3,a12,a23,a31])+1/(hg[0,0,a1,a2,a3,a12,a23,a31]hg[0,3,a1,a2,a3,a12,a23,a31])) 1/hG[3,a1,a2,a3,a12,a23,a31];

(* F function *)
hF[t_,a1_,a2_,a3_,a12_,a23_,a31_]:=2/t Log[Abs[t]]+
                                  1/hs[a1,a2,a3,a12,a23,a31] (hvs[hs[a1,a2,a3,a12,a23,a31]/t]+
                                                            2hchix[hs[a1,a2,a3,a12,a23,a31]/t]+
                                                            Log[Abs[hs[a1,a2,a3,a12,a23,a31]/t]]hls[hs[a1,a2,a3,a12,a23,a31]/t]);

Clear[nmax,hvs,hchix,hls];

nmax=10;
(* Overscript[v, _] function *)
hvs[z_]:=Module[{t,n},Sum[SeriesCoefficient[hv[t],{t,0,2n+1}]z^(2n+1),{n,1,nmax}]];
(* \[Chi] function *)
hchix[z_]:=Sum[z^(2n+1)/(2n+1)^2,{n,1,nmax}];
(* Overscript[l, _] function *)
hls[z_]:=-2Sum[z^(2n+1)/(2n+1),{n,1,nmax}];


(* ::Section:: *)
(*Working Formulae*)


(* ::Subsection:: *)
(*3-e Integrals*)


(* D. M. Fromm and R. N. Hill, Phys. Rev. A, 36, 1013 (1987) *)
Clear[Ifh];
Ifh[-1,-1,-1,-1,-1,-1,a1_,a2_,a3_,a12_,a23_,a31_]:=
	Ifh[-1,-1,-1,-1,-1,-1,a1,a2,a3,a12,a23,a31]=
	 (16\[Pi]^3)/hs[a1,a2,a3,a12,a23,a31] (Sum[hu[hb[0,0,a1,a2,a3,a12,a23,a31]hb[0,j,a1,a2,a3,a12,a23,a31]],{j,1,3}]
                              +Sum[hv[hg[k,j,a1,a2,a3,a12,a23,a31]/hs[a1,a2,a3,a12,a23,a31]],{j,0,3},{k,0,3}]);


(* E. Remiddi, Phys. Rev. A, 44, 5492 (1991) *)
Clear[Irem];

Irem[\[Alpha]1_,\[Alpha]2_,\[Alpha]3_]:=I3e[\[Alpha]1,\[Alpha]2,\[Alpha]3]=Module[{\[Zeta]},
	\[Zeta][1]=\[Alpha]1/(\[Alpha]2+\[Alpha]3);\[Zeta][2]=\[Alpha]2/(\[Alpha]1+\[Alpha]3);\[Zeta][3]=\[Alpha]3/(\[Alpha]1+\[Alpha]2);

	-((32\[Pi]^3)/(\[Alpha]1 \[Alpha]2 \[Alpha]3))\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(j = 1\), \(3\)]\((Log[\[Zeta][j]] Log[1 + \[Zeta][j]] + PolyLog[2, \(-\[Zeta][j]\)] + PolyLog[2, 1 - \[Zeta][j]])\)\)
];


(* F. E. Harris, Phys. Rev. A, 55, 1820 (1997) *)
Clear[hI0,hI1,hI2,hI3,hI4,hI5,hI];

hI0[-1,-1,-1,-1,-1,-1,a1_,a2_,a3_,a12_,a23_,a31_]:=
	hI0[-1,-1,-1,-1,-1,-1,a1,a2,a3,a12,a23,a31]=
    (16\[Pi]^3)/hs[a1,a2,a3,a12,a23,a31] (-2Sum[hv[hG[j,a1,a2,a3,a12,a23,a31]/hs[a1,a2,a3,a12,a23,a31]],{j,1,3}]
                               +Sum[hv[hg[k,j,a1,a2,a3,a12,a23,a31]/hs[a1,a2,a3,a12,a23,a31]],{j,0,3},{k,0,3}]);
                                                                      
hI1[-1,-1,-1,-1,-1,-1,a1_,a2_,a3_,a12_,a23_,a31_]:=
	hI1[-1,-1,-1,-1,-1,-1,a1,a2,a3,a12,a23,a31]=
    (16\[Pi]^3)/hs[a1,a2,a3,a12,a23,a31] (\[Pi]^2/2-2Sum[hV[hG[j,a1,a2,a3,a12,a23,a31]/hs[a1,a2,a3,a12,a23,a31]],{j,1,3}]
                                +Sum[hV[hg[k,j,a1,a2,a3,a12,a23,a31]/hs[a1,a2,a3,a12,a23,a31]],{j,0,3},{k,0,3}]);

hI2[-1,-1,-1,-1,-1,-1,a1_,a2_,a3_,a12_,a23_,a31_]:=
	hI2[-1,-1,-1,-1,-1,-1,a1,a2,a3,a12,a23,a31]=
	(16\[Pi]^3)/hs[a1,a2,a3,a12,a23,a31] (\[Pi]^2/2-2Sum[hVN[hcd[j,a1,a2,a3,a12,a23,a31],hG[j,a1,a2,a3,a12,a23,a31]/hs[a1,a2,a3,a12,a23,a31]],{j,1,3}]+
                                  Sum[hVN[hld[j,i,a1,a2,a3,a12,a23,a31],hg[i,j,a1,a2,a3,a12,a23,a31]/hs[a1,a2,a3,a12,a23,a31]],{j,0,3},{i,0,3}]);

hI3[-1,-1,-1,-1,-1,-1,a1_,a2_,a3_,a12_,a23_,a31_]:=
	hI3[-1,-1,-1,-1,-1,-1,a1,a2,a3,a12,a23,a31]=
    16\[Pi]^3 (2hD[a1,a2,a3,a12,a23,a31]hs[a1,a2,a3,a12,a23,a31]hs[a1,a2,a3,a12,a23,a31](Log[Abs[2hs[a1,a2,a3,a12,a23,a31]]]-1)
        -2Sum[hF[hG[j,a1,a2,a3,a12,a23,a31],a1,a2,a3,a12,a23,a31],{j,1,3}]
         +Sum[hF[hg[k,j,a1,a2,a3,a12,a23,a31],a1,a2,a3,a12,a23,a31],{j,0,3},{k,0,3}]);


hI4[-1,-1,-1,-1,-1,-1,a1_,a2_,a3_,a12_,a23_,a31_]:=
	hI4[-1,-1,-1,-1,-1,-1,a1,a2,a3,a12,a23,a31]=
    16\[Pi]^3 (-Sum[4/hG[j,a1,a2,a3,a12,a23,a31] Log[Abs[hG[j,a1,a2,a3,a12,a23,a31]]],{j,1,3}]
         +Sum[2/hg[k,j,a1,a2,a3,a12,a23,a31] Log[Abs[hg[k,j,a1,a2,a3,a12,a23,a31]]],{k,0,3},{j,0,3}]);

hI5[-1,-1,-1,-1,-1,-1,a1_,a2_,a3_,a12_,a23_,a31_]:=
	hI5[-1,-1,-1,-1,-1,-1,a1,a2,a3,a12,a23,a31]=
    (16\[Pi]^3)/hs[a1,a2,a3,a12,a23,a31] (\[Pi]^2/2-2Sum[hVN1[hG[j,a1,a2,a3,a12,a23,a31]/hs[a1,a2,a3,a12,a23,a31]],{j,1,3}]
                                +Sum[hVN1[hg[k,j,a1,a2,a3,a12,a23,a31]/hs[a1,a2,a3,a12,a23,a31]],{j,0,3},{k,0,3}]);

hI[K_,L_,M_,n1_,n2_,n3_,a1_,a2_,a3_,a12_,a23_,a31_]:=
	hI[K,L,M,n1,n2,n3,a1,a2,a3,a12,a23,a31]=
	Module[
	   {s=\[Sigma]2[a1,a2,a3,a12,a23,a31],
		Y1=0.0125, Y2=-0.002,
		x,y,z,a,b,c,int},

		Which[
			s>Y1,
                int=D[hI2[-1,-1,-1,-1,-1,-1,x,y,z,a,b,c],{x,K+1},{y,L+1},{z,M+1},
                      {a,n1+1},{b,n2+1},{c,n3+1}]*(-1)^(K+L+M+n1+n2+n3) /.{x->a1,y->a2,z->a3,a->a12,b->a23,c->a31}, 
            s==0, 
				int=D[hI4[-1,-1,-1,-1,-1,-1,x,y,z,a,b,c],{x,K+1},{y,L+1},{z,M+1},
                      {a,n1+1},{b,n2+1},{c,n3+1}]*(-1)^(K+L+M+n1+n2+n3) /.{x->a1,y->a2,z->a3,a->a12,b->a23,c->a31},
            s>Y2, 
				int=D[hI3[-1,-1,-1,-1,-1,-1,x,y,z,a,b,c],{x,K+1},{y,L+1},{z,M+1},
                     {a,n1+1},{b,n2+1},{c,n3+1}]*(-1)^(K+L+M+n1+n2+n3) /.{x->a1,y->a2,z->a3,a->a12,b->a23,c->a31},   
			True, 
				int=D[hI0[-1,-1,-1,-1,-1,-1,x,y,z,a,b,c],{x,K+1},{y,L+1},{z,M+1},
                      {a,n1+1},{b,n2+1},{c,n3+1}]*(-1)^(K+L+M+n1+n2+n3) /.{x->a1,y->a2,z->a3,a->a12,b->a23,c->a31}
       ];
      Re[int]
];


(* ::Subsection:: *)
(*Recurrsion Formulae*)


Clear[I0];

I0[-1,-1,-1,-1,-1,-1,a1_,a2_,a3_,a12_,a23_,a31_]=hI0[-1,-1,-1,-1,-1,-1,a1,a2,a3,a12,a23,a31];

I0[-1,-1,-1,-1,-1,n3_,a1_,a2_,a3_,a12_,a23_,a31_]:=
	I0[-1,-1,-1,-1,-1,n3,a1,a2,a3,a12,a23,a31]=
		-D[I0[-1,-1,-1,-1,-1,n3-1,a1,a2,a3,a12,a23,a31],a31];

I0[-1,-1,-1,-1,n2_,n3_,a1_,a2_,a3_,a12_,a23_,a31_]:=
	I0[-1,-1,-1,-1,n2,n3,a1,a2,a3,a12,a23,a31]=
		-D[I0[-1,-1,-1,-1,n2-1,n3,a1,a2,a3,a12,a23,a31],a23];

I0[-1,-1,-1,n1_,n2_,n3_,a1_,a2_,a3_,a12_,a23_,a31_]:=
	I0[-1,-1,-1,n1,n2,n3,a1,a2,a3,a12,a23,a31]=
		-D[I0[-1,-1,-1,n1-1,n2,n3,a1,a2,a3,a12,a23,a31],a12];

I0[-1,-1,M_,n1_,n2_,n3_,a1_,a2_,a3_,a12_,a23_,a31_]:=
	I0[-1,-1,M,n1,n2,n3,a1,a2,a3,a12,a23,a31]=
		-D[I0[-1,-1,M-1,n1,n2,n3,a1,a2,a3,a12,a23,a31],a3];

I0[-1,L_,M_,n1_,n2_,n3_,a1_,a2_,a3_,a12_,a23_,a31_]:=
	I0[-1,L,M,n1,n2,n3,a1,a2,a3,a12,a23,a31]=
		-D[I0[-1,L-1,M,n1,n2,n3,a1,a2,a3,a12,a23,a31],a2];

I0[K_,L_,M_,n1_,n2_,n3_,a1_,a2_,a3_,a12_,a23_,a31_]:=
	I0[K,L,M,n1,n2,n3,a1,a2,a3,a12,a23,a31]=
		-D[I0[K-1,L,M,n1,n2,n3,a1,a2,a3,a12,a23,a31],a1];



Clear[I2];

I2[-1,-1,-1,-1,-1,-1,a1_,a2_,a3_,a12_,a23_,a31_]=hI2[-1,-1,-1,-1,-1,-1,a1,a2,a3,a12,a23,a31];

I2[-1,-1,-1,-1,-1,n3_,a1_,a2_,a3_,a12_,a23_,a31_]:=
	I2[-1,-1,-1,-1,-1,n3,a1,a2,a3,a12,a23,a31]=
		-D[I2[-1,-1,-1,-1,-1,n3-1,a1,a2,a3,a12,a23,a31],a31];

I2[-1,-1,-1,-1,n2_,n3_,a1_,a2_,a3_,a12_,a23_,a31_]:=
	I2[-1,-1,-1,-1,n2,n3,a1,a2,a3,a12,a23,a31]=
		-D[I2[-1,-1,-1,-1,n2-1,n3,a1,a2,a3,a12,a23,a31],a23];

I2[-1,-1,-1,n1_,n2_,n3_,a1_,a2_,a3_,a12_,a23_,a31_]:=
	I2[-1,-1,-1,n1,n2,n3,a1,a2,a3,a12,a23,a31]=
		-D[I2[-1,-1,-1,n1-1,n2,n3,a1,a2,a3,a12,a23,a31],a12];

I2[-1,-1,M_,n1_,n2_,n3_,a1_,a2_,a3_,a12_,a23_,a31_]:=
	I2[-1,-1,M,n1,n2,n3,a1,a2,a3,a12,a23,a31]=
		-D[I2[-1,-1,M-1,n1,n2,n3,a1,a2,a3,a12,a23,a31],a3];

I2[-1,L_,M_,n1_,n2_,n3_,a1_,a2_,a3_,a12_,a23_,a31_]:=
	I2[-1,L,M,n1,n2,n3,a1,a2,a3,a12,a23,a31]=
		-D[I2[-1,L-1,M,n1,n2,n3,a1,a2,a3,a12,a23,a31],a2];

I2[K_,L_,M_,n1_,n2_,n3_,a1_,a2_,a3_,a12_,a23_,a31_]:=
	I2[K,L,M,n1,n2,n3,a1,a2,a3,a12,a23,a31]=
		-D[I2[K-1,L,M,n1,n2,n3,a1,a2,a3,a12,a23,a31],a1];



Clear[I3];

I3[-1,-1,-1,-1,-1,-1,a1_,a2_,a3_,a12_,a23_,a31_]=hI3[-1,-1,-1,-1,-1,-1,a1,a2,a3,a12,a23,a31];

I3[-1,-1,-1,-1,-1,n3_,a1_,a2_,a3_,a12_,a23_,a31_]:=
	I3[-1,-1,-1,-1,-1,n3,a1,a2,a3,a12,a23,a31]=
		-D[I3[-1,-1,-1,-1,-1,n3-1,a1,a2,a3,a12,a23,a31],a31];

I3[-1,-1,-1,-1,n2_,n3_,a1_,a2_,a3_,a12_,a23_,a31_]:=
	I3[-1,-1,-1,-1,n2,n3,a1,a2,a3,a12,a23,a31]=
		-D[I3[-1,-1,-1,-1,n2-1,n3,a1,a2,a3,a12,a23,a31],a23];

I3[-1,-1,-1,n1_,n2_,n3_,a1_,a2_,a3_,a12_,a23_,a31_]:=
	I3[-1,-1,-1,n1,n2,n3,a1,a2,a3,a12,a23,a31]=
		-D[I3[-1,-1,-1,n1-1,n2,n3,a1,a2,a3,a12,a23,a31],a12];

I3[-1,-1,M_,n1_,n2_,n3_,a1_,a2_,a3_,a12_,a23_,a31_]:=
	I3[-1,-1,M,n1,n2,n3,a1,a2,a3,a12,a23,a31]=
		-D[I3[-1,-1,M-1,n1,n2,n3,a1,a2,a3,a12,a23,a31],a3];

I3[-1,L_,M_,n1_,n2_,n3_,a1_,a2_,a3_,a12_,a23_,a31_]:=
	I3[-1,L,M,n1,n2,n3,a1,a2,a3,a12,a23,a31]=
		-D[I3[-1,L-1,M,n1,n2,n3,a1,a2,a3,a12,a23,a31],a2];

I3[K_,L_,M_,n1_,n2_,n3_,a1_,a2_,a3_,a12_,a23_,a31_]:=
	I3[K,L,M,n1,n2,n3,a1,a2,a3,a12,a23,a31]=
		-D[I3[K-1,L,M,n1,n2,n3,a1,a2,a3,a12,a23,a31],a1];



Clear[I4];

I4[-1,-1,-1,-1,-1,-1,a1_,a2_,a3_,a12_,a23_,a31_]=hI4[-1,-1,-1,-1,-1,-1,a1,a2,a3,a12,a23,a31];

I4[-1,-1,-1,-1,-1,n3_,a1_,a2_,a3_,a12_,a23_,a31_]:=
	I4[-1,-1,-1,-1,-1,n3,a1,a2,a3,a12,a23,a31]=
		-D[I4[-1,-1,-1,-1,-1,n3-1,a1,a2,a3,a12,a23,a31],a31];

I4[-1,-1,-1,-1,n2_,n3_,a1_,a2_,a3_,a12_,a23_,a31_]:=
	I4[-1,-1,-1,-1,n2,n3,a1,a2,a3,a12,a23,a31]=
		-D[I4[-1,-1,-1,-1,n2-1,n3,a1,a2,a3,a12,a23,a31],a23];

I4[-1,-1,-1,n1_,n2_,n3_,a1_,a2_,a3_,a12_,a23_,a31_]:=
	I4[-1,-1,-1,n1,n2,n3,a1,a2,a3,a12,a23,a31]=
		-D[I4[-1,-1,-1,n1-1,n2,n3,a1,a2,a3,a12,a23,a31],a12];

I4[-1,-1,M_,n1_,n2_,n3_,a1_,a2_,a3_,a12_,a23_,a31_]:=
	I4[-1,-1,M,n1,n2,n3,a1,a2,a3,a12,a23,a31]=
		-D[I4[-1,-1,M-1,n1,n2,n3,a1,a2,a3,a12,a23,a31],a3];

I4[-1,L_,M_,n1_,n2_,n3_,a1_,a2_,a3_,a12_,a23_,a31_]:=
	I4[-1,L,M,n1,n2,n3,a1,a2,a3,a12,a23,a31]=
		-D[I4[-1,L-1,M,n1,n2,n3,a1,a2,a3,a12,a23,a31],a2];

I4[K_,L_,M_,n1_,n2_,n3_,a1_,a2_,a3_,a12_,a23_,a31_]:=
	I4[K,L,M,n1,n2,n3,a1,a2,a3,a12,a23,a31]=
		-D[I4[K-1,L,M,n1,n2,n3,a1,a2,a3,a12,a23,a31],a1];



Clear[hIr];

hIr[K_,L_,M_,n1_,n2_,n3_,a1_,a2_,a3_,a12_,a23_,a31_]:=
	hIr[K,L,M,n1,n2,n3,a1,a2,a3,a12,a23,a31]=
	Module[
	   {s=\[Sigma]2[a1,a2,a3,a12,a23,a31],
		Y1=0.0125, Y2=-0.002,
		x,y,z,a,b,c,int},

		Which[
			s>Y1,
                int=I2[K,L,M,n1,n2,n3,x,y,z,a,b,c]/.{x->a1,y->a2,z->a3,a->a12,b->a23,c->a31}, 
            s==0, 
				int=I4[K,L,M,n1,n2,n3,x,y,z,a,b,c]/.{x->a1,y->a2,z->a3,a->a12,b->a23,c->a31}, 
            s>Y2, 
				int=I3[K,L,M,n1,n2,n3,x,y,z,a,b,c]/.{x->a1,y->a2,z->a3,a->a12,b->a23,c->a31}, 
			True, 
				int=I0[K,L,M,n1,n2,n3,x,y,z,a,b,c]/.{x->a1,y->a2,z->a3,a->a12,b->a23,c->a31}
       ];
	Re[int]   
];


(* ::Section:: *)
(*Power Series*)


(* ::Subsection:: *)
(*Auxiliary Functions*)


Clear[pA, pC, pV, pW, pU];

pA[k_,a_]:=k!/a^(k+1);

pC[v_,q_,k_]:=(2q+1)/(v+2) Binomial[v+2,2k+1]\!\(
\*UnderoverscriptBox[\(\[Product]\), \(t = 0\), \(Min[q - 1, 
\*FractionBox[\(v + 1\), \(2\)]]\)]
\*FractionBox[\(2  k + 2  t - v\), \(2  k + 2  q - 2  t + 1\)]\);

pV[k_,l_,a_,b_]/;l>=0 := pV[k,l,a,b] =\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(lp = 0\), \(l\)]\((Binomial[l, lp] pA[lp, b] pA[k + l - lp, a + b])\)\);
pV[k_,l_,a_,b_]/;l<0 := pV[k,l,a,b] = pA[k+l+1,a+b]/(k+1) Hypergeometric2F1[1,k+l+2,k+2,a/(a+b)];

pW[k_,l_,m_,a_,b_,c_]/;m>=0 := pW[k,l,m,a,b,c] = \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(mp = 0\), \(m\)]\((Binomial[m, mp] pA[mp, c] pV[k, l + m - mp, a, b + c])\)\);
pW[k_,l_,m_,a_,b_,c_]/;m<0 := pW[k,l,m,a,b,c] = pA[k+l+m+2,b]\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(lp = 0\), \(l\)]\((Binomial[l, lp] 
\*SuperscriptBox[\((\(-1\))\), \(lp\)] 
\*FractionBox[
SuperscriptBox[\(b\), \(k + m + 2\)], \(k + lp + m + 2\)]*\n\t\t\((
\*FractionBox[
SuperscriptBox[\(a\), \(lp\)], 
SuperscriptBox[\((a + b)\), \(k + lp + m + 2\)]] pU[k, lp, k + lp + m + 2, 
\*FractionBox[\(c\), \(a\)], 
\*FractionBox[\(c\), \(a + b\)]] - 
\*FractionBox[\(1\), 
SuperscriptBox[\(a\), \(k + m + 2\)]] pU[k, lp, k + lp + m + 2, 
\*FractionBox[\(c\), \(a\)], 
\*FractionBox[\(b + c\), \(a\)]])\))\)\);

pU[p_,q_,r_,x_,y_]:= pU[p,q,r,x,y] = \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(qp = 0\), \(q\)]\((Binomial[q, qp] 
\*FractionBox[
SuperscriptBox[\(x\), \(q - qp\)], 
SuperscriptBox[\((y + 1)\), \(r\)]] 
\*FractionBox[\(1\), \(p + qp + 1\)] Hypergeometric2F1[1, r, p + qp + 2, 
\*FractionBox[\(1\), \(1 + y\)]])\)\);



(* ::Subsection:: *)
(*Working Formulae*)


Clear[I3ep];

(* l, m, n are even *)
I3ep[K_,L_,M_,l_?EvenQ,m_?EvenQ,n_?EvenQ,a1_,a2_,a3_]:= 
	I3ep[K,L,M,l,m,n,a1,a2,a3]=
		64\[Pi]^3*\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(q = 0\), \(Min[l/2, m/2, n/2]\)]\((
\*FractionBox[\(1\), 
SuperscriptBox[\((2  q + 1)\), \(2\)]] \(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 0\), \(l/2 - q\)]\((pC[l, q, i]*\n\t\t\t\(
\*UnderoverscriptBox[\(\[Sum]\), \(j = 0\), \(m/2 - q\)]\((pC[m, q, j] \(
\*UnderoverscriptBox[\(\[Sum]\), \(k = 0\), \(n/2 - q\)]\((pC[n, q, k] pA[K + 2 + m - 2  j + 2  k, a1] pA[L + 2 + 2  i + n - 2  k, a2] pA[M + 2 + l - 2  i + 2  j, a3])\)\))\)\))\)\))\)\);


(* any of l, m, n is odd *)
I3ep[K_,L_,M_,l_?OddQ,m_?EvenQ,n_?EvenQ,a1_,a2_,a3_]:= 
	I3ep[K,L,M,l,m,n,a1,a2,a3]=
		64\[Pi]^3*\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(q = 0\), \(Min[m/2, n/2]\)]\((
\*FractionBox[\(1\), 
SuperscriptBox[\((2  q + 1)\), \(2\)]] \(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 0\), \(\((l + 1)\)/2\)]\((pC[l, q, i]*\n\t\t\t\t\t\t\ \ \ \t\(
\*UnderoverscriptBox[\(\[Sum]\), \(j = 0\), \(m/2 - q\)]\((pC[m, q, j] \(
\*UnderoverscriptBox[\(\[Sum]\), \(k = 0\), \(n/2 - q\)]\((pC[n, q, k] pA[K + 2 + m - 2  j + 2  k, a1]*\n\t\t\t\t\t\t\t\t\t\((pV[L + 2 + 2  i + n - 2  k, M + 2 + l - 2  i + 2  j, a2, a3] + \n\t\t\t\t\t\t\t\t\t\ pV[M + 2 + 2  q + 2  i + 2  j, L + 2 - 2  q + l - 2  i + n - 2  k, a3, a2])\))\)\))\)\))\)\))\)\);


I3ep[K_,L_,M_,l_?EvenQ,m_?OddQ,n_?EvenQ,a1_,a2_,a3_]:= I3ep[L,K,M,m,l,n,a2,a1,a3];

I3ep[K_,L_,M_,l_?EvenQ,m_?EvenQ,n_?OddQ,a1_,a2_,a3_]:= I3ep[M,L,K,n,m,l,a3,a2,a1];


(* any of l, m, n is even *)
I3ep[K_,L_,M_,l_?OddQ,m_?OddQ,n_?EvenQ,a1_,a2_,a3_]:= 
	I3ep[K,L,M,l,m,n,a1,a2,a3]=
		64\[Pi]^3*\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(q = 0\), \(n/2\)]\((
\*FractionBox[\(1\), 
SuperscriptBox[\((2  q + 1)\), \(2\)]] \(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 0\), \(\((l + 1)\)/2\)]\((pC[l, q, i] \(
\*UnderoverscriptBox[\(\[Sum]\), \(j = 0\), \(\((m + 1)\)/2\)]\((pC[m, q, j] \(
\*UnderoverscriptBox[\(\[Sum]\), \(k = 0\), \(n/2 - q\)]\((pC[n, q, k]*\n\t\t\t\((pW[K + 2 + 2  q + 2  j + 2  k, L + 2 + 2  i + n - 2  k, M + 2 - 2  q + l - 2  i + m - 2  j, a1, a2, a3] + \n\t\t\t\ pW[K + 2 + 2  q + 2  j + 2  k, M + 2 + 2  i + m - 2  j, L + 2 - 2  q + l - 2  i + n - 2  k, a1, a3, a2] + \n\t\t\t\ pW[L + 2 + 2  q + 2  i + 2  k, K + 2 + 2  j + n - 2  k, M + 2 - 2  q + l - 2  i + m - 2  j, a2, a1, a3] + \n\t\t\t\ pW[L + 2 + 2  q + 2  i + 2  k, M + 2 + 2  j + l - 2  i, K + 2 - 2  q + m - 2  j + n - 2  k, a2, a3, a1] + \n\t\t\t\ pW[M + 2 + 2  q + 2  i + 2  j, K + 2 + 2  k + m - 2  j, L + 2 - 2  q + l - 2  i + n - 2  k, a3, a1, a2] + \n\t\t\t\ pW[M + 2 + 2  q + 2  i + 2  j, L + 2 + 2  k + l - 2  i, K + 2 - 2  q + m - 2  j + n - 2  k, a3, a2, a1])\))\)\)\n\t\t)\)\))\)\))\)\);

I3ep[K_,L_,M_,l_?OddQ,m_?EvenQ,n_?OddQ,a1_,a2_,a3_]:=I3ep[K,M,L,l,n,m,a1,a3,a2];
I3ep[K_,L_,M_,l_?EvenQ,m_?OddQ,n_?OddQ,a1_,a2_,a3_]:=I3ep[M,L,K,n,m,l,a3,a2,a1];


(* l, m, n are odd
   the sum is an infinite power series.

   The infinite power series should be computed using accelerator approach, e.g.,
   Levin u-transformation (J. S. Sims and S. A. Hagstrom, J. Phys. B, 37, 1519 (2004)

   method = 0: direct sum
          = 1: Levin u transformation 
*)
I3ep[K_,L_,M_,l_?OddQ,m_?OddQ,n_?OddQ,a1_,a2_,a3_,kmax_:30,method_:1]:=
	I3ep[K,L,M,l,m,n,a1,a2,a3,kmax,method]=
	Module[
		{A,S,C,int},

		Do[A[q] = 1/(2q+1)^2 \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 0\), \(\((l + 1)\)/2\)]\((pC[l, q, i] \(
\*UnderoverscriptBox[\(\[Sum]\), \(j = 0\), \(\((m + 1)\)/2\)]\((pC[m, q, j] \(
\*UnderoverscriptBox[\(\[Sum]\), \(k = 0\), \(\((n + 1)\)/2\)]\((pC[n, q, k]*\n\t\t\t\((pW[K + 2 + 2  q + 2  j + 2  k, L + 2 + 2  i + n - 2  k, M + 2 - 2  q + l - 2  i + m - 2  j, a1, a2, a3] + \n\t\t\t\ pW[K + 2 + 2  q + 2  j + 2  k, M + 2 + 2  i + m - 2  j, L + 2 - 2  q + l - 2  i + n - 2  k, a1, a3, a2] + \n\t\t\t\ pW[L + 2 + 2  q + 2  i + 2  k, K + 2 + 2  j + n - 2  k, M + 2 - 2  q + l - 2  i + m - 2  j, a2, a1, a3] + \n\t\t\t\ pW[L + 2 + 2  q + 2  i + 2  k, M + 2 + 2  j + l - 2  i, K + 2 - 2  q + m - 2  j + n - 2  k, a2, a3, a1] + \n\t\t\t\ pW[M + 2 + 2  q + 2  i + 2  j, K + 2 + 2  k + m - 2  j, L + 2 - 2  q + l - 2  i + n - 2  k, a3, a1, a2] + \n\t\t\t\ pW[M + 2 + 2  q + 2  i + 2  j, L + 2 + 2  k + l - 2  i, K + 2 - 2  q + m - 2  j + n - 2  k, a3, a2, a1])\))\)\))\)\))\)\)
		,{q,0,kmax}];

		Which[method == 0, 
				int = 64\[Pi]^3*\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(q = 0\), \(kmax\)]\(A[q]\)\),
			  method == 1,
				S[0]=A[0];
				Do[S[j]=S[j-1]+A[j],{j,1,kmax}];

				Do[C[j]=(-1)^j Binomial[kmax,j](j+1)^(kmax-2) 1/A[j],{j,0,kmax}];
				int = 64\[Pi]^3*\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(j = 0\), \(kmax\)]\((C[j]*S[j])\)\)/\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(j = 0\), \(kmax\)]\(C[j]\)\)
		];
		int
	];



(* ::Subsection:: *)
(*Recurrence Formulae*)


Clear[pVr];
pVr[m_,0,a_,b_]/;m>=0:=pVr[m,0,a,b]=1/b pA[m,a+b];
pVr[m_,n_,a_,b_]/;(m>=0 && m+n>=0):=pVr[m,n,a,b]=1/b (n*pVr[m,n-1,a,b]+pA[m+n,a+b]);

Clear[pWr];

pWr[0,0,0,a_,b_,r_]:=pWr[0,0,0,a,b,r]=1/((a+b+r)(b+r)r);
pWr[f_,g_,h_,a_,b_,r_]/;(f>=0 && f+g>=-1 && f+g+h>=-1):=pWr[f,g,h,a,b,r]=1/r (h*pWr[f,g,h-1,a,b,r]+pVr[f,g+h,a,b+r]);
pWr[f_,g_,h_,a_,b_,r_]/;(f>=0 && f+g>=0 && f+g+h>=-1):=pWr[f,g,h,a,b,r]=1/(b+r) (g*pWr[f,g-1,h,a,b,r]+h*pWr[f,g,h-1,a,b,r]+pVr[f+g,h,a+b,r]);
pWr[f_,g_,h_,a_,b_,r_]/;(f>=1 && f+g>=0 && f+g+h>=-1):=pWr[f,g,h,a,b,r]=1/(a+b+r) (f*pWr[f-1,g,h,a,b,r]+g*pWr[f,g-1,h,a,b,r]+h*pWr[f,g,h-1,a,b,r]);



(* ::Section:: *)
(*Some analytical formulae for I*)


Clear[Dw,Dhat];
Dw[wi_,wj_,wk_]:=Log[(wj+wk)/wi]Log[(wi+wj+wk)/(wj+wk)]-PolyLog[2,-(wi/(wj+wk))]-PolyLog[2,1-wi/(wj+wk)];
Dhat[w1_,w2_,w3_]:=1/2 (Dw[w1,w2,w3]+Dw[w2,w3,w1]+Dw[w3,w1,w2]);

Clear[I3e];

(* ------------------------------------------------------------------- *)
I3e[-1,-1,-1,-1,-1,-1,w1_,w2_,w3_]:=64\[Pi]^3*Dhat[w1,w2,w3]/(w1*w2*w3);

I3e[K_,L_,M_,-1,-1,-1,w1_,w2_,w3_]/;(K+L+M)<9 :=Module[
	{a1,a2,a3},
	(-1)^(K+L+M+1)*D[I3e[-1,-1,-1,-1,-1,-1,a1,a2,a3],
				{a1,K+1},{a2,L+1},{a3,M+1}]/.{a1->w1,a2->w2,a3->w3}
];

(* ------------------------------------------------------------------- *)
I3e[-1,-1,-1,0,0,-1,w1_,w2_,w3_]:=
	64\[Pi]^3*1/(w1*w2*(w1+w2)*w3^2);

I3e[K_,L_,M_,0,0,-1,w1_,w2_,w3_]/;(K+L+M)<9 :=Module[
	{a1,a2,a3},
	(-1)^(K+L+M+1)*D[I3e[-1,-1,-1,0,0,-1,a1,a2,a3],
				{a1,K+1},{a2,L+1},{a3,M+1}]/.{a1->w1,a2->w2,a3->w3}
];

(* ------------------------------------------------------------------- *)
I3e[-1,-1,-1,1,-1,-1,w1_,w2_,w3_]:=
	64\[Pi]^3*1/(w2^2 w3^2) ((w2^2+w3^2-w1^2)/(w1*w2*w3) Dhat[w1,w2,w3]-1/(w1+w2)-1/(w1+w3)+
				(w1-w2)/(w1*w2) Log[(w1+w2)/w3]+(w1-w3)/(w1*w3) Log[(w1+w3)/w2]+(w2+w3)/(w2*w3) Log[(w2+w3)/w1]);

I3e[K_,L_,M_,1,-1,-1,w1_,w2_,w3_]/;(K+L+M)<9 :=Module[
	{a1,a2,a3},
	(-1)^(K+L+M+1)*D[I3e[-1,-1,-1,1,-1,-1,a1,a2,a3],
				{a1,K+1},{a2,L+1},{a3,M+1}]/.{a1->w1,a2->w2,a3->w3}
];

(* ------------------------------------------------------------------- *)
I3e[-1,-1,-1,1,-1,1,w1_,w2_,w3_]:=
	64\[Pi]^3*1/(w1^3 w2^5 w3^3) ((w2^2 (2w1^2+w2^2+2w3^2)-3(w1^2-w3^2)^2)Dhat[w1,w2,w3]+
					w1((w3-w2)(w2^2+3w3^2)+3w1^2 (w2+w3))Log[(w2+w3)/w1]+
					w3((w1-w2)(w2^2+3w1^2)+3w3^2 (w1+w2))Log[(w1+w2)/w3]+
					w2(w1+w3)(3(w1-w3)^2-w2^2)Log[(w1+w3)/w2]+
					2w1*w2*w3*((w1*w2^2*w3)/(w1+w3)^3-2(w1+w3)+w2+(w2^2+3w1*w3)/(w1+w3)-w1^2/(w1+w2)-w3^2/(w2+w3)));

I3e[K_,L_,M_,1,-1,1,w1_,w2_,w3_]/;(K+L+M)<9 :=Module[
	{a1,a2,a3},
	(-1)^(K+L+M+1)*D[I3e[-1,-1,-1,1,-1,1,a1,a2,a3],
				{a1,K+1},{a2,L+1},{a3,M+1}]/.{a1->w1,a2->w2,a3->w3}
];

(* ------------------------------------------------------------------- *)
I3e[-1,-1,-1,1,1,-1,w1_,w2_,w3_]:=I3e[-1,-1,-1,1,-1,1,w1,w3,w2];

I3e[K_,L_,M_,1,1,-1,w1_,w2_,w3_]/;(K+L+M)<9 :=Module[
	{a1,a2,a3},
	(-1)^(K+L+M+1)*D[I3e[-1,-1,-1,1,1,-1,a1,a2,a3],
				{a1,K+1},{a2,L+1},{a3,M+1}]/.{a1->w1,a2->w2,a3->w3}
];

(* ------------------------------------------------------------------- *)


(* ::Section:: *)
(*Final Working Formulae*)



I3e[K_,L_,M_,l_,m_,n_,a1_,a2_,a3_]:=
	I3e[K,L,M,l,m,n,a1,a2,a3]=
	I3ep[K,L,M,l,m,n,a1,a2,a3];


(* ::Section:: *)
(*Another Formulae*)


(* D. M. Fromm and R. N. Hill, Phys. Rev. A, 36, 1013 (1987) *)

Clear[fhv];

fhv[z_]:=1/2 PolyLog[2,(1-z)/2]-1/4 (Log[(1-z)/2])^2 -1/2 PolyLog[2,(1+z)/2]+1/4(Log[(1+z)/2]^2);

fhv[x_/;x>1]:=1/2 PolyLog[2,2/(1+x)]+1/2 Log[2/(1+x)]^2-1/2 PolyLog[2,2/(1-x)]-1/2 Log[2/(x-1)]^2+I/2 \[Pi]*Log[(x-1)/(x+1)];

fhv[x_/;x<-1]:=1/2 PolyLog[2,2/(1+x)]+1/2 Log[2/(-1-x)]^2-1/2 PolyLog[2,2/(1-x)]-1/2 Log[2/(1-x)]^2-I/2 \[Pi]*Log[(1-x)/(-1-x)];


(* F. E. Harris, Phys. Rev. A, 55, 1820 (1997) *)

Clear[fhI];

fhI[-1,-1,-1,-1,-1,-1,a1_,a2_,a3_,a12_,a23_,a31_]:=
	fhI[-1,-1,-1,-1,-1,-1,a1,a2,a3,a12,a23,a31]=
    (16\[Pi]^3)/hs[a1,a2,a3,a12,a23,a31] (-2Sum[fhv[hG[j,a1,a2,a3,a12,a23,a31]/hs[a1,a2,a3,a12,a23,a31]],{j,1,3}]
                               +Sum[fhv[hg[k,j,a1,a2,a3,a12,a23,a31]/hs[a1,a2,a3,a12,a23,a31]],{j,0,3},{k,0,3}]);
