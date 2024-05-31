(*$R+*)


(********************************************************
*********************************************************

          SHEAR-BAND INSTABILITY PROGRAM


   This program calculates and plots a strain-rate/
temperature map showing the domain where 'adiabatic'
shear-band instabilities are expected to develope
using both nearly-adiabatic and nearly-isothermal
criteria. Contours of critical strain are drawn on
the map

P.M.Sargent

19 May 1981
                                    Materials Laboratory,
                                      Dept.  Engineering,
                                         Trumpington St.,
                                        Cambridge CB2 1PZ

*********************************************************
********************************************************)

PROCEDURE identity(VAR id: text_string); EXTERN;

(*          This function is used to identify which    *)
(*          strength relationship is being used by     *)
(*          the program.                               *)

(*
---------------------------------------------------------
*)
FUNCTION strength(strain,s_rate,t: REAL): REAL; EXTERN;

(*          This function returns the flow stress of   *)
(*          the material at the given strain, strain   *)
(*          rate and temperature.  It completely       *)
(*          ignores history effects.                   *)

(*
---------------------------------------------------------
*)
FUNCTION s_h_e(strain,t: REAL): REAL; EXTERN;

(*          This function returns the                  *)
(*          strain-hardening exponent (it assumes      *)
(*          that the Hollomon relationship holds at    *)
(*          any temperature and strain rate). In       *)
(*          fact, all versions of this function as     *)
(*          yet in use take the strain hardening as a  *)
(*          seperable term in the routine 'strength'   *)
(*          and therefore s_h_e is set as a constant   *)
(*          for all temperatures and strain rates.     *)

(*
---------------------------------------------------------
*)
FUNCTION t_softening(strain,s_rate,t: REAL): REAL;

(*          This function returns the value of         *)
(*          d(strength)/d(temperature), the thermal    *)
(*          softening coefficient, at given values of  *)
(*          strain, strain rate and temperature. It    *)
(*          performs a simple numerical differencing.  *)

CONST
   very_small = 1.0E-10;

VAR
   incr: REAL;
   sg1,sg2 : REAL;


BEGIN

IF t>0.01 THEN incr:=0.01 ELSE incr:=t;
sg1:=strength(strain,s_rate,t-incr);
sg2:=strength(strain,s_rate,t+incr);
IF sg1<>sg2 THEN t_softening:=(sg1-sg2)/(2*incr)
ELSE t_softening:= very_small
END;

(*
---------------------------------------------------------
*)
PROCEDURE integ(z,she,epsabs,epsrel:
                 REAL; VAR ivalue: REAL); FORTRAN;

(*          This procedure is an external fortran      *)
(*          subroutine which calls Naglib routine      *)
(*          D01BDF to integrate the function           *)
(*          'exp(-z(1-y))*y**she' between 0.0 and 1.0  *)
(*          (where z=strain/(tauit*s_rate) and       *)
(*          y=x/strain; x is the variable of           *)
(*          integration varying between zero and the   *)
(*          strain at which the objective function is  *)
(*          being evaluated.                           *)

(*
---------------------------------------------------------
*)
FUNCTION objective(strain,s_rate,t,cp: REAL): REAL;

(*              The objective function is zero when    *)
(*          the variable (the strain) is equal to the  *)
(*          critical strain for 'adiabatic' shear      *)
(*          instability to occur at the given          *)
(*          temperature, strain-rate and specific      *)
(*          heat.                                      *)

EXTERN;

(*
---------------------------------------------------------
*)
FUNCTION critical_strain(s_incr,acc,guess,
                      s_rate,t,cp: REAL): REAL; EXTERN;

(*                           This function iterates    *)
(*          to find the zero of the objective          *)
(*          function given a particular strain rate    *)
(*          and temperature.                           *)

(*
---------------------------------------------------------
*)
FUNCTION sr_critical(sr_guess,fixed_strain,tau,t,cp:
                    REAL): REAL; EXTERN;

(*          This function iterates on the strain rate  *)
(*          to find the zero of the objective          *)
(*          function given a particular strain and     *)
(*          temperature.                               *)

(*
---------------------------------------------------------
*)
FUNCTION capacity(m : metal; t: REAL): REAL;

(*          This function returns the thermal          *)
(*          capacity of the metal at the given         *)
(*          temperature by making a linear             *)
(*          interpolation between the data points on   *)
(*          either side of this temperature. The       *)
(*          table of data points is supplied in        *)
(*          m.capacity .                               *)

VAR
   j: 0..100;

BEGIN
WITH m.capacity DO
   BEGIN j:=0;
   WHILE (prop_tempsj+1<t) AND (j<l) DO
      BEGIN
      j:=j+1
      END;
   IF j=0 THEN capacity:=misc_value*t;
   IF (j>0) AND (j+1<=l) THEN capacity:=propertyj
                     + (t-prop_tempsj)
                     *(propertyj+1 - propertyj)
                     /(prop_tempsj+1 - prop_tempsj);
   IF j=l THEN capacity:=propertyl
   END
END;

(*
---------------------------------------------------------
*)
FUNCTION conductivity(m : metal; t: REAL): REAL;

(*          This function returns the thermal          *)
(*          conductivity of the metal at the given     *)
(*          temperature by making a linear             *)
(*          interpolation between the data points      *)
(*          above and below this temperature. The      *)
(*          table of data points is supplied in        *)
(*          m.conductivity .                           *)

VAR
   j: 0..100;

BEGIN
WITH m.conductivity DO
   BEGIN j:=0;
   WHILE (prop_tempsj+1<t) AND (j<l) DO
      BEGIN
      j:=j+1
      END;
   IF j=0 THEN conductivity:=property1;
   IF (j>0) AND (j+1<=l) THEN conductivity:=propertyj
                     + (t-prop_tempsj)
                     *(propertyj+1 - propertyj)
                     /(prop_tempsj+1 - prop_tempsj);
   IF j=l THEN conductivity:=propertyl
   END
END;

(*
---------------------------------------------------------
*)
FUNCTION double_iter(
         FUNCTION calc_s_rate(d1,d2,d3: REAL): REAL;
         VAR guess: REAL; sr_guess,tau,t,cp: REAL): REAL;


(*              The strain_rates for                   *)
(*          nearly-isothermal and nearly-adiabatic     *)
(*          conditions are strongly dependent on the   *)
(*          critical strain values. The critical       *)
(*          strain values are weakly dependent on the  *)
(*          strain-rate. Thus the strain-rates         *)
(*          generally converge after only two          *)
(*          iterations once the critical strain has    *)
(*          been calculated using a guessed            *)
(*          strain-rate (sr_guess).                    *)

VAR
   sr0,sr1,strain: REAL;

BEGIN
IF flag_print THEN
   WRITELN('double_iter',guess:12,sr_guess:12);
strain:=critical_strain(0.01,1.0E-3,guess,
            sr_guess,t,cp);
IF strain <  twenty THEN
   BEGIN
   sr1:=calc_s_rate(strain,tau,t);
   REPEAT
      BEGIN
      sr0:=sr1;
      strain:=critical_strain(0.008,acc,strain,sr1,t,cp);
      sr1:=calc_s_rate(strain,tau,t);
      IF strain=twenty THEN
         BEGIN sr0:=max_sr; sr1:=max_sr END
      END
   UNTIL (sr0/sr1 < acc_sr) AND (sr0/sr1 > 1/acc_sr);
   double_iter:=(sr1+sr0)/2
   END
ELSE double_iter:=max_sr;

guess:=strain
END;

(*
---------------------------------------------------------
*)
FUNCTION srf_adiabatic(strain,tau,t: REAL): REAL;

(*          This function calculates the               *)
(*          nearly-adiabatic strain from only the      *)
(*          critical strain and the thermal            *)
(*          diffusivity. It ignores the strain-rate    *)
(*          sensitivity of the mechanical properties   *)
(*          (which is taken care of by the iteration   *)
(*          in double_iter).                           *)


VAR
   z: REAL;

BEGIN
z:=strain*50/tau;
IF z < m.epsilondot_zero THEN
   srf_adiabatic:=z
ELSE srf_adiabatic:=m.epsilondot_zero
END;


(*
---------------------------------------------------------
*)
FUNCTION srf_isothermal(strain,tau,t: REAL): REAL;

(*          Calculates the nearly-isothermal           *)
(*          strain-rate from only the given critical   *)
(*          strain and the thermal diffusivity. It     *)
(*          ignores the strain-rate sensitivity of     *)
(*          the mechanical properties (which is taken  *)
(*          care of by double_iter).                   *)

CONST
   e = 2.71828181;

VAR
   z: REAL;

BEGIN
z:=strain*e/tau;
IF z<m.epsilondot_zero THEN
   srf_isothermal:=z
ELSE srf_isothermal:=m.epsilondot_zero
END;
(*
---------------------------------------------------------
*)
FUNCTION truth(c: CHAR): BOOLEAN;

BEGIN
CASE c OF
   '+','1','T','t' : truth:=TRUE;
   '-','0','F','f' : truth:=FALSE;
   ELSE            : HALT
   END
END;

(*
---------------------------------------------------------
---------------------------------------------------------

          High Level graphics Routines

           (Fortran double precision)

---------------------------------------------------------
---------------------------------------------------------
*)
PROCEDURE gran2d(iflen,ifdig: INTEGER); FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE gran3d(iflen,ifdig: INTEGER); FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE gran4d(xstart,ystart,astart,
                 xstop,ystop,astop,
                 base,ainc: REAL;
                 nos: INTEGER;
                 a,b,c,d,xc,yc:REAL;
                 iflen,ifdig: INTEGER); FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE gran5d(x,y: REAL;
                 string: text_string;
                 len: INTEGER;
                 rot: REAL); FORTRAN;

(*
---------------------------------------------------------
*)
PROCEDURE gran6d(tmain: text_string; lmain: INTEGER;
                 tx   : text_string; lx   : INTEGER;
                 tl   : text_string; ly   : INTEGER);
                 FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE grfn2d(pmin,pmax: REAL; numpts: INTEGER;
                 PROCEDURE subr); FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE grfr2d; FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE grfr3d; FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE grfr4d; FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE grfr5d; FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE grfr6d(isw: INTEGER); FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE grft2d(iframe,imark,iannot,icorn: INTEGER);
                                             FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE grft3d(pmaj,pmin: REAL; nmaj,ntot: INTEGER;
                 val: REAL; nsub: INTEGER); FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE grft4d(pmaj,pmin: REAL; nmaj,ntot: INTEGER;
                 val: REAL; nsub: INTEGER); FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE grgr6d(x,y: array0; n: INTEGER);  FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE grgr7d(x,y: array0; n: INTEGER); FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE grgr8d(x,y: array0; n: INTEGER); FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE grgr9d(x,y: array0; n: INTEGER); FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE grlm2d(xlow,xhigh,ylow,yhigh: REAL); FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE grlm3d(xlow,xhigh,ylow,yhigh: REAL); FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE grlm4d(VAR xlow,xhigh,ylow,yhigh: REAL);
                                          FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE grlm5d(VAR xlow,xhigh,ylow,yhigh: REAL);
                                          FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE grlm6d; FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE grlm7d(irot: INTEGER); FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE grlm8d(VAR irot: INTEGER); FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE grms2d(xbase,ybase: REAL); FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE grms3d(alen: REAL ); FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE grms4d(obslow,obshigh: REAL; num: INTEGER;
                 base: REAL; VAR rlow,rhigh,rinc: REAL);
                                                FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE grms5d(alen: REAL); FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE grpn2d(icol: INTEGER); FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE grpn3d(icol: INTEGER); FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE grpn4d(icol: INTEGER); FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE grpn5d(icol: INTEGER); FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE grty4d(errgap: REAL; minpts: INTEGER); FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE grty5d(ifreq,ichar: INTEGER); FORTRAN;
(*
---------------------------------------------------------
---------------------------------------------------------

          Low Level Graphics Routines

---------------------------------------------------------
---------------------------------------------------------
*)
PROCEDURE option(i,j,k: INTEGER); FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE moveto(x,y: REAL); FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE moveby(x,y: REAL); FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE drawto(x,y: REAL); FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE drawby(x,y: REAL); FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE pltchr(n: INTEGER); FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE trmchr(a,b,c,d: REAL); FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE sclchr(x,y: REAL); FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE rotchr(a: REAL); FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE brkplt(n: INTEGER); FORTRAN;
(*
---------------------------------------------------------
*)
PROCEDURE origin(x,y: REAL; n: INTEGER); FORTRAN;
(*
---------------------------------------------------------
---------------------------------------------------------

          Local Style Graphics Procedures

---------------------------------------------------------
---------------------------------------------------------
*)
PROCEDURE pl_init;

(*          Does a brkplt while supressing the         *)
(*          'Fortran Unit 8' message                   *)

BEGIN
trmchr(0.0,0.0,0.0,0.0);
brkplt(8);
trmchr(3.5,0.0,0.0,5.0)
END;

(*
---------------------------------------------------------
*)
PROCEDURE pl_identity;

VAR
   id: text_string;

BEGIN
trmchr(2.8,1.4,0.0,3.5);
grpn5d(3);
identity(id);
gran5d(70.0,25.0,id,n3,0.0);
trmchr(3.5,0.0,0.0,5.0)
END;

(*
---------------------------------------------------------
*)
PROCEDURE pl_label;

BEGIN
trmchr(2.8,1.4,0.0,3.5);
grpn5d(3);
gran5d(70.0,25.0,m.name2,n3,0.0);
trmchr(3.5,0.0,0.0,5.0)
END;

(*
---------------------------------------------------------
*)
PROCEDURE pl_set;

(*          Sets paper limitis (in mm) for the graph.  *)

BEGIN
grlm2d(70.0,240.0,50.0,220.0)
END;

(*
---------------------------------------------------------
*)
PROCEDURE pl_reset;

(*          Sets user limits for the graph to be the   *)
(*          same as the paper limits set in 'pl_set'   *)
(*          - effectively restoring the scales to mm.  *)

BEGIN
gran2d(0,0);
gran3d(0,0);
grlm3d(70.0,240.0,50.0,220.0)
END;

(*
---------------------------------------------------------
*)
PROCEDURE pl_marks;

(*          Sets lengs and sense of the major and      *)
(*          minor axes marks and sets the frame        *)
(*          colour to black.                           *)

BEGIN
grpn3d(1);
grft3d(-3.0,-1.5,5,10,0.0,0);
grft4d(-3.0,-1.5,5,20,0.0,0)
END;

(*
---------------------------------------------------------
*)
PROCEDURE rmxmn(a: array0;
                n: index;
            VAR mx,mn: REAL);

(*          Procedure returns the maximum and minimum  *)
(*          values from an array of type               *)
(*          'array0'.                             *)

VAR
   j: index;

BEGIN
mx:=a0;
mn:=mx;
FOR j:=1 TO n DO
   BEGIN
   IF mn>aj THEN mn:=aj;
   IF mx<aj THEN mx:=aj
   END
END;

(*
---------------------------------------------------------
*)

PROCEDURE pl_asterisk;

VAR
   t: text_string;

BEGIN
trmchr(6.5,0.0,0.0,10.0);
grpn5d(2);
t1:='*';
gran5d(7.0,110.0,t,1,0.0)
END;

(*
---------------------------------------------------------
*)

PROCEDURE pl_array(z: array0);

VAR
   mx,mn : REAL;

BEGIN
rmxmn(z,t_n,mx,mn);
grlm3d(t_min,t_max,mn,mx);
grgr6d(t_array,z,t_n+1)
END;

(*
---------------------------------------------------------
*)
PROCEDURE pl_a4;

(*          Plots four characters (crosses for c=78)   *)
(*          delineating an A4 tectangle.               *)

CONST
   c = 78;

BEGIN
moveto(10.0,10.0);
pltchr(c);
moveto(10.0,271.0);
pltchr(c);
moveto(213.0,271.0);
pltchr(c);
moveto(213.0,10.0);
pltchr(c)
END;

(*
---------------------------------------------------------
*)
PROCEDURE pl_r_title;

(*          To write out a text string on the plot     *)
(*          the string has to be copied into a         *)
(*          standard length array, the same type as    *)
(*          declared in the gran5d or gran6d           *)
(*          declaration.                               *)

CONST
   la = 26;
   lb = 11;
   lc = 13;

VAR
   a,b,c : text_string;
   i     : index3;
   a2    : ARRAY1..la OF CHAR;
   b2    : ARRAY1..lb OF CHAR;
   c2    : ARRAY1..lc OF CHAR;

BEGIN
a2:='Shear Band Instability Map';
b2:='strain-rate';
c2:='temperature/K';

FOR i:=1 TO la DO ai:=a2i;
FOR i:=1 TO lb DO bi:=b2i;
FOR i:=1 TO lc DO ci:=c2i;

gran6d(a,la,c,lc,b,lb)

END;
(*
---------------------------------------------------------
*)
PROCEDURE pl_tau_title;

CONST
   la = 35;
   lb = 11;
   lc = 13;

VAR
   a,b,c : text_string;
   i     : index3;
   a2    : ARRAY1..la OF CHAR;
   b2    : ARRAY1..lb OF CHAR;
   c2    : ARRAY1..lc OF CHAR;

BEGIN
a2:='Time Constant for Thermal Diffusion';
b2:='tau/seconds';
c2:='temperature/K';

FOR i:=1 TO la DO ai:=a2i;
FOR i:=1 TO lb DO bi:=b2i;
FOR i:=1 TO lc DO ci:=c2i;

gran6d(a,la,c,lc,b,lb)

END;
(*
---------------------------------------------------------
*)
PROCEDURE pl_y_title;

CONST
   la = 43;
   lb = 29;
   lc = 13;

VAR
   a,b,c : text_string;
   i     : index3;
   a2    : ARRAY1..la OF CHAR;
   b2    : ARRAY1..lb OF CHAR;
   c2    : ARRAY1..lc OF CHAR;

BEGIN
a2:='Flow Stress at Adiabatic Instability Strain';
b2:='(strength/thermal capacity)/K';
c2:='temperature/K';

FOR i:=1 TO la DO ai:=a2i;
FOR i:=1 TO lb DO bi:=b2i;
FOR i:=1 TO lc DO ci:=c2i;

gran6d(a,la,c,lc,b,lb)

END;
(*
---------------------------------------------------------
*)
PROCEDURE pl_parms;

(*          There is no easy way to write numbers      *)
(*          onto a plot using pascal. Pascal can only  *)
(*          send arrays of characters to be written    *)
(*          as strings. The technique used here is to  *)
(*          write the numbers to an internal file,     *)
(*          then to read from that file into           *)
(*          character strings. These strings are then  *)
(*          written onto the plot using gran5d.        *)

VAR
   y       : REAL;
   a,b,c,d : text_string;
   id      : text_string;
   r       : ARRAY1..n1 OF text_string;

BEGIN
RESET(buff);
WRITELN(buff,'acc_sr = ',acc_sr:5:2);
WRITELN(buff,'acc    = ',acc:9);
WRITELN(buff,'t_n    = ',t_n:5);
WRITELN(buff,'ar     = ',ar:9);


FOR i:=1 TO strains_n DO
   WRITELN(buff,strains_arrayi:5:2);

RESET(buff);
READLN(buff,a);
READLN(buff,b);
READLN(buff,c);
READLN(buff,d);

FOR i:=1 TO strains_n DO
   READLN(buff,ri);

trmchr(6.5,0.0,0.0,10.0);
grpn5d(2);
gran5d(70.0,300.0,m.name,10,0.0);
trmchr(3.5,0.0,0.0,5.0);
identity(id);
gran5d(70.0,15.0,id,n3,0.0);

grpn5d(1);
gran5d(70.0,267.0,d,n3,0);
gran5d(70.0,260.0,a,n3,0);
gran5d(70.0,253.0,b,n3,0);
gran5d(70.0,246.0,c,n3,0);

grpn5d(2);
FOR i:=1 TO strains_n DO
   BEGIN
   y:=300.0 - 8*i;
   gran5d(205.0,y,ri,n3,0)
   END
END;

(*
---------------------------------------------------------
*)

PROCEDURE pl_results;

VAR
   a,b : array0;
   j,k : index;
   i   : index1;
   r   : REAL;

BEGIN
FOR i:=strains_n DOWNTO 1 DO
   BEGIN
   j:=0;
   WHILE (results_s_ratei,j=0.0)
   AND   (j < t_n) DO j:=j+1;
   IF j<t_n THEN
      BEGIN
      FOR k:=0 TO j DO
         BEGIN
         bk:=t_arrayj;
         ak:=max_sr
         END;
      REPEAT
         BEGIN
         bj:=t_arrayj;
         aj:=results_s_ratei,j;
         r:=results_s_ratei,j/isothermal_s_ratej;
         j:=j+1
         END
      UNTIL ((r < acc_sr) AND (r > 1/acc_sr))
      OR    (j=t_n);
      IF j<t_n THEN FOR k:=j TO t_n DO
         BEGIN
         ak:=isothermal_s_ratek;
         bk:=t_arrayk
         END;
      grgr7d(b,a,t_n+1)
      END
   END
END;

