
#include <iostream>
#include <armadillo>
#include <math.h>
#include "UsbMX.h"

#define MOVE_OFFSET 2048 
#define MULTI_OFFSET 6144

using namespace std;
using namespace arma;

/* Inputs: Obe_d, phi_d, q, Abi */
/* Outputs: qd, posErr, oriErr */

double getSgn(double x);
mat sixDofJ(cube Abi);
cube getAb_im1(vec q);
vec rad2enc(mat q);
vec FK(vec q);

//UsbMX control;

////////////////g++ Dock.cpp -larmadillo UsbMX.cpp -o Dock///////////////////

int main(void)
{
    // Communication Setup
    /*control.begin("/dev/ttyUSB0", B1000000);
	control.setEndless(2, OFF); // Sets the servo to "Servo" mode
	usleep(10*1000);
	control.setEndless(3, OFF); // Sets the servo to "Servo" mode
	usleep(10*1000);
	control.move(2, 0);  
    usleep(1000);
	control.move(3, 0);
	usleep(5000*1000);*/
	// End COMM Setup

	double offset, dt, eta_d, eta, qi, q2pi, mag;

	int a, i, attempts = 10000, reach=0;
	    
	vec serv(6,1),
        q(6,1),
 	    qd(6,1),
 	    qdeg(6,1),
 	    init(6,1),
	    obe_d(3,1),
	    obe(3,1),
	    eps_d(3,1),		 
	    eps(3,1),
	    ep(3,1),
	    eo(3,1),
	    obc(3,1),
	    oc_d(3,1),
	    csi_d(6,1),
	    err(6,1),
	    G(6,1),
	    fkOut(7,1); //contains obe_d, eta_d & eps_d
	    
	mat  K(6,6), 
	     Rz1(3,3),
	     Ry(3,3),
	     Rz2(3,3),
	     Rd(3,3),
	     Jg(6,6),
	     Abe(4,4),
	     posErr(1, attempts),
	     oriErr(1, attempts),
	     Rbe(3,3),
	     Rbe_d(3,3);
	     
	cube Abi(4,4,8);
	     
	 ////////* From base cam *///////////// 
	 /*obe_d << -1.4995 << endr// upper limit: 32 (when z=8)
	       << 25.4758 << endr // lower limit: based on assembly
	       << 1.7648 << endr;// lower limit: < -45 (w/10 on xy)*/
	 //////// must be reachable with Rbe_d/////////
	 
 while(reach == 0)
  {
	/* oc_d << 0 << endr // base cam 
	      << 10 << endr// to     (THIS WORKS)
	      << -5.5 << endr;// ring*/
	      
	   oc_d << 4.268 << endr // base cam 
	      << 18 << endr// to
	      << -30 << endr;// ring
	       
	 obc << -5 << endr
	     << 10 << endr //base -> cam (EYEBALL)
	     << 5.5 << endr;
	      
	  obe_d = oc_d + obc;  
	  obe_d.print("obe_d: ");
	  
	 // Want Z straigt down & Y away from center of Target//
	  Rbe_d << 1 << 0 << 0 << endr
	        << 0 << -1 << 0 << endr
	        << 0 << 0 << -1 << endr;

	eta_d = 0.5*sqrt(Rbe_d(0,0)+Rbe_d(1,1)+Rbe_d(2,2)+1);
    eps_d << 0.5*getSgn(Rbe_d(2,1)-Rbe_d(1,2))*sqrt(Rbe_d(0,0)-Rbe_d(1,1)-Rbe_d(2,2)+1) << endr
        << 0.5*getSgn(Rbe_d(0,2)-Rbe_d(2,0))*sqrt(Rbe_d(1,1)-Rbe_d(0,0)-Rbe_d(2,2)+1) << endr
	    << 0.5*getSgn(Rbe_d(1,0)-Rbe_d(0,1))*sqrt(Rbe_d(2,2)-Rbe_d(0,0)-Rbe_d(1,1)+1) << endr;
	           
	// Arbitrary Start. Will read this from servos //
	q << 0.1 << endr // desired should be +.0027
	  << 0.1 << endr
	  << 0.1 << endr
	  << 0.1 << endr
	  << 0.1 << endr
	  << 0.1 << endr;
	dt = 0.3; 
	K = eye(6,6);

	/* Vector Initialization */
	posErr.zeros(1,attempts);
	oriErr.zeros(1,attempts);
	eps.zeros(3,1);

	/* Jacobian */
	for(a=0; a<attempts; a++)
	{
		Abi = getAb_im1(q);

		/* sub q into Abe for FK */
		Abe = Abi.slice(7);
		obe = Abe.submat(0,3,2,3);

		/* Convert Rbe to quaternion */
		eta = 0.5*sqrt(Abe(0,0)+Abe(1,1)+Abe(2,2)+1);
    eps << 0.5*getSgn(Abe(2,1)-Abe(1,2))*sqrt(Abe(0,0)-Abe(1,1)-Abe(2,2)+1) << endr
        << 0.5*getSgn(Abe(0,2)-Abe(2,0))*sqrt(Abe(1,1)-Abe(0,0)-Abe(2,2)+1) << endr
	    << 0.5*getSgn(Abe(1,0)-Abe(0,1))*sqrt(Abe(2,2)-Abe(0,0)-Abe(1,1)+1) << endr;
        
		/* calculate error */
		ep = obe_d - obe;
		eo = -(eta_d*eps) + (eta*eps_d) - cross(eps_d,eps); 

		posErr(0,a) = sqrt(ep(0,0)*ep(0,0)+ep(1,0)*ep(1,0)+ep(2,0)*ep(2,0)); 
		oriErr(0,a) = sqrt(eo(0,0)*eo(0,0)+eo(1,0)*eo(1,0)+eo(2,0)*eo(2,0));

		/* determine when error is acceptable */

		if((posErr(0,a) < 10e-4) && (oriErr(0,a) < 10e-4))
		{
			cout << "Attempts: "<< a << endl;
			q.print("q: ");
			reach = 1;
			for(i=0; i<6; i++)
			{
			    // Convert (q>pi) --> (-pi<q<pi) //
			    q2pi = q(i)/(2*M_PI); 
			    qi = int(q2pi);
			    qd(i) = q2pi-qi;
			    if (qd(i) <= 0)
			    { 
			        mag = abs(qd(i));
			        if (mag > .5)
			        {
			            qd(i) = 2*M_PI*(1-mag);
			        }
			        else
		   	        {
			            qd(i) = 2*M_PI*qd(i);
			        }
			    }
			    else
			    { 
			        if (qd(i) > .5)
			        {
			            qd(i) = -2*M_PI*(1-qd(i));
			        }
			        else
		   	        {
			            qd(i) = 2*M_PI*qd(i);
			        }
			    }
			}
			qd.print("qd: ");
			/////////rad -> deg////////////////
			qdeg = qd * 180 / M_PI;
			qdeg.print("qdeg: ");
			posErr = posErr.submat(0,0, 0,a);
			oriErr = oriErr.submat(0,0, 0,a);
			
			/////// q -> encoder values/////
			serv = (qd*651.7395);
			init << 475 << endr
			     << 3700 << endr // subject to slip! 
			     << 2230 << endr
			     << 3033 << endr
			     << 1950 << endr
			     << 1900 << endr;
			

			
			for(i = 1; i < 6; i++)
			{
			    serv(i) = int(serv(i) + init(i)); //Joint offset!!!
			}
		
			if (serv(0) < 0)
			{
			    serv(0) = (qd(0)*3)+(2*M_PI);
			    serv(0) = 1425+(serv(0)*651.7395);
			    serv(0) = 4095+serv(0);
			}
			
		    if (serv(3) > 4095)
            {
                serv(3) = init(3) + serv(3) - 4095;
			}

			serv.print("Encoder Values: ");
		    
		    // Check joint limits //
		    if (serv(1) > 7173 || serv(1) < 77) //Subject to slip!
		    {
		        cout << "JOINT 2 O.B.!!!" << endl;
	            break;
	            reach = 0;
		    }
		    if (serv(2) < 920 || serv(2) > 3555)
		    {
		        cout << "JOINT 3 O.B.!!!" << endl;
		        break;
		        reach = 0;
		    }
		    if (serv(3) < 0 || serv(3) > 4095)
		    {
		        cout << "JOINT 4 O.B.!!!" << endl;
		        break;
		        reach = 0;
		    }
		    if (serv(4) < 600 || serv(4) > 3250)
		    {
		        cout << "JOINT 5 O.B.!!!" << endl;
		        break;
		        reach = 0;
		    }
		    if (serv(5) < 630 || serv(5) > 3266)
		    {
		        cout << "JOINT 6 O.B.!!!" << endl;
		        break;
		        reach = 0;
		    }
		    if (reach)
		    {
			   /* control.move(2, int(serv(4)+MOVE_OFFSET));
			    usleep(2000000);
			    control.move(3, int(serv(5)+MOVE_OFFSET));
			    usleep(2000000);*/
			}
		   
			break;
		}
                      
		/* Iterate q using Jg */
		Jg = sixDofJ(Abi); 
		for(i=0; i<3; i++)
		{
			err(i,0)   = ep(i,0);
			err(i+3,0) = eo(i,0);		
		}
		//err.print("err: ");
		q = q + dt*solve(Jg,K*err);
    }
    
     if(a == attempts)
   	 {
	    cout << "Out of Workspace for Desired Orientation" << endl;
	 }
  }
    return 0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
/* sgn Function */
double getSgn(double x)
{ 
	if(x >= 0)
		return 1;
	else
		return -1;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Jacobian Function */
mat sixDofJ(cube Abi)
{
	int j;
	vec B(3,1), Oe(3,1), Oj(3,1), E(3,1);
	mat A(4,4), C(4,4), D(4,4), Jg(6,6);
	Jg.zeros(6,6);
	for(j = 0; j < 6; j = j + 1)
	{
		A = Abi.slice(j);
		B = A.submat(0,2, 2,2);

		C = Abi.slice(7);
		Oe = C.submat(0,3, 2,3);

		D = Abi.slice(j);
		Oj = D.submat(0,3, 2,3);

		E = cross(B, (Oe-Oj));

		Jg.submat(0,j, 2,j) = E;
		Jg.submat(3,j, 5,j) = B;
	}
	return Jg;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
cube getAb_im1(vec q)
{
	/* Generates transformation matrices between joint frames of a manipulator
	% 'i' refers to index of matrix depth (:,:,i)
	%  Aim1_i: aka each A from (i-1) frame to (i) for i = 1:n
	%    A01, A12, A13, ... A(n-1)n
	%  Ab_im1:   aka each A from base frame to (i-1) for i = [1:n, e]
	%    Ab0, Ab1, Ab2, ... Abn, Abe */
	mat d(6,1),theta(6,1),a(6,1),alpha(6,1);
	mat DHTable(6,4), A6e(4,4);
	cube Aim1_i(4,4,6), Abi(4,4,8);
	double thet,D,A,alph, L0 = 3.5, L6 = 2, n, i, j;

    DHTable << 7.3076 << q(0) - M_PI/2.0 <<   0.0 << M_PI/2 << endr
		    << 3.232  << q(1) + M_PI     <<  8.96 << 0 << endr 
		    <<   0.0  << q(2) + M_PI/2.0 <<   0.0 << M_PI/2 << endr
		    << 17.36  << q(3) - M_PI/2.0 <<   0.0 << M_PI/2 << endr
		    <<   0.0  << q(4) + M_PI/2.0 << 11.96 << M_PI/2 << endr 
		    <<   0.0  << q(5)            << 18.5  << 0  << endr;

	/* Separate DH table into columns */

	d     = DHTable.submat(0,0, 5,0);
	theta = DHTable.submat(0,1, 5,1);
	a     = DHTable.submat(0,2, 5,2);
	alpha = DHTable.submat(0,3, 5,3);

	/* n: number of joints */
	n = 6;

	/* Preallocate matrices */
	Aim1_i.zeros(4,4,6);
	Abi.zeros(4,4,8);

	/* Calculate each A from (i-1) to (i) */
	for (i = 0; i < 6; i++)
	{
	 	thet = theta(i,0); //'operator=' error w/out these
		alph = alpha(i,0);
		A = a(i,0);
		D = d(i,0); 
	    Aim1_i.slice(i) << cos(thet) << -sin(thet)*cos(alph) <<  sin(thet)*sin(alph) << A*cos(thet) << endr
 			            << sin(thet) <<  cos(thet)*cos(alph) << -cos(thet)*sin(alph) << A*sin(thet) << endr
 			            << 0.0 << sin(alph) << cos(alph) << D   << endr
 			            << 0.0 << 0.0       << 0.0       << 1.0 << endr;
	}

	// Calculate each Ab(i-1)
	Abi.slice(0) << 1.0 << 0.0  << 0.0 << 2.5 << endr
		         << 0.0  << 1.0 << 0.0 << 3.5    << endr
	             << 0.0  << 0.0  << 1.0 << 2 << endr
		         << 0.0  << 0.0  << 0.0 << 1.0    << endr;
	A6e << 0.0 << 0   << 1   << 0.0 << endr
	    << 1   << 0   << 0.0 << 0.0 << endr
	    << 0   << 1   << 0   << 0.0 << endr
	    << 0.0 << 0.0 << 0.0 << 1.0 << endr;

	for (j = 0; j < 6; j++)
	{
		Abi.slice(j+1) = Abi.slice(j) * Aim1_i.slice(j);
	}
	Abi.slice(n+1) = Abi.slice(n) * A6e;
	return Abi;
}

//////////////////////////////////////////////////////////////////////////////////

// Notes: 
// converges when q_init = [0.1, 0.1 ..... 0.1] but not zeros
//////////////////////////////////////////////////////////////////////
// Test values:
// X from -10 -> 10cm (ideal = 0, may have to hardcode for -.7320)
// Y from 23 -> 33cm (ideal = 28.8)
// Z from -6 -> -10cm
// based on assumption that CPG gets plates w/in 10cm







