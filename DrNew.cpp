
#include <iostream>
#include <armadillo>
#include <math.h>
#include "UsbMX.h"
#include <sys/time.h>
#include <stdio.h>

#define MOVE_OFFSET 2048 
#define MULTI_OFFSET 6144

using namespace std;
using namespace arma;

/* Inputs: Obe_d, phi_d, q, Abi */
/* Outputs: qd, posErr, oriErr */

double getSgn(double x);
mat sixDofJ(cube Abi);
cube getAbe(vec q);
vec rad2enc(mat q);
cube getAbeach(vec q);

//UsbMX control;

void drControl(double *Obe, int *controlPos);

////////////////g++ DrNew.cpp -larmadillo UsbMX.cpp -o Dock///////////////////
int main()
{
    double Obe[3] = {280, 0 ,-50};
    int controlPos[6] = {0, 0, 0, 0, 0, 0}, 
    i;
    
    drControl(Obe, controlPos);
    
    cout << endl << "Main servo values: " << endl;
    for(i=0; i<6; i++)
    {
        cout << controlPos[i] << endl;
    }
    
    return 0;
}


void drControl(double *Obe, int *controlPos)
{
    struct timeval start,stop;
    gettimeofday(&start,NULL);
    
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

	int a, i, attempts = 10000, reach=0, phase = 1;
	double lambda = .1;
	    
	vec serv(6,1),
        q(6,1),
        sinq(6,1),
        cosq(6,1),
 	    qd(6,1),
 	    qdeg(6,1),
 	    init(6,1),
	    obe_d(3,1),
	    phi_d(3,1),
	    obe(3,1),
	    dobe_d(3,1),
	    ocp(3,1),
	    oce(3,1),
	    oep(3,1),
	    eps_d(3,1),		 
	    eps(3,1),
	    ep(3,1),
	    eo(3,1),
	    obc(3,1),
	    o_ct(3,1), /////////////////////////////
	    oc_d(3,1),
	    oep_b(3,1),
	    csi_d(6,1),
	    err(6,1),
	    G(6,1),
	    app(3,1), //approach offset
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
	     Rp(3,3),
	     Rt(3,3),
	     Rs(3,3),
	     Rbe_d(3,3);
	   
	     
	cube Abi(4,4,8);

 while(reach == 0)
  {	      
	   if (phase == 1) //Docking
	   {
        //This is CV//
	    /*oc_d << 35    << endr // base cam 
	         << -10   << endr // to
	         << -10   << endr;// ring
	       
	    obc << -5 << endr
	        << 10 << endr //base -> cam (EYEBALL)
	        << 5.5 << endr;
	      
	    obe_d = oc_d + obc;*/ 
	    
	    /* 
	    cout << "Input Xbe(mm): " << endl;
	    cin >> obe_d(0);
	    cout << "Input Ybe(mm): " << endl;
	    cin >> obe_d(1);
	    cout << "Input Zbe(mm): " << endl;
	    cin >> obe_d(2);
	    */
	    
	    obe_d(0) = Obe[0];
	    obe_d(1) = Obe[1];
	    obe_d(2) = Obe[2];
	    
	    app << 0 << endr
	        << 0 << endr
	        << 50 << endr;
	    obe_d += app;
	   
	  
	    // Want Z straigt down & Y away from center of Target//
	    Rbe_d << 0 << 1 << 0 << endr
	          << 1 << 0 << 0 << endr
	          << 0 << 0 << -1 << endr;
	   }
	   else if (phase == 2) //Refueling
	   {
	    	 ////////* From E-E cam *///////////// 
         /*ocp << -20 << endr// E-E cam
	         << 0   << endr // to
	         << -30 << endr;//port

	     oce << -2.5 << endr
	         << 0 << endr //cam -> E-E (EYEBALL)
	         << 7 << endr;
	        //Translate camera -> base frame//
	     oep = ocp - oce;//EE -> port
	     Abi = getAbe(q);
	     Abe = Abi.slice(7);
	     obe = Abe.submat(0,3, 2,3);
	     Rbe = Abe.submat(0,0, 2,2);
	     oep_b = Rbe * oep;//w.r.t. base frame
	     obe_d = obe + oep_b;*/
	   }
	          
	// Arbitrary Start. Will read this from servos //
	q << 1.2781 << endr // desired should be +.0027
	  << 1.4653 << endr // should be - for MATLAB
	  << .0261 << endr
	  << .8408 << endr
	  << 1.4837 << endr
	  << 3.0411 << endr;
	  
	dt = 0.3; 
	K = eye(6,6);
	K(1,1) = 2;
	K(2,2) = 2;
	K(3,3) = 2;
	
	/*obe_d << 20 <<  endr
	      << 0 << endr
	      << 5 << endr;*/
    phi_d << M_PI << endr
          << 0 << endr
          << 0 << endr;
    
	/*Rp << 1  <<    0           << 0.0            << endr
	    << 0 <<  cos(phi_d(0)) << -sin(phi_d(0)) << endr
	    << 0 <<  sin(phi_d(0)) << cos(phi_d(0))  << endr;

	Rt <<  cos(phi_d(1)) << 0.0 << sin(phi_d(1)) << endr
	   <<              0.0 << 1.0 <<             0.0 << endr
	   << -sin(phi_d(1)) << 0.0 << cos(phi_d(1)) << endr;

	Rs << cos(phi_d(2)) << -sin(phi_d(2)) << 0.0 << endr
	    << sin(phi_d(2)) <<  cos(phi_d(2)) << 0.0 << endr
	    <<           0.0 <<            0.0 << 1.0 << endr;
	    
	Rd = Rp*Rt*Rs;*/
	
	eta_d = 0.5*sqrt(Rbe_d(0,0)+Rbe_d(1,1)+Rbe_d(2,2)+1);
    eps_d << 0.5*getSgn(Rbe_d(2,1)-Rbe_d(1,2))*sqrt(Rbe_d(0,0)-Rbe_d(1,1)-Rbe_d(2,2)+1) << endr
        << 0.5*getSgn(Rbe_d(0,2)-Rbe_d(2,0))*sqrt(Rbe_d(1,1)-Rbe_d(0,0)-Rbe_d(2,2)+1) << endr
	    << 0.5*getSgn(Rbe_d(1,0)-Rbe_d(0,1))*sqrt(Rbe_d(2,2)-Rbe_d(0,0)-Rbe_d(1,1)+1) << endr;
	
	/* Vector Initialization */
	posErr.zeros(1,attempts);
	oriErr.zeros(1,attempts);
	eps.zeros(3,1);

	/* Jacobian */
	for(a=0; a<attempts; a++)
	{
		// Abi = getAbeach(q);
		Abi = getAbe(q);

		/* sub q into Abe for FK */
		Abe = Abi.slice(7);
		obe = Abe.submat(0,3,2,3);
        
		/* Convert Rbe to quaternion */
		eta = 0.5*sqrt(Abe(0,0)+Abe(1,1)+Abe(2,2)+1);
        eps << 0.5*getSgn(Abe(2,1)-Abe(1,2))*sqrt(Abe(0,0)-Abe(1,1)-Abe(2,2)+1) << endr
            << 0.5*getSgn(Abe(0,2)-Abe(2,0))*sqrt(Abe(1,1)-Abe(0,0)-Abe(2,2)+1) << endr
            << 0.5*getSgn(Abe(1,0)-Abe(0,1))*sqrt(Abe(2,2)-Abe(0,0)-Abe(1,1)+1) << endr;
        //cout << eta << endl;
        //eps.print("eps: ");
        
        /* calculate error */
        if (phase == 1) 
        {
		    ep = obe_d - obe;
		}
		else if (phase == 2)
		{
		    o_ct = trans(Rbe)*(obe_d-obe)-oce; //Check this math
		    ep = Rbe*(oce + o_ct);
		}
		eo = -(eta_d*eps) + (eta*eps_d) - cross(eps_d,eps);
		posErr(0,a) = sqrt(ep(0,0)*ep(0,0)+ep(1,0)*ep(1,0)+ep(2,0)*ep(2,0)); 
		oriErr(0,a) = sqrt(eo(0,0)*eo(0,0)+eo(1,0)*eo(1,0)+eo(2,0)*eo(2,0));

		/* determine when error is acceptable */

		if((posErr(0,a) < 10e-2) && (oriErr(0,a) < 10e-2))
		{
			//cout << "Attempts: "<< a << endl;
			q(1) = -q(1); //joint 2 reversed (flip for MATLAB sim)
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
			//qd.print("qd: ");
			
			/////////rad -> deg////////////////
			qdeg = qd * 180 / M_PI;
			posErr = posErr.submat(0,0, 0,a);
			oriErr = oriErr.submat(0,0, 0,a);
			
			/////// q -> encoder values/////
			serv = (qd*651.7395);
			init << 1987 << endr
			     << 2293 << endr // subject to slip! 
			     << 2230 << endr
			     << 1810 << endr
			     << 1950 << endr
			     << 2430 << endr;
			
			for(i = 0; i < 6; i++)
			{
			    serv(i) = int(serv(i) + init(i)); //Joint offset!!!
			    if (serv(i) < 0)
			    {
			        serv(i) = serv(i) + 4095;
			    }
			}
		
			//serv.print("Encoder Values: ");
			//cout << endl << serv(0) << endl;
			for (i=0; i<6; i++)
			{
			    controlPos[i] = serv(i);
			}
	
		    if (serv(2) < 940  || serv(2) > 3330)
		    {
		        cout << "JOINT 3 O.B.!!!" << endl;
		        break;
		        reach = 0;
		    }
		    if (serv(3) < 515 || serv(3) > 3160)
		    {
		        cout << "JOINT 4 O.B.!!!" << endl;
		        break;
		        reach = 0;
		    }
		    if (serv(4) < 600 || serv(4) > 3260)
		    {
		        cout << "JOINT 5 O.B.!!!" << endl;
		        break;
		        reach = 0;
		    }
		    if (serv(5) < 1665 || serv(5) > 3715)
		    {
		        cout << "JOINT 6 O.B.!!!" << endl;
		        break;
		        reach = 0;
		    }
		    if (reach)
		    {
			   /* control.move(2, int(serv(2)));
			    usleep(2000000);
			    control.move(3, int(serv(3)));
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
		    
		    q = q + dt*solve((Jg+eye(6,6)*lambda),K*err);
		   for(i=0; i<6; i++)
		   {
		        sinq(i) = sin(q(i));
		        cosq(i) = cos(q(i));
		        q(i) = atan2(sinq(i),cosq(i));
           }
    
     if(a == attempts)
   	 {
	    cout << "Out of Workspace for Desired Orientation" << endl;
	 }
        }
  }
   gettimeofday(&stop,NULL);
   double secs;
   secs = (double)(stop.tv_usec-start.tv_usec)/1000000 + (double)(stop.tv_sec-start.tv_sec);
   //printf("time: %f\n",secs);
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
cube getAbe(vec q)
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
	double thet,D,A,a1,alph,L1,L2,L3,L4,L5,n,i,j;
    L1 = 128.8;
    L2 = 194.8;
    L3 = 118.5 + 30.0;
    L4 = 118.5;
    L5 = 261.5;
    a1 = 0; // actually 1.8573cm
    
    DHTable << L1     << -M_PI+q(0)      <<   0.0 << M_PI/2  << endr
		    << a1     << q(1) + M_PI/2   <<   L2  << -M_PI/2 << endr 
		    <<    0   << q(2)            <<   L3  << M_PI/2  << endr
		    <<    0   << q(3)            <<   L4  << 0       << endr
		    <<    0   << q(4) + M_PI/2.0 <<   0   << M_PI/2  << endr 
		    <<   L5   << q(5)            <<   0   << 0       << endr;

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
	Abi.slice(0)  = eye(4,4);
	A6e = eye(4,4);

	for (j = 0; j < 6; j++)
	{
		Abi.slice(j+1) = Abi.slice(j) * Aim1_i.slice(j);
	}
	Abi.slice(n+1) = Abi.slice(n) * A6e;
	return Abi;
}
//////////////////////////////////////////////////////////////////////////////////
cube getAbeach(vec q)
{
// Jonathan Buchholz' auto function
// Sparta Robotics FK function of q
// Last updated: 4/23/2018

// Joint position vector
float  q1,  q2,  q3,  q4,  q5,  q6;
cube Abeach_q(4,4,8);
q1 = q(0);
q2 = q(1);
q3 = q(2);
q4 = q(3);
q5 = q(4);
q6 = q(5);

// Abeach_q: transformation from base to each frame
Abeach_q.zeros(4,4,8);

Abeach_q.slice(0)  << 1.0 << 0 << 0 << 0 << endr
 << 0 << 1.0 << 0 << 0 << endr
 << 0 << 0 << 1.0 << 0 << endr
 << 0 << 0 << 0 << 1.0 << endr;

Abeach_q.slice(1)  << -1.0*cos(q1) << 0 << -1.0*sin(q1) << 0 << endr
 << -1.0*sin(q1) << 0 << cos(q1) << 0 << endr
 << 0 << 1.0 << 0 << 12.880000000000000004440892098501 << endr
 << 0 << 0 << 0 << 1.0 << endr;

Abeach_q.slice(2)  << cos(q1)*sin(q2) << sin(q1) << cos(q1)*cos(q2) << 19.479999999999999982236431605997*cos(q1)*sin(q2) << endr
 << sin(q1)*sin(q2) << -1.0*cos(q1) << cos(q2)*sin(q1) << 19.479999999999999982236431605997*sin(q1)*sin(q2) << endr
 << cos(q2) << 0 << -1.0*sin(q2) << 19.479999999999999982236431605997*cos(q2) + 12.880000000000000004440892098501 << endr
 << 0 << 0 << 0 << 1.0 << endr;

Abeach_q.slice(3)  << sin(q1)*sin(q3) + cos(q1)*cos(q3)*sin(q2) << cos(q1)*cos(q2) << cos(q1)*sin(q2)*sin(q3) - 1.0*cos(q3)*sin(q1) << 19.479999999999999982236431605997*cos(q1)*sin(q2) + 11.850000000000000005551115123126*sin(q1)*sin(q3) + 11.850000000000000005551115123126*cos(q1)*cos(q3)*sin(q2) << endr
 << cos(q3)*sin(q1)*sin(q2) - 1.0*cos(q1)*sin(q3) << cos(q2)*sin(q1) << cos(q1)*cos(q3) + sin(q1)*sin(q2)*sin(q3) << 19.479999999999999982236431605997*sin(q1)*sin(q2) - 11.850000000000000005551115123126*cos(q1)*sin(q3) + 11.850000000000000005551115123126*cos(q3)*sin(q1)*sin(q2) << endr
 << cos(q2)*cos(q3) << -1.0*sin(q2) << cos(q2)*sin(q3) << 19.479999999999999982236431605997*cos(q2) + 11.850000000000000005551115123126*cos(q2)*cos(q3) + 12.880000000000000004440892098501 << endr
 << 0 << 0 << 0 << 1.0 << endr;

Abeach_q.slice(4)  << cos(q4)*(sin(q1)*sin(q3) + cos(q1)*cos(q3)*sin(q2)) + cos(q1)*cos(q2)*sin(q4) << cos(q1)*cos(q2)*cos(q4) - 1.0*sin(q4)*(sin(q1)*sin(q3) + cos(q1)*cos(q3)*sin(q2)) << cos(q1)*sin(q2)*sin(q3) - 1.0*cos(q3)*sin(q1) << 19.479999999999999982236431605997*cos(q1)*sin(q2) + 11.850000000000000005551115123126*sin(q1)*sin(q3) + 11.850000000000000005551115123126*cos(q4)*(sin(q1)*sin(q3) + cos(q1)*cos(q3)*sin(q2)) + 11.850000000000000005551115123126*cos(q1)*cos(q3)*sin(q2) + 11.850000000000000005551115123126*cos(q1)*cos(q2)*sin(q4) << endr
 << cos(q2)*sin(q1)*sin(q4) - 1.0*cos(q4)*(cos(q1)*sin(q3) - 1.0*cos(q3)*sin(q1)*sin(q2)) << sin(q4)*(cos(q1)*sin(q3) - 1.0*cos(q3)*sin(q1)*sin(q2)) + cos(q2)*cos(q4)*sin(q1) << cos(q1)*cos(q3) + sin(q1)*sin(q2)*sin(q3) << 19.479999999999999982236431605997*sin(q1)*sin(q2) - 11.850000000000000005551115123126*cos(q1)*sin(q3) - 11.850000000000000005551115123126*cos(q4)*(cos(q1)*sin(q3) - 1.0*cos(q3)*sin(q1)*sin(q2)) + 11.850000000000000005551115123126*cos(q3)*sin(q1)*sin(q2) + 11.850000000000000005551115123126*cos(q2)*sin(q1)*sin(q4) << endr
 << cos(q2)*cos(q3)*cos(q4) - 1.0*sin(q2)*sin(q4) << - 1.0*cos(q4)*sin(q2) - 1.0*cos(q2)*cos(q3)*sin(q4) << cos(q2)*sin(q3) << 19.479999999999999982236431605997*cos(q2) + 11.850000000000000005551115123126*cos(q2)*cos(q3) - 11.850000000000000005551115123126*sin(q2)*sin(q4) + 11.850000000000000005551115123126*cos(q2)*cos(q3)*cos(q4) + 12.880000000000000004440892098501 << endr
 << 0 << 0 << 0 << 1.0 << endr;

Abeach_q.slice(5)  << - 1.0*sin(q5)*(cos(q4)*(sin(q1)*sin(q3) + cos(q1)*cos(q3)*sin(q2)) + cos(q1)*cos(q2)*sin(q4)) - 1.0*cos(q5)*(sin(q4)*(sin(q1)*sin(q3) + cos(q1)*cos(q3)*sin(q2)) - 1.0*cos(q1)*cos(q2)*cos(q4)) << cos(q1)*sin(q2)*sin(q3) - 1.0*cos(q3)*sin(q1) << cos(q5)*(cos(q4)*(sin(q1)*sin(q3) + cos(q1)*cos(q3)*sin(q2)) + cos(q1)*cos(q2)*sin(q4)) - 1.0*sin(q5)*(sin(q4)*(sin(q1)*sin(q3) + cos(q1)*cos(q3)*sin(q2)) - 1.0*cos(q1)*cos(q2)*cos(q4)) << 19.479999999999999982236431605997*cos(q1)*sin(q2) + 11.850000000000000005551115123126*sin(q1)*sin(q3) + 11.850000000000000005551115123126*cos(q4)*(sin(q1)*sin(q3) + cos(q1)*cos(q3)*sin(q2)) + 11.850000000000000005551115123126*cos(q1)*cos(q3)*sin(q2) + 11.850000000000000005551115123126*cos(q1)*cos(q2)*sin(q4) << endr
 << cos(q5)*(sin(q4)*(cos(q1)*sin(q3) - 1.0*cos(q3)*sin(q1)*sin(q2)) + cos(q2)*cos(q4)*sin(q1)) + sin(q5)*(cos(q4)*(cos(q1)*sin(q3) - 1.0*cos(q3)*sin(q1)*sin(q2)) - 1.0*cos(q2)*sin(q1)*sin(q4)) << cos(q1)*cos(q3) + sin(q1)*sin(q2)*sin(q3) << sin(q5)*(sin(q4)*(cos(q1)*sin(q3) - 1.0*cos(q3)*sin(q1)*sin(q2)) + cos(q2)*cos(q4)*sin(q1)) - 1.0*cos(q5)*(cos(q4)*(cos(q1)*sin(q3) - 1.0*cos(q3)*sin(q1)*sin(q2)) - 1.0*cos(q2)*sin(q1)*sin(q4)) << 19.479999999999999982236431605997*sin(q1)*sin(q2) - 11.850000000000000005551115123126*cos(q1)*sin(q3) - 11.850000000000000005551115123126*cos(q4)*(cos(q1)*sin(q3) - 1.0*cos(q3)*sin(q1)*sin(q2)) + 11.850000000000000005551115123126*cos(q3)*sin(q1)*sin(q2) + 11.850000000000000005551115123126*cos(q2)*sin(q1)*sin(q4) << endr
 << sin(q5)*(sin(q2)*sin(q4) - 1.0*cos(q2)*cos(q3)*cos(q4)) - 1.0*cos(q5)*(cos(q4)*sin(q2) + cos(q2)*cos(q3)*sin(q4)) << cos(q2)*sin(q3) << - 1.0*sin(q5)*(cos(q4)*sin(q2) + cos(q2)*cos(q3)*sin(q4)) - 1.0*cos(q5)*(sin(q2)*sin(q4) - 1.0*cos(q2)*cos(q3)*cos(q4)) << 19.479999999999999982236431605997*cos(q2) + 11.850000000000000005551115123126*cos(q2)*cos(q3) - 11.850000000000000005551115123126*sin(q2)*sin(q4) + 11.850000000000000005551115123126*cos(q2)*cos(q3)*cos(q4) + 12.880000000000000004440892098501 << endr
 << 0 << 0 << 0 << 1.0 << endr;

Abeach_q.slice(6)  << - 1.0*cos(q6)*(sin(q5)*(cos(q4)*(sin(q1)*sin(q3) + cos(q1)*cos(q3)*sin(q2)) + cos(q1)*cos(q2)*sin(q4)) + cos(q5)*(sin(q4)*(sin(q1)*sin(q3) + cos(q1)*cos(q3)*sin(q2)) - 1.0*cos(q1)*cos(q2)*cos(q4))) - 1.0*sin(q6)*(cos(q3)*sin(q1) - 1.0*cos(q1)*sin(q2)*sin(q3)) << sin(q6)*(sin(q5)*(cos(q4)*(sin(q1)*sin(q3) + cos(q1)*cos(q3)*sin(q2)) + cos(q1)*cos(q2)*sin(q4)) + cos(q5)*(sin(q4)*(sin(q1)*sin(q3) + cos(q1)*cos(q3)*sin(q2)) - 1.0*cos(q1)*cos(q2)*cos(q4))) - 1.0*cos(q6)*(cos(q3)*sin(q1) - 1.0*cos(q1)*sin(q2)*sin(q3)) << cos(q5)*(cos(q4)*(sin(q1)*sin(q3) + cos(q1)*cos(q3)*sin(q2)) + cos(q1)*cos(q2)*sin(q4)) - 1.0*sin(q5)*(sin(q4)*(sin(q1)*sin(q3) + cos(q1)*cos(q3)*sin(q2)) - 1.0*cos(q1)*cos(q2)*cos(q4)) << 19.479999999999999982236431605997*cos(q1)*sin(q2) + 11.850000000000000005551115123126*sin(q1)*sin(q3) + 26.150000000000000022204460492503*cos(q5)*(cos(q4)*(sin(q1)*sin(q3) + cos(q1)*cos(q3)*sin(q2)) + cos(q1)*cos(q2)*sin(q4)) - 26.150000000000000022204460492503*sin(q5)*(sin(q4)*(sin(q1)*sin(q3) + cos(q1)*cos(q3)*sin(q2)) - 1.0*cos(q1)*cos(q2)*cos(q4)) + 11.850000000000000005551115123126*cos(q4)*(sin(q1)*sin(q3) + cos(q1)*cos(q3)*sin(q2)) + 11.850000000000000005551115123126*cos(q1)*cos(q3)*sin(q2) + 11.850000000000000005551115123126*cos(q1)*cos(q2)*sin(q4) << endr
 << sin(q6)*(cos(q1)*cos(q3) + sin(q1)*sin(q2)*sin(q3)) + cos(q6)*(cos(q5)*(sin(q4)*(cos(q1)*sin(q3) - 1.0*cos(q3)*sin(q1)*sin(q2)) + cos(q2)*cos(q4)*sin(q1)) + sin(q5)*(cos(q4)*(cos(q1)*sin(q3) - 1.0*cos(q3)*sin(q1)*sin(q2)) - 1.0*cos(q2)*sin(q1)*sin(q4))) << cos(q6)*(cos(q1)*cos(q3) + sin(q1)*sin(q2)*sin(q3)) - 1.0*sin(q6)*(cos(q5)*(sin(q4)*(cos(q1)*sin(q3) - 1.0*cos(q3)*sin(q1)*sin(q2)) + cos(q2)*cos(q4)*sin(q1)) + sin(q5)*(cos(q4)*(cos(q1)*sin(q3) - 1.0*cos(q3)*sin(q1)*sin(q2)) - 1.0*cos(q2)*sin(q1)*sin(q4))) << sin(q5)*(sin(q4)*(cos(q1)*sin(q3) - 1.0*cos(q3)*sin(q1)*sin(q2)) + cos(q2)*cos(q4)*sin(q1)) - 1.0*cos(q5)*(cos(q4)*(cos(q1)*sin(q3) - 1.0*cos(q3)*sin(q1)*sin(q2)) - 1.0*cos(q2)*sin(q1)*sin(q4)) << 19.479999999999999982236431605997*sin(q1)*sin(q2) - 11.850000000000000005551115123126*cos(q1)*sin(q3) + 26.150000000000000022204460492503*sin(q5)*(sin(q4)*(cos(q1)*sin(q3) - 1.0*cos(q3)*sin(q1)*sin(q2)) + cos(q2)*cos(q4)*sin(q1)) - 26.150000000000000022204460492503*cos(q5)*(cos(q4)*(cos(q1)*sin(q3) - 1.0*cos(q3)*sin(q1)*sin(q2)) - 1.0*cos(q2)*sin(q1)*sin(q4)) - 11.850000000000000005551115123126*cos(q4)*(cos(q1)*sin(q3) - 1.0*cos(q3)*sin(q1)*sin(q2)) + 11.850000000000000005551115123126*cos(q3)*sin(q1)*sin(q2) + 11.850000000000000005551115123126*cos(q2)*sin(q1)*sin(q4) << endr
 << cos(q2)*sin(q3)*sin(q6) - 1.0*cos(q6)*(cos(q5)*(cos(q4)*sin(q2) + cos(q2)*cos(q3)*sin(q4)) - 1.0*sin(q5)*(sin(q2)*sin(q4) - 1.0*cos(q2)*cos(q3)*cos(q4))) << sin(q6)*(cos(q5)*(cos(q4)*sin(q2) + cos(q2)*cos(q3)*sin(q4)) - 1.0*sin(q5)*(sin(q2)*sin(q4) - 1.0*cos(q2)*cos(q3)*cos(q4))) + cos(q2)*cos(q6)*sin(q3) << - 1.0*sin(q5)*(cos(q4)*sin(q2) + cos(q2)*cos(q3)*sin(q4)) - 1.0*cos(q5)*(sin(q2)*sin(q4) - 1.0*cos(q2)*cos(q3)*cos(q4)) << 19.479999999999999982236431605997*cos(q2) + 11.850000000000000005551115123126*cos(q2)*cos(q3) - 11.850000000000000005551115123126*sin(q2)*sin(q4) - 26.150000000000000022204460492503*sin(q5)*(cos(q4)*sin(q2) + cos(q2)*cos(q3)*sin(q4)) - 26.150000000000000022204460492503*cos(q5)*(sin(q2)*sin(q4) - 1.0*cos(q2)*cos(q3)*cos(q4)) + 11.850000000000000005551115123126*cos(q2)*cos(q3)*cos(q4) + 12.880000000000000004440892098501 << endr
 << 0 << 0 << 0 << 1.0 << endr;

Abeach_q.slice(7)  << - 1.0*cos(q6)*(sin(q5)*(cos(q4)*(sin(q1)*sin(q3) + cos(q1)*cos(q3)*sin(q2)) + cos(q1)*cos(q2)*sin(q4)) + cos(q5)*(sin(q4)*(sin(q1)*sin(q3) + cos(q1)*cos(q3)*sin(q2)) - 1.0*cos(q1)*cos(q2)*cos(q4))) - 1.0*sin(q6)*(cos(q3)*sin(q1) - 1.0*cos(q1)*sin(q2)*sin(q3)) << sin(q6)*(sin(q5)*(cos(q4)*(sin(q1)*sin(q3) + cos(q1)*cos(q3)*sin(q2)) + cos(q1)*cos(q2)*sin(q4)) + cos(q5)*(sin(q4)*(sin(q1)*sin(q3) + cos(q1)*cos(q3)*sin(q2)) - 1.0*cos(q1)*cos(q2)*cos(q4))) - 1.0*cos(q6)*(cos(q3)*sin(q1) - 1.0*cos(q1)*sin(q2)*sin(q3)) << cos(q5)*(cos(q4)*(sin(q1)*sin(q3) + cos(q1)*cos(q3)*sin(q2)) + cos(q1)*cos(q2)*sin(q4)) - 1.0*sin(q5)*(sin(q4)*(sin(q1)*sin(q3) + cos(q1)*cos(q3)*sin(q2)) - 1.0*cos(q1)*cos(q2)*cos(q4)) << 19.479999999999999982236431605997*cos(q1)*sin(q2) + 11.850000000000000005551115123126*sin(q1)*sin(q3) + 26.150000000000000022204460492503*cos(q5)*(cos(q4)*(sin(q1)*sin(q3) + cos(q1)*cos(q3)*sin(q2)) + cos(q1)*cos(q2)*sin(q4)) - 26.150000000000000022204460492503*sin(q5)*(sin(q4)*(sin(q1)*sin(q3) + cos(q1)*cos(q3)*sin(q2)) - 1.0*cos(q1)*cos(q2)*cos(q4)) + 11.850000000000000005551115123126*cos(q4)*(sin(q1)*sin(q3) + cos(q1)*cos(q3)*sin(q2)) + 11.850000000000000005551115123126*cos(q1)*cos(q3)*sin(q2) + 11.850000000000000005551115123126*cos(q1)*cos(q2)*sin(q4) << endr
 << sin(q6)*(cos(q1)*cos(q3) + sin(q1)*sin(q2)*sin(q3)) + cos(q6)*(cos(q5)*(sin(q4)*(cos(q1)*sin(q3) - 1.0*cos(q3)*sin(q1)*sin(q2)) + cos(q2)*cos(q4)*sin(q1)) + sin(q5)*(cos(q4)*(cos(q1)*sin(q3) - 1.0*cos(q3)*sin(q1)*sin(q2)) - 1.0*cos(q2)*sin(q1)*sin(q4))) << cos(q6)*(cos(q1)*cos(q3) + sin(q1)*sin(q2)*sin(q3)) - 1.0*sin(q6)*(cos(q5)*(sin(q4)*(cos(q1)*sin(q3) - 1.0*cos(q3)*sin(q1)*sin(q2)) + cos(q2)*cos(q4)*sin(q1)) + sin(q5)*(cos(q4)*(cos(q1)*sin(q3) - 1.0*cos(q3)*sin(q1)*sin(q2)) - 1.0*cos(q2)*sin(q1)*sin(q4))) << sin(q5)*(sin(q4)*(cos(q1)*sin(q3) - 1.0*cos(q3)*sin(q1)*sin(q2)) + cos(q2)*cos(q4)*sin(q1)) - 1.0*cos(q5)*(cos(q4)*(cos(q1)*sin(q3) - 1.0*cos(q3)*sin(q1)*sin(q2)) - 1.0*cos(q2)*sin(q1)*sin(q4)) << 19.479999999999999982236431605997*sin(q1)*sin(q2) - 11.850000000000000005551115123126*cos(q1)*sin(q3) + 26.150000000000000022204460492503*sin(q5)*(sin(q4)*(cos(q1)*sin(q3) - 1.0*cos(q3)*sin(q1)*sin(q2)) + cos(q2)*cos(q4)*sin(q1)) - 26.150000000000000022204460492503*cos(q5)*(cos(q4)*(cos(q1)*sin(q3) - 1.0*cos(q3)*sin(q1)*sin(q2)) - 1.0*cos(q2)*sin(q1)*sin(q4)) - 11.850000000000000005551115123126*cos(q4)*(cos(q1)*sin(q3) - 1.0*cos(q3)*sin(q1)*sin(q2)) + 11.850000000000000005551115123126*cos(q3)*sin(q1)*sin(q2) + 11.850000000000000005551115123126*cos(q2)*sin(q1)*sin(q4) << endr
 << cos(q2)*sin(q3)*sin(q6) - 1.0*cos(q6)*(cos(q5)*(cos(q4)*sin(q2) + cos(q2)*cos(q3)*sin(q4)) - 1.0*sin(q5)*(sin(q2)*sin(q4) - 1.0*cos(q2)*cos(q3)*cos(q4))) << sin(q6)*(cos(q5)*(cos(q4)*sin(q2) + cos(q2)*cos(q3)*sin(q4)) - 1.0*sin(q5)*(sin(q2)*sin(q4) - 1.0*cos(q2)*cos(q3)*cos(q4))) + cos(q2)*cos(q6)*sin(q3) << - 1.0*sin(q5)*(cos(q4)*sin(q2) + cos(q2)*cos(q3)*sin(q4)) - 1.0*cos(q5)*(sin(q2)*sin(q4) - 1.0*cos(q2)*cos(q3)*cos(q4)) << 19.479999999999999982236431605997*cos(q2) + 11.850000000000000005551115123126*cos(q2)*cos(q3) - 11.850000000000000005551115123126*sin(q2)*sin(q4) - 26.150000000000000022204460492503*sin(q5)*(cos(q4)*sin(q2) + cos(q2)*cos(q3)*sin(q4)) - 26.150000000000000022204460492503*cos(q5)*(sin(q2)*sin(q4) - 1.0*cos(q2)*cos(q3)*cos(q4)) + 11.850000000000000005551115123126*cos(q2)*cos(q3)*cos(q4) + 12.880000000000000004440892098501 << endr
 << 0 << 0 << 0 << 1.0 << endr;

return Abeach_q;
}
/////////////////////////////////////////////////////////////////////////////////////////////
// Notes: 
// converges when q_init = [0.1, 0.1 ..... 0.1] but not zeros
//////////////////////////////////////////////////////////////////////
// Test values:
// Y from -10 -> 10cm (ideal = 0, may have to hardcode for -.7320)
// X from 23 -> 33cm (ideal = 28.8)
// Z from -6 -> -10cm
// based on assumption that CPG gets plates w/in 10cm
// measure a from catia & generate new FK
// add offset vectors for linear moves to target
//2 reversed, 1 wonky
