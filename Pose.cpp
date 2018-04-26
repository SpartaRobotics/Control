////// Pose Test //////

#include <iostream>
#include <armadillo>
#include <math.h>
#include "UsbMX.h"

#define MOVE_OFFSET 2048 
#define MULTI_OFFSET 6144

using namespace std;
using namespace arma;

UsbMX control;

////////////////g++ Pose.cpp -larmadillo UsbMX.cpp -o Pose///////////////////

int main(void)
{
    vec serv(6,1);
    
    // Communication Setup
    control.begin("/dev/ttyUSB0", B1000000);
	control.setEndless(1, OFF); // Sets the servo to "Servo" mode
	usleep(10*1000);
	control.setEndless(2, OFF); // Sets the servo to "Servo" mode
    usleep(10*1000);
	control.setEndless(3, OFF); // Sets the servo to "Servo" mode
    usleep(10*1000);
	control.setEndless(4, OFF); // Sets the servo to "Servo" mode
    usleep(10*1000);
	control.setEndless(5, OFF); // Sets the servo to "Servo" mode
    usleep(10*1000);
	control.setEndless(6, OFF); // Sets the servo to "Servo" mode
	usleep(10*1000);
	control.move(1, 0);  
    usleep(1000);
	control.move(2, 0);
	usleep(1000);
	control.move(3, 0);
	usleep(1000);
	control.move(4, 0);
	usleep(1000);
	control.move(5, 0);
	usleep(1000);
	control.move(6, 0);
	usleep(5000*1000);
	// End COMM Setup
	
	serv << 52 << endr
	     << 3372 << endr
	     << 2759 << endr
	     << 2208 << endr
	     << 3123 << endr
	     << 2910 << endr;
	
	//Move//
	 control.moveSpeed(1, int(serv(0)+MOVE_OFFSET),16);
	 usleep(2000000);
	 control.moveSpeed(2, int(serv(1)+MOVE_OFFSET),16);
	 usleep(2000000);
	 control.moveSpeed(3, int(serv(2)+MOVE_OFFSET),16);
	 usleep(2000000);
	 control.moveSpeed(4, int(serv(3)+MOVE_OFFSET),16);
	 usleep(2000000);
	 control.moveSpeed(5, int(serv(4)+MOVE_OFFSET),16);
	 usleep(2000000);
	 control.moveSpeed(6, int(serv(5)+MOVE_OFFSET),16);
	 usleep(2000000);
	 
}
