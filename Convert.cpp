#include <iostream>
#include <armadillo>
#include <math.h>
#include "UsbMX.h"
#include <sys/time.h>
#include <stdio.h>

using namespace std;
using namespace arma;

////////////////g++ Convert.cpp -larmadillo -o Convert///////////////////

int main(void)
{
    vec serv(3,1),
        init(3,1);
    int i, conv=1;
        
	//cout << "(1) q->enc or (2) enc->q: " << endl;
	//cin >> conv;
	
	serv << 1633 << endr
	<< 1750 << endr
	<< 2247 << endr
	<< 2358 << endr
	<< 2917 << endr
	<< 2628 << endr;
	
	/*if (conv == 1)
	{
	    serv = (serv*651.7395);
		init << 800 << endr
			 << 795 << endr // subject to slip! 
			 << 2230 << endr
			 << 1810 << endr
			 << 1950 << endr
			 << 646 << endr;
			
		for(i = 0; i < 6; i++)
		{
		    serv(i) = int(serv(i) + init(i)); //Joint offset!!!
			if (serv(i) < 0)
			{
			    serv(i) = serv(i) + 4095;
			}
		}
	}*/
		if (conv == 1)
	    {
	        init << 800 << endr
			 << 795 << endr // subject to slip! 
			 << 2230 << endr
			 << 1810 << endr
			 << 1950 << endr
			 << 646 << endr;
	        serv = serv - init;
	        serv = (serv/651.7395);
		}
		serv.print();
		return 0;
}

