/*
********************************************************************************************
    Library for Dynamixel MX28AT on the Jetson TK1
    
    The purpose of this library is to simplify operation on the Jetson TK1 microcontroller
    through either the GPIO UART or USB UART. 
    
    Special thanks to Josue Alejandro Savage for the Dynamixel Arduino library for
    which this library was derived from as well as the Dynamixel SDK. In addition to
    Jetson Hacks for providing the Jetson TK1 GPIO library.
    
    MODIFICATIONS:
    2/18/2018 - Created the library with read and write functions
    2/23/2018 - Added the USB UART comptability
    
    TODO:
    - Fix Read functions to packets in order and to avoid faulty packets
    - Fix Read to read data majority of the time
    
********************************************************************************************

AUTHOR: Bruce Nelson
ORGANIZATION: Sparta Robotics

*/

#include "UsbMX.h"


void UsbMX::begin(const char *stream, speed_t baud)
{    
    uart0_filestream = open(stream, O_RDWR| O_NOCTTY );
    
    struct termios tty;
    struct termios tty_old;
    memset (&tty, 0, sizeof tty);

    if ( tcgetattr ( uart0_filestream, &tty ) != 0 ) {
       std::cout << "Error " << errno << " from tcgetattr: " << strerror(errno) << std::endl;
    }

    tty_old = tty;

    cfsetospeed (&tty, baud);
    cfsetispeed (&tty, baud);

    tty.c_cflag     &=  ~PARENB;
    tty.c_cflag     &=  ~CSTOPB;
    tty.c_cflag     &=  ~CSIZE;
    tty.c_cflag     |=  CS8;

    tty.c_cflag     &=  ~CRTSCTS;           // no flow control
    tty.c_cc[VMIN]   =  1;                  // read doesn't block
    tty.c_cc[VTIME]  =  5;                  // 0.5 seconds read timeout
    tty.c_cflag     |=  CREAD | CLOCAL;     // turn on READ & ignore ctrl lines

    cfmakeraw(&tty);

    tcflush( uart0_filestream, TCIFLUSH );
    if ( tcsetattr ( uart0_filestream, TCSANOW, &tty ) != 0) {
       std::cout << "Error " << errno << " from tcsetattr" << std::endl;
    }
}

void UsbMX::disconnect()
{	
    close(uart0_filestream);
    printf("DISCONNECTED!!! \n");
}

int UsbMX::reset(unsigned char ID)
{
    memset(tx_buffer, 0, sizeof(tx_buffer) );

	Checksum = (~(ID + MX_RESET_LENGTH + MX_RESET))&0xFF;
    
	
    tx_buffer[0] = MX_START;
    tx_buffer[1] = MX_START;
    tx_buffer[2] = ID;
    tx_buffer[3] = MX_RESET_LENGTH;
    tx_buffer[4] = MX_RESET;
    tx_buffer[5] = Checksum;
    
	count = write(uart0_filestream, &tx_buffer, 6);
	if (count < 0)
	{
		printf("UART TX error\n");
	}
	
	usleep(TX_DELAY_TIME);
	

    return 0;
}

int UsbMX::ping(unsigned char ID)
{
    memset(tx_buffer, 0, sizeof(tx_buffer) );

	Checksum = (~(ID + MX_READ_DATA + MX_PING))&0xFF;
    
	
    tx_buffer[0] = MX_START;
    tx_buffer[1] = MX_START;
    tx_buffer[2] = ID;
    tx_buffer[3] = MX_PING;
    tx_buffer[4] = Checksum;
    
	count = write(uart0_filestream, &tx_buffer, 5);
	if (count < 0)
	{
		printf("UART TX error\n");
	}
	
	usleep(TX_DELAY_TIME);
	

    return 0;
}

int UsbMX::setID(unsigned char ID, unsigned char newID)
{
    memset(tx_buffer, 0, sizeof(tx_buffer) );

	Checksum = (~(ID + MX_ID_LENGTH + MX_WRITE_DATA + MX_ID + newID))&0xFF;
    
	
    tx_buffer[0] = MX_START;
    tx_buffer[1] = MX_START;
    tx_buffer[2] = ID;
    tx_buffer[3] = MX_ID_LENGTH;
    tx_buffer[4] = MX_WRITE_DATA;
    tx_buffer[5] = MX_ID;
    tx_buffer[6] = newID;
    tx_buffer[7] = Checksum;
    
	count = write(uart0_filestream, &tx_buffer, 8);
	if (count < 0)
	{
		printf("UART TX error\n");
	}
	
	usleep(TX_DELAY_TIME);
	

    return 0;
}

int UsbMX::setBD(unsigned char ID, long baud)
{
    unsigned char Baud_Rate = (2000000/baud) - 1;
    
    memset(tx_buffer, 0, sizeof(tx_buffer) );

	Checksum = (~(ID + MX_BD_LENGTH + MX_WRITE_DATA + MX_BAUD_RATE + Baud_Rate))&0xFF;
    
	
    tx_buffer[0] = MX_START;
    tx_buffer[1] = MX_START;
    tx_buffer[2] = ID;
    tx_buffer[3] = MX_BD_LENGTH;
    tx_buffer[4] = MX_WRITE_DATA;
    tx_buffer[5] = MX_BAUD_RATE;
    tx_buffer[6] = Baud_Rate;
    tx_buffer[7] = Checksum;
    
	count = write(uart0_filestream, &tx_buffer, 8);
	if (count < 0)
	{
		printf("UART TX error\n");
	}
	
	usleep(TX_DELAY_TIME);
	

    return 0;
}

int UsbMX::move(unsigned char ID, int Position)
{
    char Position_H,Position_L;
    
    memset(tx_buffer, 0, sizeof(tx_buffer) );
   
    Position_H = Position >> 8;
    Position_L = Position;
	Checksum = (~(ID + MX_GOAL_LENGTH + MX_WRITE_DATA + MX_GOAL_POSITION_L + Position_L + Position_H))&0xFF;
    
	
    tx_buffer[0] = MX_START;
    tx_buffer[1] = MX_START;
    tx_buffer[2] = ID;
    tx_buffer[3] = MX_GOAL_LENGTH;
    tx_buffer[4] = MX_WRITE_DATA;
    tx_buffer[5] = MX_GOAL_POSITION_L;
    tx_buffer[6] = Position_L;
    tx_buffer[7] = Position_H;
    tx_buffer[8] = Checksum;
    
	count = write(uart0_filestream, &tx_buffer, 9);
	if (count < 0)
	{
		printf("UART TX error\n");
	}
	
	usleep(TX_DELAY_TIME);
	

    return 0;
}

int UsbMX::moveSpeed(unsigned char ID, int Position, int Speed)
{
    char Position_H,Position_L,Speed_H,Speed_L;
    
    memset(tx_buffer, 0, sizeof(tx_buffer) );
   
    Position_H = Position >> 8;
    Position_L = Position;
    Speed_H = Speed >> 8;
    Speed_L = Speed;
	Checksum = (~(ID + MX_GOAL_SP_LENGTH + MX_WRITE_DATA + MX_GOAL_POSITION_L + Position_L + Position_H + Speed_L + Speed_H))&0xFF;
    
	
    tx_buffer[0] = MX_START;
    tx_buffer[1] = MX_START;
    tx_buffer[2] = ID;
    tx_buffer[3] = MX_GOAL_SP_LENGTH;
    tx_buffer[4] = MX_WRITE_DATA;
    tx_buffer[5] = MX_GOAL_POSITION_L;
    tx_buffer[6] = Position_L;
    tx_buffer[7] = Position_H;
    tx_buffer[8] = Speed_L;
    tx_buffer[9] = Speed_H;
    tx_buffer[10] = Checksum;
    
	count = write(uart0_filestream, &tx_buffer, 11);
	if (count < 0)
	{
		printf("UART TX error\n");
	}
	
	usleep(TX_DELAY_TIME);
	

    return 0;
}

int UsbMX::moveDeg(unsigned char ID, int Degrees)
{
	int Position;
	
	Degrees += (CENTER);
	Position = ( float(Degrees) / 360) * 4095;
	
	move(ID, Position);
	
	return 0;
}

int UsbMX::moveSpeedDeg(unsigned char ID, int Degrees, int Speed)
{
	int Position;
	
	Degrees += CENTER;
	Position = ( float(Degrees) / 360) * 4095;
	
	moveSpeed(ID, Position, Speed);
	
	return 0;
}

int UsbMX::setEndless(unsigned char ID,bool Status)
{

    if ( Status ) {	// for continous mode
	    char MX_CCW_AL_LT = 0;
	    Checksum = (~(ID + MX_GOAL_LENGTH + MX_WRITE_DATA + MX_CCW_ANGLE_LIMIT_L))&0xFF;
    
        memset(tx_buffer, 0, sizeof(tx_buffer) );
        
	    
        tx_buffer[0] = MX_START;
        tx_buffer[1] = MX_START;
        tx_buffer[2] = ID;
        tx_buffer[3] = MX_GOAL_SP_LENGTH;
        tx_buffer[4] = MX_WRITE_DATA;
        tx_buffer[5] = MX_CW_ANGLE_LIMIT_L;
        tx_buffer[6] = MX_CCW_AL_LT;
        tx_buffer[7] = MX_CCW_AL_LT;
        tx_buffer[8] = MX_CCW_AL_LT;
        tx_buffer[9] = MX_CCW_AL_LT;
        tx_buffer[10]= Checksum;
        
	    count = write(uart0_filestream, &tx_buffer, 11);
	    if (count < 0)
	    {
		    printf("UART TX error\n");
	    }
	
	    usleep(TX_DELAY_TIME);
	    

        return 0;
    
    }
    else // for servo mode
    {
	    turn(ID,0,0);
	    Checksum = (~(ID + MX_GOAL_LENGTH + MX_WRITE_DATA + MX_CCW_ANGLE_LIMIT_L + MX_CCW_AL_L + MX_CCW_AL_H))&0xFF;
    
        memset(tx_buffer, 0, sizeof(tx_buffer) );
        
	    

        tx_buffer[0] = MX_START;
        tx_buffer[1] = MX_START;
        tx_buffer[2] = ID;
        tx_buffer[3] = MX_GOAL_LENGTH;
        tx_buffer[4] = MX_WRITE_DATA;
        tx_buffer[5] = MX_CCW_ANGLE_LIMIT_L;
        tx_buffer[6] = MX_CCW_AL_L;
        tx_buffer[7] = MX_CCW_AL_H;
        tx_buffer[8] = Checksum;
        
	    count = write(uart0_filestream, &tx_buffer, 9);
	    if (count < 0)
	    {
		    printf("UART TX error\n");
	    }
	
	    usleep(TX_DELAY_TIME);
	    

        return 0;
    }
}

int UsbMX::turn(unsigned char ID, bool SIDE, int Speed)
{
    
	if (SIDE == 0){
		char Speed_H,Speed_L;
		Speed_H = Speed >> 8;
		Speed_L = Speed;
		Checksum = (~(ID + MX_SPEED_LENGTH + MX_WRITE_DATA + MX_GOAL_SPEED_L + Speed_L + Speed_H))&0xFF;
    
	    
        tx_buffer[0] = MX_START;
        tx_buffer[1] = MX_START;
        tx_buffer[2] = ID;
        tx_buffer[3] = MX_SPEED_LENGTH;
        tx_buffer[4] = MX_WRITE_DATA;
        tx_buffer[5] = MX_GOAL_SPEED_L;
        tx_buffer[6] = Speed_L;
        tx_buffer[7] = Speed_H;
        tx_buffer[8] = Checksum;
        
	    count = write(uart0_filestream, &tx_buffer, 9);
	    if (count < 0)
	    {
		    printf("UART TX error\n");
	    }
	
	    usleep(TX_DELAY_TIME);
	    

        return 0;
	}
	else
	{
		char Speed_H,Speed_L;
		Speed_H = (Speed >> 8) + 4;
		Speed_L = Speed;
		Checksum = (~(ID + MX_SPEED_LENGTH + MX_WRITE_DATA + MX_GOAL_SPEED_L + Speed_L + Speed_H))&0xFF;
    
	    
        tx_buffer[0] = MX_START;
        tx_buffer[1] = MX_START;
        tx_buffer[2] = ID;
        tx_buffer[3] = MX_SPEED_LENGTH;
        tx_buffer[4] = MX_WRITE_DATA;
        tx_buffer[5] = MX_GOAL_SPEED_L;
        tx_buffer[6] = Speed_L;
        tx_buffer[7] = Speed_H;
        tx_buffer[8] = Checksum;
        
	    count = write(uart0_filestream, &tx_buffer, 9);
	    if (count < 0)
	    {
		    printf("UART TX error\n");
	    }
	
	    usleep(TX_DELAY_TIME);
	    

        return 0;
		}
}

int UsbMX::moveRW(unsigned char ID, int Position)
{

    char Position_H,Position_L;
    
    Position_H = Position >> 8;
    Position_L = Position;
    
    Checksum = (~(ID + MX_GOAL_LENGTH + MX_REG_WRITE + MX_GOAL_POSITION_L + Position_L + Position_H))&0xFF;
        
	
    tx_buffer[0] = MX_START;
    tx_buffer[1] = MX_START;
    tx_buffer[2] = ID;
    tx_buffer[3] = MX_GOAL_LENGTH;
    tx_buffer[4] = MX_REG_WRITE;
    tx_buffer[5] = MX_GOAL_POSITION_L;
    tx_buffer[6] = Position_L;
    tx_buffer[7] = Position_H;
    tx_buffer[8] = Checksum;
    
	count = write(uart0_filestream, &tx_buffer, 9);
	if (count < 0)
	{
		printf("UART TX error\n");
	}
	
	usleep(TX_DELAY_TIME);
	

    return 0;
}

int UsbMX::moveSpeedRW(unsigned char ID, int Position, int Speed)
{

    char Position_H,Position_L,Speed_H,Speed_L;
    
    memset(tx_buffer, 0, sizeof(tx_buffer) );
   
    Position_H = Position >> 8;
    Position_L = Position;
    Speed_H = Speed >> 8;
    Speed_L = Speed;
	Checksum = (~(ID + MX_GOAL_SP_LENGTH + MX_REG_WRITE + MX_GOAL_POSITION_L + Position_L + Position_H + Speed_L + Speed_H))&0xFF;
    
	
    tx_buffer[0] = MX_START;
    tx_buffer[1] = MX_START;
    tx_buffer[2] = ID;
    tx_buffer[3] = MX_GOAL_SP_LENGTH;
    tx_buffer[4] = MX_REG_WRITE;
    tx_buffer[5] = MX_GOAL_POSITION_L;
    tx_buffer[6] = Position_L;
    tx_buffer[7] = Position_H;
    tx_buffer[8] = Speed_L;
    tx_buffer[9] = Speed_H;
    tx_buffer[10] = Checksum;
    
	count = write(uart0_filestream, &tx_buffer, 11);
	if (count < 0)
	{
		printf("UART TX error\n");
	}
	
	usleep(TX_DELAY_TIME);
	

    return 0;
}

void UsbMX::action()
{
    
	
    tx_buffer[0] = MX_START;
    tx_buffer[1] = MX_START;
    tx_buffer[2] = BROADCAST_ID;
    tx_buffer[3] = MX_ACTION_LENGTH;
    tx_buffer[4] = MX_ACTION;
    tx_buffer[5] = MX_ACTION_CHECKSUM;
    
	count = write(uart0_filestream, &tx_buffer, 6);
	if (count < 0)
	{
		printf("UART TX error\n");
	}
	
	usleep(TX_DELAY_TIME);
	
}

int UsbMX::torqueStatus( unsigned char ID, bool Status)
{
    
    Checksum = (~(ID + MX_TORQUE_LENGTH + MX_WRITE_DATA + MX_TORQUE_ENABLE + Status))&0xFF;
    
	
    tx_buffer[0] = MX_START;
    tx_buffer[1] = MX_START;
    tx_buffer[2] = ID;
    tx_buffer[3] = MX_TORQUE_LENGTH;
    tx_buffer[4] = MX_WRITE_DATA;
    tx_buffer[5] = MX_TORQUE_ENABLE;
    tx_buffer[6] = Status;
    tx_buffer[7] = Checksum;
    
	count = write(uart0_filestream, &tx_buffer, 8);
	if (count < 0)
	{
		printf("UART TX error\n");
	}
	
	usleep(TX_DELAY_TIME);
	

    return 0;
}

int UsbMX::ledStatus( unsigned char ID, bool Status)
{
    
    Checksum = (~(ID + MX_LED_LENGTH + MX_WRITE_DATA + MX_LED + Status))&0xFF;
    
	
    tx_buffer[0] = MX_START;
    tx_buffer[1] = MX_START;
    tx_buffer[2] = ID;
    tx_buffer[3] = MX_LED_LENGTH;
    tx_buffer[4] = MX_WRITE_DATA;
    tx_buffer[5] = MX_LED;
    tx_buffer[6] = Status;
    tx_buffer[7] = Checksum;
    
	count = write(uart0_filestream, &tx_buffer, 8);
	if (count < 0)
	{
		printf("UART TX error\n");
	}
	
	usleep(TX_DELAY_TIME);
	

    return 0;
}

int UsbMX::setTempLimit(unsigned char ID, unsigned char Temperature)
{
    
    Checksum = (~(ID + MX_TL_LENGTH +MX_WRITE_DATA+ MX_LIMIT_TEMPERATURE + Temperature))&0xFF;
    
	
    tx_buffer[0] = MX_START;
    tx_buffer[1] = MX_START;
    tx_buffer[2] = ID;
    tx_buffer[3] = MX_TL_LENGTH;
    tx_buffer[4] = MX_WRITE_DATA;
    tx_buffer[5] = MX_LIMIT_TEMPERATURE;
    tx_buffer[6] = Temperature;
    tx_buffer[7] = Checksum;
    
	count = write(uart0_filestream, &tx_buffer, 8);
	if (count < 0)
	{
		printf("UART TX error\n");
	}
	
	usleep(TX_DELAY_TIME);
	

    return 0;
}

int UsbMX::setVoltageLimit(unsigned char ID, unsigned char DVoltage, unsigned char UVoltage)
{
    
    Checksum = (~(ID + MX_VL_LENGTH +MX_WRITE_DATA+ MX_DOWN_LIMIT_VOLTAGE + DVoltage + UVoltage))&0xFF;
    
	
    tx_buffer[0] = MX_START;
    tx_buffer[1] = MX_START;
    tx_buffer[2] = ID;
    tx_buffer[3] = MX_VL_LENGTH;
    tx_buffer[4] = MX_WRITE_DATA;
    tx_buffer[5] = MX_DOWN_LIMIT_VOLTAGE;
    tx_buffer[6] = DVoltage;
	tx_buffer[7] = UVoltage;
    tx_buffer[8] = Checksum;
    
	count = write(uart0_filestream, &tx_buffer, 9);
	if (count < 0)
	{
		printf("UART TX error\n");
	}
	
	usleep(TX_DELAY_TIME);
	

    return 0;
}

int UsbMX::setAngleLimit(unsigned char ID, int CWLimit, int CCWLimit)
{
    char CW_H,CW_L,CCW_H,CCW_L;
    
    CW_H = CWLimit >> 8;    
    CW_L = CWLimit;
    CCW_H = CCWLimit >> 8;
    CCW_L = CCWLimit;  
    
	Checksum = (~(ID + MX_VL_LENGTH +MX_WRITE_DATA+ MX_CW_ANGLE_LIMIT_L + CW_H + CW_L + MX_CCW_ANGLE_LIMIT_L + CCW_H + CCW_L))&0xFF;
    
	
    tx_buffer[0] = MX_START;
    tx_buffer[1] = MX_START;
    tx_buffer[2] = ID;
    tx_buffer[3] = MX_CCW_CW_LENGTH;
    tx_buffer[4] = MX_WRITE_DATA;
    tx_buffer[5] = MX_CW_ANGLE_LIMIT_L;
    tx_buffer[6] = CW_L;
    tx_buffer[7] = CW_H;
    tx_buffer[8] = MX_CCW_ANGLE_LIMIT_L;
    tx_buffer[9] = CCW_L;
    tx_buffer[10] = CCW_H;
    tx_buffer[11] = Checksum;
    
	count = write(uart0_filestream, &tx_buffer, 12);
	if (count < 0)
	{
		printf("UART TX error\n");
	}
	
	usleep(TX_DELAY_TIME);
	

    return 0;
}

int UsbMX::setMaxTorque(unsigned char ID, int MaxTorque)
{
    char MaxTorque_H,MaxTorque_L;
    
    MaxTorque_H = MaxTorque >> 8;
    MaxTorque_L = MaxTorque;
    
	Checksum = (~(ID + MX_MT_LENGTH + MX_WRITE_DATA + MX_MAX_TORQUE_L + MaxTorque_L + MaxTorque_H))&0xFF;
    
	
    tx_buffer[0] = MX_START;
    tx_buffer[1] = MX_START;
    tx_buffer[2] = ID;
    tx_buffer[3] = MX_MT_LENGTH;
    tx_buffer[4] = MX_WRITE_DATA;
    tx_buffer[5] = MX_MAX_TORQUE_L;
    tx_buffer[6] = MaxTorque_L;
    tx_buffer[7] = MaxTorque_H;
    tx_buffer[8] = Checksum;
    
	count = write(uart0_filestream, &tx_buffer, 9);
	if (count < 0)
	{
		printf("UART TX error\n");
	}
	
	usleep(TX_DELAY_TIME);
	

    return 0;
}

int UsbMX::setSRL(unsigned char ID, unsigned char SRL)
{
    
    Checksum = (~(ID + MX_SRL_LENGTH + MX_WRITE_DATA + MX_RETURN_LEVEL + SRL))&0xFF;
    
	
    tx_buffer[0] = MX_START;
    tx_buffer[1] = MX_START;
    tx_buffer[2] = ID;
    tx_buffer[3] = MX_SRL_LENGTH;
    tx_buffer[4] = MX_WRITE_DATA;
    tx_buffer[5] = MX_RETURN_LEVEL;
    tx_buffer[6] = SRL;
    tx_buffer[7] = Checksum;
    
	count = write(uart0_filestream, &tx_buffer, 8);
	if (count < 0)
	{
		printf("UART TX error\n");
	}
	
	usleep(TX_DELAY_TIME);
	

    return 0;
}

int UsbMX::setRDT(unsigned char ID, unsigned char RDT)
{
    
   Checksum = (~(ID + MX_RDT_LENGTH + MX_WRITE_DATA + MX_RETURN_DELAY_TIME + (RDT/2)))&0xFF;
    
	
    tx_buffer[0] = MX_START;
    tx_buffer[1] = MX_START;
    tx_buffer[2] = ID;
    tx_buffer[3] = MX_RDT_LENGTH;
    tx_buffer[4] = MX_WRITE_DATA;
    tx_buffer[5] = MX_RETURN_DELAY_TIME;
    tx_buffer[6] = RDT/2;
    tx_buffer[7] = Checksum;
    
	count = write(uart0_filestream, &tx_buffer, 8);
	if (count < 0)
	{
		printf("UART TX error\n");
	}
	
	usleep(TX_DELAY_TIME);
	

    return 0;
}

int UsbMX::setLEDAlarm(unsigned char ID, unsigned char LEDAlarm)
{
    
    Checksum = (~(ID + MX_LEDALARM_LENGTH + MX_WRITE_DATA + MX_ALARM_LED + LEDAlarm))&0xFF;
    
	
    tx_buffer[0] = MX_START;
    tx_buffer[1] = MX_START;
    tx_buffer[2] = ID;
    tx_buffer[3] = MX_LEDALARM_LENGTH;
    tx_buffer[4] = MX_WRITE_DATA;
    tx_buffer[5] = MX_ALARM_LED;
    tx_buffer[6] = LEDAlarm;
    tx_buffer[7] = Checksum;
    
	count = write(uart0_filestream, &tx_buffer, 8);
	if (count < 0)
	{
		printf("UART TX error\n");
	}
	
	usleep(TX_DELAY_TIME);
	

    return 0;
}

int UsbMX::setShutdownAlarm(unsigned char ID, unsigned char SALARM)
{
    
    Checksum = (~(ID + MX_SALARM_LENGTH + MX_ALARM_SHUTDOWN + MX_ALARM_LED + SALARM))&0xFF;
    
	
    tx_buffer[0] = MX_START;
    tx_buffer[1] = MX_START;
    tx_buffer[2] = ID;
    tx_buffer[3] = MX_SALARM_LENGTH;
    tx_buffer[4] = MX_WRITE_DATA;
    tx_buffer[5] = MX_ALARM_SHUTDOWN;
    tx_buffer[6] = SALARM;
    tx_buffer[7] = Checksum;
    
	count = write(uart0_filestream, &tx_buffer, 8);
	if (count < 0)
	{
		printf("UART TX error\n");
	}
	
	usleep(TX_DELAY_TIME);
	

    return 0;
}

int UsbMX::setCMargin(unsigned char ID, unsigned char CWCMargin, unsigned char CCWCMargin)
{
    
    Checksum = (~(ID + MX_CM_LENGTH +MX_WRITE_DATA+ MX_CW_COMPLIANCE_MARGIN + CWCMargin + MX_CCW_COMPLIANCE_MARGIN + CCWCMargin))&0xFF;
    
	
    tx_buffer[0] = MX_START;
    tx_buffer[1] = MX_START;
    tx_buffer[2] = ID;
    tx_buffer[3] = MX_CM_LENGTH;
    tx_buffer[4] = MX_WRITE_DATA;
    tx_buffer[5] = MX_CW_COMPLIANCE_MARGIN;
    tx_buffer[6] = CWCMargin;
    tx_buffer[7] = MX_CCW_COMPLIANCE_MARGIN;
    tx_buffer[8] = CCWCMargin;
    tx_buffer[9] = Checksum;
    
	count = write(uart0_filestream, &tx_buffer, 10);
	if (count < 0)
	{
		printf("UART TX error\n");
	}
	
	usleep(TX_DELAY_TIME);
	

    return 0;
}

int UsbMX::setCSlope(unsigned char ID, unsigned char CWCSlope, unsigned char CCWCSlope)
{
    
    Checksum = (~(ID + MX_CS_LENGTH +MX_WRITE_DATA+ MX_CW_COMPLIANCE_SLOPE + CWCSlope + MX_CCW_COMPLIANCE_SLOPE + CCWCSlope))&0xFF;
    
	
    tx_buffer[0] = MX_START;
    tx_buffer[1] = MX_START;
    tx_buffer[2] = ID;
    tx_buffer[3] = MX_CS_LENGTH;
    tx_buffer[4] = MX_WRITE_DATA;
    tx_buffer[5] = MX_CW_COMPLIANCE_SLOPE;
    tx_buffer[6] = CWCSlope;
    tx_buffer[7] = MX_CCW_COMPLIANCE_SLOPE;
    tx_buffer[8] = CCWCSlope;
    tx_buffer[9] = Checksum;
    
	count = write(uart0_filestream, &tx_buffer, 10);
	if (count < 0)
	{
		printf("UART TX error\n");
	}
	
	usleep(TX_DELAY_TIME);
	

    return 0;
}

int UsbMX::setPunch(unsigned char ID, int Punch)
{

    char Punch_H,Punch_L;
    
    Punch_H = Punch >> 8;
    Punch_L = Punch;
    
	Checksum = (~(ID + MX_PUNCH_LENGTH + MX_WRITE_DATA + MX_PUNCH_L + Punch_L + Punch_H))&0xFF;
    
	
    tx_buffer[0] = MX_START;
    tx_buffer[1] = MX_START;
    tx_buffer[2] = ID;
    tx_buffer[3] = MX_PUNCH_LENGTH;
    tx_buffer[4] = MX_WRITE_DATA;
    tx_buffer[5] = MX_PUNCH_L;
    tx_buffer[6] = Punch_L;
    tx_buffer[7] = Punch_H;
    tx_buffer[8] = Checksum;
    
	count = write(uart0_filestream, &tx_buffer, 9);
	if (count < 0)
	{
		printf("UART TX error\n");
	}
	
	usleep(TX_DELAY_TIME);
	

    return 0;
}

int UsbMX::moving(unsigned char ID)
{

    Checksum = (~(ID + MX_MOVING_LENGTH  + MX_READ_DATA + MX_MOVING + MX_BYTE_READ))&0xFF;
    
	
    tx_buffer[0] = MX_START;
    tx_buffer[1] = MX_START;
    tx_buffer[2] = ID;
    tx_buffer[3] = MX_MOVING_LENGTH;
    tx_buffer[4] = MX_READ_DATA;
    tx_buffer[5] = MX_MOVING;
    tx_buffer[6] = MX_BYTE_READ;
    tx_buffer[7] = Checksum;
    
	count = write(uart0_filestream, &tx_buffer, 8);
	if (count < 0)
	{
		printf("UART TX error\n");
	}
	
	usleep(TX_DELAY_TIME);
	

	Moving_Byte = -1;
	Time_Counter = 0;
	while((bytesToRead() < 7) & (Time_Counter < TIME_OUT))
	{
	    Time_Counter++;
	    usleep(1000);
	}
	
	Read_Byte = 0;
	Read_Byte = read(uart0_filestream, &rx_buffer, bytesToRead() );
#if 0
    printf("STUFF(%d): ", ret);
    for(int i = 0; i < ret; i++)
        printf("%d ", rx_buffer[i]);
    printf("\n");
#endif	
	if(Read_Byte < 0)
	{
        printf("ERROR! NOTHING READ!\n");
        return -1;
    }
    else
    {
        for(int iter = 0; iter < (Read_Byte-1); iter++)
        {
            if( (rx_buffer[iter] == 255) & (rx_buffer[iter+1] == 255) )
            {
                if( (Error_Byte = rx_buffer[iter+4]) != 0 )
                {
                    printf("ERROR!\n");
                    return (Error_Byte * (-1));
                }
                    
                Moving_Byte = rx_buffer[iter+5];
                return Moving_Byte;  
            }
        }
    }
	    
    return Moving_Byte;
}

int UsbMX::lockRegister(unsigned char ID)
{
    
    Checksum = (~(ID + MX_LR_LENGTH + MX_WRITE_DATA + MX_LOCK + LOCK))&0xFF;

	
    tx_buffer[0] = MX_START;
    tx_buffer[1] = MX_START;
    tx_buffer[2] = ID;
    tx_buffer[3] = MX_LR_LENGTH;
    tx_buffer[4] = MX_WRITE_DATA;
    tx_buffer[5] = MX_LOCK;
    tx_buffer[6] = LOCK;
    tx_buffer[7] = Checksum;
    
	count = write(uart0_filestream, &tx_buffer, 8);
	if (count < 0)
	{
		printf("UART TX error\n");
	}
	
	usleep(TX_DELAY_TIME);
	

    return 0;                 // Return the read error
}

int UsbMX::RWStatus(unsigned char ID)
{

    Checksum = (~(ID + MX_RWS_LENGTH  + MX_READ_DATA + MX_REGISTERED_INSTRUCTION + MX_BYTE_READ))&0xFF;
    
	
    tx_buffer[0] = MX_START;
    tx_buffer[1] = MX_START;
    tx_buffer[2] = ID;
    tx_buffer[3] = MX_RWS_LENGTH;
    tx_buffer[4] = MX_READ_DATA;
    tx_buffer[5] = MX_REGISTERED_INSTRUCTION;
    tx_buffer[6] = MX_BYTE_READ;
    tx_buffer[7] = Checksum;
    
	count = write(uart0_filestream, &tx_buffer, 8);
	if (count < 0)
	{
		printf("UART TX error\n");
	}
	
	usleep(TX_DELAY_TIME);
	

	RWS_Byte = -1;
	Time_Counter = 0;
	while((bytesToRead() < 7) & (Time_Counter < TIME_OUT))
	{
	    Time_Counter++;
	    usleep(1000);
	}
	
	Read_Byte = 0;
	Read_Byte = read(uart0_filestream, &rx_buffer, bytesToRead() );
#if 0
    printf("STUFF(%d): ", ret);
    for(int i = 0; i < ret; i++)
        printf("%d ", rx_buffer[i]);
    printf("\n");
#endif	
	if(Read_Byte < 0)
	{
        printf("ERROR! NOTHING READ!\n");
        return -1;
    }
    else
    {
        for(int iter = 0; iter < (Read_Byte-1); iter++)
        {
            if( (rx_buffer[iter] == 255) & (rx_buffer[iter+1] == 255) )
            {
                if( (Error_Byte = rx_buffer[iter+4]) != 0 )
                {
                    printf("ERROR!\n");
                    return (Error_Byte * (-1));
                }
                    
                RWS_Byte = rx_buffer[iter+5];
                return RWS_Byte;  
            }
        }
    }
	    
    return RWS_Byte;
}

int UsbMX::readTemperature(unsigned char ID)
{

    memset(tx_buffer, 0, sizeof(tx_buffer) );
    memset(tx_buffer, 0, sizeof(rx_buffer) );
   
    Checksum = (~(ID + MX_TEM_LENGTH  + MX_READ_DATA + MX_PRESENT_TEMPERATURE + MX_BYTE_READ))&0xFF;
    
	
    tx_buffer[0] = MX_START;
    tx_buffer[1] = MX_START;
    tx_buffer[2] = ID;
    tx_buffer[3] = MX_TEM_LENGTH;
    tx_buffer[4] = MX_READ_DATA;
    tx_buffer[5] = MX_PRESENT_TEMPERATURE;
    tx_buffer[6] = MX_BYTE_READ;
    tx_buffer[7] = Checksum;
    
	count = write(uart0_filestream, &tx_buffer, 8);
	if (count < 0)
	{
		printf("UART TX error\n");
	}
	
	usleep(TX_DELAY_TIME);
	

    
	Temperature_Byte = -1;
	Time_Counter = 0;
	while((bytesToRead() < 7) & (Time_Counter < TIME_OUT))
	{
	    Time_Counter++;
	    usleep(1000);
	}
	
	Read_Byte = 0;
	Read_Byte = read(uart0_filestream, &rx_buffer, bytesToRead() );
#if 0
    printf("STUFF(%d): ", ret);
    for(int i = 0; i < ret; i++)
        printf("%d ", rx_buffer[i]);
    printf("\n");
#endif	
	if(Read_Byte < 0)
	{
        printf("ERROR! NOTHING READ!\n");
        return -1;
    }
    else
    {
        for(int iter = 0; iter < (Read_Byte-1); iter++)
        {
            if( (rx_buffer[iter] == 255) & (rx_buffer[iter+1] == 255) )
            {
                if( (Error_Byte = rx_buffer[iter+4]) != 0 )
                {
                    printf("ERROR!\n");
                    return (Error_Byte * (-1));
                }
                    
                Temperature_Byte = rx_buffer[iter+6];
                return Temperature_Byte;  
            }
        }
    }
	    
    return Temperature_Byte;
}

int UsbMX::readVoltage(unsigned char ID)
{

    memset(tx_buffer, 0, sizeof(tx_buffer) );
    memset(tx_buffer, 0, sizeof(rx_buffer) );
   
    Checksum = (~(ID + MX_VOLT_LENGTH  + MX_READ_DATA + MX_PRESENT_VOLTAGE + MX_BYTE_READ))&0xFF;
    
	
    tx_buffer[0] = MX_START;
    tx_buffer[1] = MX_START;
    tx_buffer[2] = ID;
    tx_buffer[3] = MX_VOLT_LENGTH;
    tx_buffer[4] = MX_READ_DATA;
    tx_buffer[5] = MX_PRESENT_VOLTAGE;
    tx_buffer[6] = MX_BYTE_READ;
    tx_buffer[7] = Checksum;
    
	count = write(uart0_filestream, &tx_buffer, 8);
	if (count < 0)
	{
		printf("UART TX error\n");
	}
	
	usleep(TX_DELAY_TIME);
	

    
	Voltage_Byte = -1;
	Time_Counter = 0;
	while((bytesToRead() < 7) & (Time_Counter < TIME_OUT))
	{
	    Time_Counter++;
	    usleep(1000);
	}
	
	Read_Byte = 0;
	Read_Byte = read(uart0_filestream, &rx_buffer, bytesToRead() );
#if 0
    printf("STUFF(%d): ", ret);
    for(int i = 0; i < ret; i++)
        printf("%d ", rx_buffer[i]);
    printf("\n");
#endif	
	if(Read_Byte < 0)
	{
        printf("ERROR! NOTHING READ!\n");
        return -1;
    }
    else
    {
        for(int iter = 0; iter < (Read_Byte-1); iter++)
        {
            if( (rx_buffer[iter] == 255) & (rx_buffer[iter+1] == 255) )
            {
                if( (Error_Byte = rx_buffer[iter+4]) != 0 )
                {
                    printf("ERROR!\n");
                    return (Error_Byte * (-1));
                }
                    
                Voltage_Byte = rx_buffer[iter+6];
                return Voltage_Byte;  
            }
        }
    }
	    
    return Voltage_Byte;
}

int UsbMX::readPosition(unsigned char ID)
{

    memset(tx_buffer, 0, sizeof(tx_buffer) );
    memset(tx_buffer, 0, sizeof(rx_buffer) );
   
    Checksum = (~(ID + MX_POS_LENGTH  + MX_READ_DATA + MX_PRESENT_POSITION_L + MX_BYTE_READ_POS))&0xFF;
    
	
    tx_buffer[0] = MX_START;
    tx_buffer[1] = MX_START;
    tx_buffer[2] = ID;
    tx_buffer[3] = MX_POS_LENGTH;
    tx_buffer[4] = MX_READ_DATA;
    tx_buffer[5] = MX_PRESENT_POSITION_L;
    tx_buffer[6] = MX_BYTE_READ_POS;
    tx_buffer[7] = Checksum;
    
	count = write(uart0_filestream, &tx_buffer, 8);
	if (count < 0)
	{
		printf("UART TX error\n");
	}
	
	usleep(TX_DELAY_TIME);
	

    
	Position_Byte = -1;
	Time_Counter = 0;
	while((bytesToRead() < 7) & (Time_Counter < TIME_OUT))
	{
	    Time_Counter++;
	    usleep(1000);
	}
	
	Read_Byte = 0;
	Read_Byte = read(uart0_filestream, &rx_buffer, bytesToRead() );
#if 0
    printf("STUFF(%d): ", ret);
    for(int i = 0; i < ret; i++)
        printf("%d ", rx_buffer[i]);
    printf("\n");
#endif	
	if(Read_Byte < 0)
	{
        printf("ERROR! NOTHING READ!\n");
        return -1;
    }
    else
    {
        for(int iter = 0; iter < (Read_Byte-1); iter++)
        {
            if( (rx_buffer[iter] == 255) & (rx_buffer[iter+1] == 255) )
            {
                if( (Error_Byte = rx_buffer[iter+4]) != 0 )
                {
                    printf("ERROR!\n");
                    return (Error_Byte * (-1));
                }
                    
                Position_Byte = rx_buffer[iter+5] + (rx_buffer[iter+6] << 8);
                return Position_Byte;  
            }
        }
    }
	    
    return Position_Byte;
}

int UsbMX::readSpeed(unsigned char ID)
{

    memset(tx_buffer, 0, sizeof(tx_buffer) );
    memset(tx_buffer, 0, sizeof(rx_buffer) );
   
    Checksum = (~(ID + MX_POS_LENGTH  + MX_READ_DATA + MX_PRESENT_SPEED_L + MX_BYTE_READ_POS))&0xFF;
    
	
    tx_buffer[0] = MX_START;
    tx_buffer[1] = MX_START;
    tx_buffer[2] = ID;
    tx_buffer[3] = MX_POS_LENGTH;
    tx_buffer[4] = MX_READ_DATA;
    tx_buffer[5] = MX_PRESENT_SPEED_L;
    tx_buffer[6] = MX_BYTE_READ_POS;
    tx_buffer[7] = Checksum;
    
	count = write(uart0_filestream, &tx_buffer, 8);
	if (count < 0)
	{
		printf("UART TX error\n");
	}
	
	usleep(TX_DELAY_TIME);
	

    
	Speed_Byte = -1;
	Time_Counter = 0;
	while((bytesToRead() < 7) & (Time_Counter < TIME_OUT))
	{
	    Time_Counter++;
	    usleep(1000);
	}
	
	Read_Byte = 0;
	Read_Byte = read(uart0_filestream, &rx_buffer, bytesToRead() );
#if 0
    printf("STUFF(%d): ", ret);
    for(int i = 0; i < ret; i++)
        printf("%d ", rx_buffer[i]);
    printf("\n");
#endif	
	if(Read_Byte < 0)
	{
        printf("ERROR! NOTHING READ!\n");
        return -1;
    }
    else
    {
        for(int iter = 0; iter < (Read_Byte-1); iter++)
        {
            if( (rx_buffer[iter] == 255) & (rx_buffer[iter+1] == 255) )
            {
                if( (Error_Byte = rx_buffer[iter+4]) != 0 )
                {
                    printf("ERROR!\n");
                    return (Error_Byte * (-1));
                }
                    
                Speed_Byte = rx_buffer[iter+5] + (rx_buffer[iter+6] << 8);
                return Speed_Byte;  
            }
        }
    }
	    
    return Speed_Byte;
}

int UsbMX::readLoad(unsigned char ID)
{

    memset(tx_buffer, 0, sizeof(tx_buffer) );
    memset(tx_buffer, 0, sizeof(rx_buffer) );
   
    Checksum = (~(ID + MX_POS_LENGTH  + MX_READ_DATA + MX_PRESENT_LOAD_L + MX_BYTE_READ_POS))&0xFF;
    
	
    tx_buffer[0] = MX_START;
    tx_buffer[1] = MX_START;
    tx_buffer[2] = ID;
    tx_buffer[3] = MX_POS_LENGTH;
    tx_buffer[4] = MX_READ_DATA;
    tx_buffer[5] = MX_PRESENT_LOAD_L;
    tx_buffer[6] = MX_BYTE_READ_POS;
    tx_buffer[7] = Checksum;
    
	count = write(uart0_filestream, &tx_buffer, 8);
	if (count < 0)
	{
		printf("UART TX error\n");
	}
	
	usleep(TX_DELAY_TIME);
	

    
	Load_Byte = -1;
	Time_Counter = 0;
	while((bytesToRead() < 7) & (Time_Counter < TIME_OUT))
	{
	    Time_Counter++;
	    usleep(1000);
	}
	
	Read_Byte = 0;
	Read_Byte = read(uart0_filestream, &rx_buffer, bytesToRead() );
#if 0
    printf("STUFF(%d): ", ret);
    for(int i = 0; i < ret; i++)
        printf("%d ", rx_buffer[i]);
    printf("\n");
#endif	
	if(Read_Byte < 0)
	{
        printf("ERROR! NOTHING READ!\n");
        return -1;
    }
    else
    {
        for(int iter = 0; iter < (Read_Byte-1); iter++)
        {
            if( (rx_buffer[iter] == 255) & (rx_buffer[iter+1] == 255) )
            {
                if( (Error_Byte = rx_buffer[iter+4]) != 0 )
                {
                    printf("ERROR!\n");
                    return (Error_Byte * (-1));
                }
                    
                Load_Byte = rx_buffer[iter+5] + (rx_buffer[iter+6] << 8);
                return Load_Byte;  
            }
        }
    }
	    
    return Load_Byte;
}

int UsbMX::bytesToRead()
{
    int bytes = 0;
	ioctl(uart0_filestream, FIONREAD, &bytes);
    return bytes;
}


