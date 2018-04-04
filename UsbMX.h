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
    - Adjust for user input UART and baud rates
    
********************************************************************************************

AUTHOR: Bruce Nelson
ORGANIZATION: Sparta Robotics

*/


#ifndef UsbMX_h
#define UsbMX_h

	// EEPROM AREA  ///////////////////////////////////////////////////////////
#define MX_MODEL_NUMBER_L           0
#define MX_MODEL_NUMBER_H           1
#define MX_VERSION                  2
#define MX_ID                       3
#define MX_BAUD_RATE                4
#define MX_RETURN_DELAY_TIME        5
#define MX_CW_ANGLE_LIMIT_L         6
#define MX_CW_ANGLE_LIMIT_H         7
#define MX_CCW_ANGLE_LIMIT_L        8
#define MX_CCW_ANGLE_LIMIT_H        9
#define MX_SYSTEM_DATA2             10
#define MX_LIMIT_TEMPERATURE        11
#define MX_DOWN_LIMIT_VOLTAGE       12
#define MX_UP_LIMIT_VOLTAGE         13
#define MX_MAX_TORQUE_L             14
#define MX_MAX_TORQUE_H             15
#define MX_RETURN_LEVEL             16
#define MX_ALARM_LED                17
#define MX_ALARM_SHUTDOWN           18
#define MX_OPERATING_MODE           19
#define MX_DOWN_CALIBRATION_L       20
#define MX_DOWN_CALIBRATION_H       21
#define MX_UP_CALIBRATION_L         22
#define MX_UP_CALIBRATION_H         23

	// RAM AREA  //////////////////////////////////////////////////////////////
#define MX_TORQUE_ENABLE            24
#define MX_LED                      25
#define MX_CW_COMPLIANCE_MARGIN     26
#define MX_CCW_COMPLIANCE_MARGIN    27
#define MX_CW_COMPLIANCE_SLOPE      28
#define MX_CCW_COMPLIANCE_SLOPE     29
#define MX_GOAL_POSITION_L          30
#define MX_GOAL_POSITION_H          31
#define MX_GOAL_SPEED_L             32
#define MX_GOAL_SPEED_H             33
#define MX_TORQUE_LIMIT_L           34
#define MX_TORQUE_LIMIT_H           35
#define MX_PRESENT_POSITION_L       36
#define MX_PRESENT_POSITION_H       37
#define MX_PRESENT_SPEED_L          38
#define MX_PRESENT_SPEED_H          39
#define MX_PRESENT_LOAD_L           40
#define MX_PRESENT_LOAD_H           41
#define MX_PRESENT_VOLTAGE          42
#define MX_PRESENT_TEMPERATURE      43
#define MX_REGISTERED_INSTRUCTION   44
#define MX_PAUSE_TIME               45
#define MX_MOVING                   46
#define MX_LOCK                     47
#define MX_PUNCH_L                  48
#define MX_PUNCH_H                  49

    // Instruction Set ///////////////////////////////////////////////////////////////
#define MX_PING                     1
#define MX_READ_DATA                2
#define MX_WRITE_DATA               3
#define MX_REG_WRITE                4
#define MX_ACTION                   5
#define MX_RESET                    6
#define MX_SYNC_WRITE               131

	// Specials ///////////////////////////////////////////////////////////////
#define OFF                         0
#define ON                          1
#define LEFT						0
#define RIGHT                       1
#define MX_BYTE_READ                1
#define MX_BYTE_READ_POS            2
#define MX_RESET_LENGTH				2
#define MX_ACTION_LENGTH			2
#define MX_ID_LENGTH                4
#define MX_LR_LENGTH                4
#define MX_SRL_LENGTH               4
#define MX_RDT_LENGTH               4
#define MX_LEDALARM_LENGTH          4
#define MX_SALARM_LENGTH            4
#define MX_TL_LENGTH                4
#define MX_VL_LENGTH                6
#define MX_CM_LENGTH                6
#define MX_CS_LENGTH                6
#define MX_CCW_CW_LENGTH            8
#define MX_BD_LENGTH                4
#define MX_TEM_LENGTH               4
#define MX_MOVING_LENGTH            4
#define MX_RWS_LENGTH               4
#define MX_VOLT_LENGTH              4
#define MX_LED_LENGTH               4
#define MX_TORQUE_LENGTH            4
#define MX_POS_LENGTH               4
#define MX_MT_LENGTH                5
#define MX_PUNCH_LENGTH             5
#define MX_SPEED_LENGTH             5
#define MX_GOAL_LENGTH              5
#define MX_GOAL_SP_LENGTH           7
#define MX_ACTION_CHECKSUM			250
#define BROADCAST_ID                254
#define MX_START                    255
#define MX_CW_AL_L                  255 
#define MX_CW_AL_H                  15
#define MX_CCW_AL_L                 255 
#define MX_CCW_AL_H                 15
#define TIME_OUT                    10         
#define TX_DELAY_TIME				400
#define Tx_MODE                     1
#define Rx_MODE                     0
#define LOCK                        1
#define CENTER						180
#define FULL_SPEED					1020
#define TRIQUARTER_SPEED			756
#define HALF_SPEED					512
#define QUARTER_SPEED				256
#define STOP						0

#include <iostream>
#include <stdio.h>
#include <unistd.h>     // Used for UART
#include <fcntl.h>      // Used for UART
#include <termios.h>    // Used for UART
#include <inttypes.h>   // Types
#include <sys/ioctl.h>  // UART Read
#include <errno.h>
#include <string.h>

class UsbMX {
private:

    struct termios options;
    
	unsigned char tx_buffer[14];
	unsigned char rx_buffer[8];
    
	unsigned char Checksum; 
	unsigned char Direction_Pin;
	unsigned char Time_Counter;
	unsigned char Incoming_Byte;               
	unsigned char Position_High_Byte;
	unsigned char Position_Low_Byte;
	unsigned char Speed_High_Byte;
	unsigned char Speed_Low_Byte;
	unsigned char Load_High_Byte;
	unsigned char Load_Low_Byte;
	
	int uart0_filestream;
	int count;
	int Read_Byte;
	int Moving_Byte;
	int RWS_Byte;
	int Speed_Byte;
	int Load_Byte;
	int Position_Byte;
	int Temperature_Byte;
	int Voltage_Byte;
	int Error_Byte; 

public:
    void begin(const char *stream, speed_t baud);
    void disconnect();
    
    int reset(unsigned char ID);
	int ping(unsigned char ID);
	
	int setID(unsigned char ID, unsigned char newID);
	int setBD(unsigned char ID, long baud);
    
    int move(unsigned char ID, int Position);
	int moveSpeed(unsigned char ID, int Position, int Speed);
	int moveDeg(unsigned char ID, int Degrees);
	int moveSpeedDeg(unsigned char ID, int Degrees, int Speed);
	int setEndless(unsigned char ID,bool Status);
	int turn(unsigned char ID, bool SIDE, int Speed);
	int moveRW(unsigned char ID, int Position);
	int moveSpeedRW(unsigned char ID, int Position, int Speed);
	
	void action(void);
    
	int torqueStatus(unsigned char ID, bool Status);
	int ledStatus(unsigned char ID, bool Status);

	int readTemperature(unsigned char ID);
	int readVoltage(unsigned char ID);
	int readPosition(unsigned char ID);
	int readSpeed(unsigned char ID);
	int readLoad(unsigned char ID);
	
	int setTempLimit(unsigned char ID, unsigned char Temperature);
	int setAngleLimit(unsigned char ID, int CWLimit, int CCWLimit);
	int setVoltageLimit(unsigned char ID, unsigned char DVoltage, unsigned char UVoltage);
	int setMaxTorque(unsigned char ID, int MaxTorque);
	int setSRL(unsigned char ID, unsigned char SRL);
	int setRDT(unsigned char ID, unsigned char RDT);
	int setLEDAlarm(unsigned char ID, unsigned char LEDAlarm);
	int setShutdownAlarm(unsigned char ID, unsigned char SALARM);
	int setCMargin(unsigned char ID, unsigned char CWCMargin, unsigned char CCWCMargin);
	int setCSlope(unsigned char ID, unsigned char CWCSlope, unsigned char CCWCSlope);
	int setPunch(unsigned char ID, int Punch);
    
	int moving(unsigned char ID);
	int lockRegister(unsigned char ID);
	int RWStatus(unsigned char ID);
	
    int bytesToRead();
};

extern UsbMX UMX;

#endif
