# build an executable for UsbMX

CC = g++
CFLAGS = -g -Wall -I

UMX = UsbMX
TARGET = Control

all: $(TARGET)

$(TARGET): $(TARGET).o $(UMX).o
	$(CC) $(CFLAGS) $< $(UMX).o -larmadillo -o $@
		
$(TARGET).o: $(TARGET).cpp
	$(CC) $(CFLAGS) -c $< -o $@
	
$(UMX).o: $(UMX).cpp $(UMX).h
	$(CC) $(CFLAGS) -c $< -o $@


target: $(TARGET)

$(TARGET): $(TARGET).o $(UMX).o
	$(CC) $(CFLAGS) $< $(UMX).o -larmadillo -o $@

clean:
	$(RM) -f core *.o $(TARGET)

cleanall:
		$(RM) -f core *.o $*.o $(TARGET) *.cpp~ *.cpp~ *.h~ 
