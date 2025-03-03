#include <SPI.h>

#define GPIOA A0
#define DAC_CS 10  
#define BUFFER_SIZE 100000  

void setup() {
    Serial.begin(115200);
    SPI.begin();
    pinMode(DAC_CS, OUTPUT);
    digitalWrite(DAC_CS, HIGH);
}

void loop() {
    static int buffer[BUFFER_SIZE];  
    int count = 0;

    while (Serial.available()) {
        String data = Serial.readStringUntil('\n');  
        if (data == "END") {
            Serial.flush();   
            outputSignal(buffer, count);  
            return;
        }

        char *ptr = strtok((char*)data.c_str(), ",");
        while (ptr != NULL && count < BUFFER_SIZE) {
            buffer[count++] = atoi(ptr);
            ptr = strtok(NULL, ",");
        }
    }
}

void outputSignal(int *buffer, int count) {
    for (int i = 0; i < count; i++) {
        writeDAC(buffer[i]);
        delayMicroseconds(100);  // Adjust for sampling rate
    }
}

void writeDAC(uint16_t value) {
    digitalWrite(DAC_CS, LOW);
    SPI.transfer16(0x3000 | (value & 0x0FFF));  
    digitalWrite(DAC_CS, HIGH);
}
