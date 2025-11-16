#include <SPI.h>
#include <IntervalTimer.h>

#define DAC_CS 10
#define BUFFER_SIZE 100000

int buffer[BUFFER_SIZE];
volatile int sampleIndex = 0;
volatile int totalSamples = 0;

IntervalTimer dacTimer; 

void setup() {
    Serial.begin(115200);
    SPI.begin();
    pinMode(DAC_CS, OUTPUT);
    digitalWrite(DAC_CS, HIGH);
}

void loop() {
    int count = 0;

    while (Serial.available()) {
        String data = Serial.readStringUntil('\n');
        delay(2000);
        if (data == "END") {
            Serial.flush();  
            delay(500);  // Ensure all serial data is processed

            totalSamples = count;
            sampleIndex = 0;

            // Start timer to output samples at 10kHz (every 100us)
            dacTimer.begin(outputSampleISR, 100);  // in microseconds
            while (sampleIndex < totalSamples);    // Wait until playback is done
            dacTimer.end();  // Stop timer when done

            return;
        }

        char *ptr = strtok((char*)data.c_str(), ",");
        while (ptr != NULL && count < BUFFER_SIZE) {
            buffer[count++] = atoi(ptr);
            ptr = strtok(NULL, ",");
        }
    }
}

void outputSampleISR() {
    if (sampleIndex < totalSamples) {
        writeDAC(buffer[sampleIndex++]);
    }
}

void writeDAC(uint16_t value) {
    digitalWrite(DAC_CS, LOW);
    SPI.transfer16(0x3000 | (value & 0x0FFF));
    digitalWrite(DAC_CS, HIGH);
}
