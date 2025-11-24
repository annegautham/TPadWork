#include "MCP492X.h"
#include <math.h>
#include <Wire.h>
#include <Adafruit_INA219.h>

#define PIN_SPI_CHIP_SELECT_DAC 10
#define DAC_MAX 4095
#define DAC_REF_VOLTAGE 3.3

MCP492X myDac(PIN_SPI_CHIP_SELECT_DAC);
Adafruit_INA219 ina219;

// Your shunt resistor (in ohms)
const float R_SHUNT = 0.33;

// How often to print current (microseconds)
const unsigned long CURRENT_PRINT_INTERVAL = 1000;   // 2 ms = 500 Hz print rate

void setup() {
  Serial.begin(115200);
  Wire.begin();

  // Initialize DAC
  myDac.begin();

  // Initialize INA219
  if (!ina219.begin()) {
    Serial.println("INA219 not found!");
    while (1);
  }

  // Optional: increase resolution
  ina219.setCalibration_16V_400mA();

  Serial.println("System ready!");
}


// Generate a sine wave
// void sinGen(float frequency, float amplitude) {
//   const int maxSamples = 1000;
//   float period = 1.0e6 / frequency;
//   float sampleInterval = period / maxSamples;

//   uint16_t sineTable[maxSamples];
  
//   for (int i = 0; i < maxSamples; i++) {
//     float rawSine = sin(2 * PI * i / maxSamples);
//     float scaledSine = (rawSine + 1.0) / 2.0;
//     float val = (uint16_t)(scaledSine * (amplitude / DAC_REF_VOLTAGE) * DAC_MAX);
//     sineTable[i] = val;
//   }

//   unsigned long lastTime = micros();
//   int index = 0;
//   unsigned long startTime = micros();
//   while (true) {
//     unsigned long now = micros();
//     if (now - lastTime >= sampleInterval) {
//       lastTime += sampleInterval;
//       myDac.analogWrite(sineTable[index]);
//       index = (index + 1) % maxSamples;
//     }
//   }
// }

void sinGen(float frequency, float amplitude) {
  const int maxSamples = 1000;
  float period = 1.0e6 / frequency;
  float sampleInterval = period / maxSamples;

  uint16_t sineTable[maxSamples];

  // Precompute sine table
  for (int i = 0; i < maxSamples; i++) {
    float rawSine = sin(2 * PI * i / maxSamples);
    float scaledSine = (rawSine + 1.0) / 2.0;
    float val = scaledSine * (amplitude / DAC_REF_VOLTAGE) * DAC_MAX;
    sineTable[i] = (uint16_t)val;
  }

  unsigned long lastTime = micros();
  unsigned long lastCurrentTime = micros();
  int index = 0;

  while (true) {
    unsigned long now = micros();

    // ----- DAC WAVEFORM OUTPUT -----
    if (now - lastTime >= sampleInterval) {
      lastTime += sampleInterval;
      myDac.analogWrite(sineTable[index]);
      index = (index + 1) % maxSamples;
    }

    // ----- CURRENT MEASUREMENT -----
    if (now - lastCurrentTime >= CURRENT_PRINT_INTERVAL) {
      lastCurrentTime = now;

      float current_mA = ina219.getCurrent_mA();
      float current_A = current_mA / 1000.0;

      // Convert to force later using F = B*L*I

      Serial.print("I = ");
      Serial.print(current_A, 4);
      Serial.println(" A");
    }
  }
}

// Generate a chirp signal
void chirpGen(float freq1, float finalFreq, float timeDelay, float amplitude) {
  const int maxSamples = 100000;
  unsigned long durationMicros = timeDelay * 1e6;
  float sampleInterval = durationMicros / maxSamples;

  uint16_t sineTable[maxSamples];
  for (int i = 0; i < maxSamples; i++) {
    float t = (float)i / maxSamples;
    float currentFreq = freq1 + t * (finalFreq - freq1);
    float phase = sin(2 * PI * currentFreq * t * timeDelay);
    float scaledSine = (phase + 1.0) / 2.0;
    sineTable[i] = (uint16_t)(scaledSine * (amplitude / DAC_REF_VOLTAGE) * DAC_MAX);
  }
  
  int index = 0;
  unsigned long lastTime = micros();
  unsigned long startTime = lastTime;
  while (micros() - startTime < durationMicros) {
    unsigned long now = micros();
    if (now - lastTime >= sampleInterval) {
      lastTime += sampleInterval;
      myDac.analogWrite(sineTable[index]);
      index = (index + 1) % maxSamples;
    }
  }
}

void loop() {
  sinGen(20, 0.8);
}
