#include "MCP492X.h"
#include <math.h>

#define PIN_SPI_CHIP_SELECT_DAC 10
#define DAC_MAX 4095
#define DAC_REF_VOLTAGE 3.3

MCP492X myDac(PIN_SPI_CHIP_SELECT_DAC);

void setup() {
  myDac.begin();
  Serial.begin(9600);
}

// Generate a sine wave
void sinGen(float frequency, float amplitude) {
  const int maxSamples = 1000;
  float period = 1.0e6 / frequency;
  float sampleInterval = period / maxSamples;

  uint16_t sineTable[maxSamples];
  
  for (int i = 0; i < maxSamples; i++) {
    float rawSine = sin(2 * PI * i / maxSamples);
    float scaledSine = (rawSine + 1.0) / 2.0;
    float val = (uint16_t)(scaledSine * (amplitude / DAC_REF_VOLTAGE) * DAC_MAX);
    sineTable[i] = val;
  }

  unsigned long lastTime = micros();
  int index = 0;
  unsigned long startTime = micros();
  while (true) {
    unsigned long now = micros();
    if (now - lastTime >= sampleInterval) {
      lastTime += sampleInterval;
      myDac.analogWrite(sineTable[index]);
      index = (index + 1) % maxSamples;
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
  sinGen(250, 2);
}
