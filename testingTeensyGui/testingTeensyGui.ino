#include <SPI.h>
#include <Wire.h>
#include <IntervalTimer.h>
#include <Adafruit_INA219.h>
#include "HX711.h"

// ========================================================
// PIN DEFINITIONS
// ========================================================
#define DAC_CS   10
#define SYNC_PIN 5

// HX711 PINS
#define HX_DT   6
#define HX_SCK  7

// INA219 I2C address
const uint8_t INA_ADDR = 0x40;
const uint8_t INA_BYTES = 2;

// ========================================================
// CHIRP BUFFER
// ========================================================
#define BUFFER_SIZE 60000

DMAMEM uint16_t dacBuffer[BUFFER_SIZE];
volatile uint32_t sampleIndex = 0;
volatile uint32_t totalSamples = 0;
volatile bool     playing      = false;

IntervalTimer dacTimer;
IntervalTimer inaTimer;

// ========================================================
// INA219 BUFFERS
// ========================================================
Adafruit_INA219 ina219;

#define INA_BUFFER_SIZE 20000
DMAMEM double currentBuffer[INA_BUFFER_SIZE];
DMAMEM double timeBuffer[INA_BUFFER_SIZE];
volatile uint32_t inaCount = 0;

volatile bool collectingINA  = false;
volatile bool inaSampleFlag  = false;
uint32_t chirpStart_us       = 0;

const uint32_t INA_INTERVAL_US = 400;  // 2.5 kHz target

// ========================================================
// HX711 BUFFERS (80 Hz)
// ========================================================
HX711 hx;

#define HX_BUFFER_SIZE 400
DMAMEM double hxTime[HX_BUFFER_SIZE];
DMAMEM double hxVal[HX_BUFFER_SIZE];
volatile uint32_t hxCount = 0;

long hx_offset = 0;

// ========================================================
// SPI DAC WRITE
// ========================================================
void writeDAC(uint16_t value) {
  digitalWrite(DAC_CS, LOW);
  SPI.transfer16(0x3000 | (value & 0x0FFF));
  digitalWrite(DAC_CS, HIGH);
}

// ========================================================
// TIMER ISR — OUTPUT ONE DAC SAMPLE
// ========================================================
void outputSampleISR() {
  if (sampleIndex < totalSamples) {
    writeDAC(dacBuffer[sampleIndex++]);
  } else {
    // END OF CHIRP
    dacTimer.end();
    writeDAC(2048);
    digitalWrite(SYNC_PIN, LOW);

    collectingINA = false;
    inaTimer.end();
    playing = false;

    // ----------------------------
    // OUTPUT INA DATA
    // ----------------------------
    Serial.print("INA ");
    Serial.println(inaCount);
    for (uint32_t i = 0; i < inaCount; i++) {
      Serial.print(timeBuffer[i], 6);
      Serial.print(",");
      Serial.println(currentBuffer[i], 6);
    }

    // ----------------------------
    // OUTPUT HX DATA
    // ----------------------------
    Serial.print("HX ");
    Serial.println(hxCount);
    for (uint32_t i = 0; i < hxCount; i++) {
      Serial.print(hxTime[i], 6);
      Serial.print(",");
      Serial.println(hxVal[i], 4);
    }

    Serial.println("DONE");
  }
}

// ========================================================
// INA TIMER ISR
// ========================================================
void inaTimerISR() {
  if (collectingINA) {
    inaSampleFlag = true;
  }
}

// ========================================================
// CHIRP GENERATION
// ========================================================
float f1  = 30.0f;
float f2  = 350.0f;
float dur = 5.0f;
float env = 0.8f;
float fs  = 10000.0f;

void generateChirp() {
  totalSamples = (uint32_t)(dur * fs);
  if (totalSamples > BUFFER_SIZE)
    totalSamples = BUFFER_SIZE;

  float k = (f2 - f1) / dur;

  for (uint32_t n = 0; n < totalSamples; n++) {
    float t = (float)n / fs;
    float phase = 2.0f * PI * (f1 * t + 0.5f * k * t * t);
    float s = env * sinf(phase);

    float x = 0.5f * (s + 1.0f);
    if (x < 0) x = 0;
    if (x > 1) x = 1;

    dacBuffer[n] = (uint16_t)(x * 4095);
  }
}

// ========================================================
// INA219 HIGH-SPEED CONFIG
// ========================================================
void configureINA219HighSpeed() {
  ina219.begin();
  ina219.setCalibration_32V_2A();

  const uint8_t REG_CONFIG = 0x00;
  uint16_t config =
      (1 << 13) |
      (0 << 11) |
      (0x0 << 7) |
      (0x0 << 3) |
      0x7;

  Wire.beginTransmission(INA_ADDR);
  Wire.write(REG_CONFIG);
  Wire.write(config >> 8);
  Wire.write(config & 0xFF);
  Wire.endTransmission();
}

// ========================================================
// INA219 RAW
// ========================================================
bool readINA219Once(double &current_A) {
  const uint8_t REG_SHUNT = 0x01;

  Wire.beginTransmission(INA_ADDR);
  Wire.write(REG_SHUNT);
  Wire.endTransmission(false);

  Wire.requestFrom(INA_ADDR, INA_BYTES);
  if (Wire.available() < 2) return false;

  int16_t raw = (Wire.read() << 8) | Wire.read();
  float shunt_mV = raw * 0.01f;
  float I = (shunt_mV / 1000.0f) / 0.33f;
  current_A = I;
  return true;
}

// ========================================================
// HANDLE CHIRP COMMAND
// ========================================================
void handleChirpCommand(const String &line) {
  char cmd[16];
  float pf1, pf2, pdur, penv, pfs;

  int nParsed = sscanf(line.c_str(), "%15s %f %f %f %f %f",
                       cmd, &pf1, &pf2, &pdur, &penv, &pfs);

  if (nParsed != 6 || String(cmd) != "CHIRP") {
    Serial.println("ERR BAD_CMD");
    return;
  }

  f1 = pf1; f2 = pf2; dur = pdur; env = penv; fs = pfs;

  if (fs <= 0) fs = 10000;
  if (dur <= 0) dur = 5.0;
  if (env < 0) env = 0;
  if (env > 1) env = 1;

  generateChirp();

  sampleIndex = 0;
  inaCount = 0;
  hxCount  = 0;
  playing = true;
  collectingINA = true;
  inaSampleFlag = false;

  chirpStart_us = micros();

  float dt_us = 1e6f / fs;
  digitalWrite(SYNC_PIN, HIGH);
  dacTimer.begin(outputSampleISR, (unsigned long)dt_us);
  inaTimer.begin(inaTimerISR, INA_INTERVAL_US);

  Serial.println("STARTED");
}

// ========================================================
// SETUP
// ========================================================
void setup() {
  Serial.begin(115200);
  delay(200);

  Wire.begin();
  SPI.begin();

  pinMode(DAC_CS, OUTPUT);
  digitalWrite(DAC_CS, HIGH);

  pinMode(SYNC_PIN, OUTPUT);
  digitalWrite(SYNC_PIN, LOW);

  writeDAC(2048);

  configureINA219HighSpeed();

  // -------------------------
  // HX711 SETUP
  // -------------------------
  hx.begin(HX_DT, HX_SCK);

  // Force 80 Hz mode
  pinMode(HX_SCK, OUTPUT);
  digitalWrite(HX_SCK, HIGH);
  delay(1);
  digitalWrite(HX_SCK, LOW);
  delay(5);

  hx.set_gain(128);
  delay(50);

  long sum = 0;
  for (int i = 0; i < 50; i++) {
    while (!hx.is_ready());
    sum += hx.read();
  }
  hx_offset = sum / 50;

  Serial.print("HX711 offset=");
  Serial.println(hx_offset);

  Serial.println("READY");
}

// ========================================================
// MAIN LOOP
// ========================================================
void loop() {

  // -------------------------------
  // Command parser
  // -------------------------------
  if (Serial.available()) {
    String line = Serial.readStringUntil('\n');
    line.trim();
    if (!playing && line.startsWith("CHIRP"))
      handleChirpCommand(line);
  }

  // -------------------------------
  // INA219 samples
  // -------------------------------
  if (collectingINA && inaSampleFlag) {
    inaSampleFlag = false;

    if (inaCount < INA_BUFFER_SIZE) {
      double current_A;
      if (readINA219Once(current_A)) {
        double t = (double)(micros() - chirpStart_us) * 1e-6;
        timeBuffer[inaCount]    = t;
        currentBuffer[inaCount] = current_A;
        inaCount++;
      }
    }
  }

  // -------------------------------
  // HX711 samples (~80 Hz)
  // -------------------------------
  if (collectingINA && hx.is_ready()) {
    if (hxCount < HX_BUFFER_SIZE) {
      long raw = hx.read() - hx_offset;
      double t = (double)(micros() - chirpStart_us) * 1e-6;
      hxTime[hxCount] = t;
      hxVal[hxCount]  = (double)raw;
      hxCount++;
    }
  }
}
