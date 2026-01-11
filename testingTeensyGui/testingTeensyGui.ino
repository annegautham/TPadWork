#include <SPI.h>
#include <Wire.h>
#include <IntervalTimer.h>
#include <Adafruit_INA219.h>

// ========================================================
// PIN DEFINITIONS
// ========================================================
#define DAC_CS   10
#define SYNC_PIN 5

// HX711 pins
#define HX_DT   6
#define HX_SCK  7

// INA219
const uint8_t INA_ADDR  = 0x40;

// ========================================================
// CHIRP BUFFER
// ========================================================
#define BUFFER_SIZE 60000

DMAMEM uint16_t dacBuffer[BUFFER_SIZE];
volatile uint32_t sampleIndex   = 0;
volatile uint32_t totalSamples  = 0;
volatile bool     playing       = false;

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

const uint32_t INA_INTERVAL_US = 400; // ~2.5 kHz


// ========================================================
// WRITE TO MCP4921 DAC
// ========================================================
void writeDAC(uint16_t value) {
  digitalWrite(DAC_CS, LOW);
  SPI.transfer16(0x3000 | (value & 0x0FFF));
  digitalWrite(DAC_CS, HIGH);
}

// ========================================================
// DAC TIMER ISR — OUTPUT ONE SAMPLE
// ========================================================
void outputSampleISR() {
  if (sampleIndex < totalSamples) {

    writeDAC(dacBuffer[sampleIndex++]);

  } else {
    // End chirp
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

    Serial.println("DONE");
  }
}

// ========================================================
// INA TIMER ISR — Set flag for main loop
// ========================================================
void inaTimerISR() {
  if (collectingINA)
    inaSampleFlag = true;
}

// ========================================================
// GENERATE CHIRP
// ========================================================
float f1  = 30;
float f2  = 350;
float dur = 5;
float env = 0.8;
float fs  = 10000;

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
// INA219 HIGH-SPEED CONFIG FIXED
// ========================================================
void configureINA219HighSpeed() {

  ina219.begin();
  ina219.setCalibration_32V_2A();   // defines Current_LSB = 100 µA

  const uint8_t REG_CONFIG = 0x00;

  // Correct settings:
  // BRNG = 32V   → 1
  // PG   = ±320mV → 3
  // BADC = 12-bit → 3
  // SADC = 12-bit → 3
  // MODE = shunt continuous → 7

  uint16_t config =
      (1 << 13) |   // BRNG = 32V
      (3 << 11) |   // PG = ±320mV
      (3 << 7)  |   // BADC = 12-bit
      (3 << 3)  |   // SADC = 12-bit
      0x7;          // MODE = continuous shunt

  Wire.beginTransmission(INA_ADDR);
  Wire.write(REG_CONFIG);
  Wire.write(config >> 8);
  Wire.write(config & 0xFF);
  Wire.endTransmission();
}

// ========================================================
// READ CURRENT REGISTER (0x04) — FIXED
// ========================================================
bool readINA219Once(double &current_A) {

  const uint8_t REG_CURRENT = 0x04;  // <-- IMPORTANT

  Wire.beginTransmission(INA_ADDR);
  Wire.write(REG_CURRENT);
  Wire.endTransmission(false);

  Wire.requestFrom(INA_ADDR, (uint8_t)2);
  if (Wire.available() < 2)
    return false;

  int16_t raw = (Wire.read() << 8) | Wire.read();

  // From calibration_32V_2A: Current_LSB = 100 µA
  current_A = raw * 0.0001;

  return true;
}

// ========================================================
// HANDLE 'CHIRP ...' COMMAND
// ========================================================
void handleChirpCommand(const String &line) {
  char cmd[16];
  float pf1, pf2, pdur, penv, pfs;

  int n = sscanf(line.c_str(), "%15s %f %f %f %f %f",
                 cmd, &pf1, &pf2, &pdur, &penv, &pfs);

  if (n != 6 || String(cmd) != "CHIRP") {
    Serial.println("ERR BAD_CMD");
    return;
  }

  f1 = pf1;
  f2 = pf2;
  dur = pdur;
  env = penv;
  fs = pfs;

  generateChirp();

  sampleIndex   = 0;
  inaCount      = 0;
  playing       = true;
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


  Serial.println("READY");
}

// ========================================================
// MAIN LOOP
// ========================================================
void loop() {

  // Receive commands
  if (Serial.available()) {
    String line = Serial.readStringUntil('\n');
    line.trim();
    if (!playing && line.startsWith("CHIRP"))
      handleChirpCommand(line);
  }

  // INA219 sampling
  if (collectingINA && inaSampleFlag) {
    inaSampleFlag = false;

    if (inaCount < INA_BUFFER_SIZE) {
      double I;
      if (readINA219Once(I)) {
        double t = (double)(micros() - chirpStart_us) * 1e-6;
        timeBuffer[inaCount]    = t;
        currentBuffer[inaCount] = I;
        inaCount++;
      }
    }
  }

}
