#include <SPI.h>

const int CS_PIN = 10; // Chip Select pin for MCP4921

void setup() {
  Serial.begin(115200);
  SPI.begin();
  pinMode(CS_PIN, OUTPUT);
  digitalWrite(CS_PIN, HIGH);
}

void loop() {
  if (Serial.available() >= 2) {
    // Read 16-bit value from serial
    uint16_t value = Serial.read() << 8 | Serial.read();
    
    // Send to MCP4921
    digitalWrite(CS_PIN, LOW);
    SPI.transfer16(0x3000 | (value & 0x0FFF)); // 0x3000 for MCP4921 config bits
    digitalWrite(CS_PIN, HIGH);
  }
}
  