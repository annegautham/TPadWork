#include "HX711.h"

#define DT_PIN   6   // DT
#define SCK_PIN  7   // SCK

HX711 scale;

void setup() {
  Serial.begin(115200);
  delay(500);

  Serial.println("Starting HX711...");

  pinMode(SCK_PIN, OUTPUT);
  digitalWrite(SCK_PIN, LOW);

  scale.begin(DT_PIN, SCK_PIN);
  scale.set_gain(128);

  Serial.println("Initialized. Reading...");
}

void loop() {
  bool rdy = scale.is_ready();
  // Serial.print("Ready? ");
  // Serial.println(rdy);

  if (rdy) {
    long val = scale.read();
    // Serial.print("Raw: ");
    Serial.println(val);
  }

  delay(100);
}
