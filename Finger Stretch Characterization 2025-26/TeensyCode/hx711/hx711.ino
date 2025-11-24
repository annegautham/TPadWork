#include "HX711.h"

#define DT_PIN 6
#define SCK_PIN 7

HX711 scale;

void setup() {
  Serial.begin(115200);
  scale.begin(DT_PIN, SCK_PIN);
  Serial.println("HX711 ready");
}

void loop() {

  // -------- Read Channel A (4-wire load cell) --------
  scale.set_gain(128);     // Channel A @ 128x
  delayMicroseconds(10);   // allow mux to settle
  long readingA = scale.read();

  // -------- Read Channel B (3-wire load cell) --------
  scale.set_gain(32);      // Channel B @ 32x
  delayMicroseconds(10);   // allow mux to settle
  long readingB = scale.read();

  // -------- Print --------
  Serial.print("A = ");
  Serial.print(readingA);
  Serial.print("    B = ");
  Serial.println(readingB);

  delay(10);
}

