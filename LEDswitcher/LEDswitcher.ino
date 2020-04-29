
// Teensy3.2 code to trigger blue and violet LEDs for widefield imaging.
// Requires three digital inputs for operation:
// Exposure signal from the imaging camera (IN_STROBE).
// Three separate digital lines to switch between all blue, all violet or alternating illumination. Only one input should be ON at any time.
// Provides two outputs to trigger blue and violet LEDs
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Define digital lines
#define IN_STROBE 3 // strobe signal from camera
#define OUT_BLUE 5 // digial line to control blue LED
#define OUT_VIOLET 6 // digital line to control violet LED
#define IN_STROBESWITCH 10 // trigger signal for strobe mode - uses strobe signal to switch between blue and violet LED
#define IN_BLUESWITCH 11 // trigger signal to switch on blue LED
#define IN_VIOLETSWITCH 12 // trigger signal to switch on violet LED

// other variables
int blueCounter = 0; // counter for 'blue frames'
int blueFrames = 1; // number of 'blue frames' before switching to violet. When using more than 1 blue frame, loss of LED output power due to heating should be controlled to avoid artifacts.
int violetCounter = 0; // counter for 'violet frames'
int violetFrames = 1; // number of 'violet frames' before switching to blue
int deadDuration = 1; // duration of dead-time in ms.
bool noSwitch = false; // flag for dead-time
bool switchLight = true; // flag when expecting strobe signal. Set to false after switching LEDs and gets refreshed when strobe reads low.
bool currentBlue = true; // flag to indicate LED state. true if blue is on, false if violet is on.
unsigned long clocker = millis(); // timer to create dead-time after receiving a strobe signal.
bool strobeOn = false; // flag for strobe signal. useful to avoid multiple reads from trigger line.
bool blueSwitch = false; // flag for blue switch.
bool violSwitch = false; // flag for violet switch.
bool strobeSwitch = false; // flag for strobe switch.

//////////////////////////////////////////////////////////////////////////////////////////////////

void setup() {
  // set up digital lines
  pinMode(OUT_BLUE, OUTPUT); // signal for blue LED
  pinMode(OUT_VIOLET, OUTPUT); // signal for violet LED
  
  pinMode(IN_STROBE, INPUT); 
  pinMode(IN_STROBESWITCH, INPUT);
  pinMode(IN_BLUESWITCH, INPUT);
  pinMode(IN_VIOLETSWITCH, INPUT);
}

//////////////////////////////////////////////////////////////////////////////////////////////////

void loop() {

  strobeOn = digitalReadFast(IN_STROBE); // only produce output when strobe signal is high

  if (!strobeOn) { // LEDs should only be on if strobe signal is high
    digitalWriteFast(OUT_BLUE,LOW);
    digitalWriteFast(OUT_VIOLET,LOW);
    blueSwitch = digitalReadFast(IN_BLUESWITCH);
    violSwitch = digitalReadFast(IN_VIOLETSWITCH);
    strobeSwitch = digitalReadFast(IN_STROBESWITCH);
  }

  // check blue, violet and strobe switch
  if (blueSwitch) { // switch on blue light
    if (strobeOn) {
      digitalWriteFast(OUT_BLUE,HIGH);
      digitalWriteFast(OUT_VIOLET,LOW);
    }
    
    // reset counters and ensure blue light is first when switching to strobe mode
    violetCounter = 0;
    blueCounter = 0;
    currentBlue = true;
  }

  else if (violSwitch) { // switch on violet light
    if (strobeOn) {
      digitalWriteFast(OUT_BLUE,LOW);
      digitalWriteFast(OUT_VIOLET,HIGH);
    }
      
    // reset counters and ensure blue light is first when switching to strobe mode
    violetCounter = 0;
    blueCounter = 0;
    currentBlue = true;
  }
  
  else if (strobeSwitch) { // run strobe mode
    // control state of output signals
      if (strobeOn) { // received strobe signal
          if (currentBlue) {
            digitalWriteFast(OUT_BLUE,HIGH);
            digitalWriteFast(OUT_VIOLET,LOW);
          }
          else {
            digitalWriteFast(OUT_BLUE,LOW);
            digitalWriteFast(OUT_VIOLET,HIGH);
          }
      }
      
    // determine switch to different light color
    if (noSwitch) { // wait for duation of dead-time before reading strobe again.
        if ((millis()-clocker) >= deadDuration) {
          noSwitch = false;
        }
    }
    
    else { // dead-time has passed. check strobe signal.
      if (strobeOn) { // received strobe signal
               
        if (switchLight) { // switch LEDs or increase counter
          clocker = millis(); // reset clocker
          switchLight = false; // dont switch again until strobe signal has gone low
          noSwitch = true; // dont switch again until dead-time has passed
          
          if (currentBlue) {
            ++ blueCounter;            
            if (blueCounter == blueFrames) { // check if requested amount of blue illuminated frames has been reached
              blueCounter = 0;
              currentBlue = false;
            }
          }
          
          else {
            ++ violetCounter;            
            if (violetCounter == violetFrames) { // check if requested amount of violet illuminated frames has been reached
              violetCounter = 0;
              currentBlue = true;
            }
          }
        }
      }
      
      else { // strobe line is low
        switchLight = true;
      }
    }
  }

  else { // no trigger signals are present. Both LEDs are off and reset counters.
    digitalWriteFast(OUT_BLUE, LOW);
    digitalWriteFast(OUT_VIOLET, LOW);
    violetCounter = 0;
    blueCounter = 0;
    currentBlue = true;
  }
}
