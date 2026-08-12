#include "EventClock.h"
using namespace std;

int main() {
    EventClock<TimeTicks::Seconds> eventClock;
    eventClock.startClock("test");
    for(unsigned int i = 0; i < 4000000000; ++i) ;
    eventClock.stopClock("test");
    eventClock.printEvent("test");

    eventClock.startClock("test");
    for(unsigned int i = 0; i < 4000000000; ++i) ;
    eventClock.stopClock("test");
    eventClock.printEvent("test");

    return 0;
}
