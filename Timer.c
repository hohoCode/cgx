/* system includes */
#include <sys/time.h>
#include "Timer.h"
#include "ComTypes.h"


double 
get_time_ms() {

   static struct timeval tval;
   gettimeofday(&tval, NULL);
   return (tval.tv_sec * 1000 + tval.tv_usec/1000.0);

}


void
timer_start(mytimer_t * timer) {

   timer->start_ms = get_time_ms();
   timer->stop_ms = 0;
}


double
timer_stop(mytimer_t * timer) {

   timer->stop_ms = get_time_ms();
   return (timer->stop_ms - timer->start_ms);
}


double
timer_elapsed(mytimer_t * timer) {

   return (timer->stop_ms - timer->start_ms);

}
