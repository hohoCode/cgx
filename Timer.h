#ifndef MYTIMER_H
#define MYTIMER_H

typedef struct mytimer_s {

   double start_ms;
   double stop_ms;

} mytimer_t;


void
timer_start(mytimer_t * timer);

double
timer_stop(mytimer_t * timer);

double
timer_elapsed(mytimer_t * timer);

#endif /* MYTIMER_H */
