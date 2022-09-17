/* G77 fortran callable usleep */
#include <unistd.h>
void usleep_(usecs)
     long *usecs;
{
  usleep(  (unsigned long) *usecs);
}
