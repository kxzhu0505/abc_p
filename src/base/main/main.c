#include <misc/util/abc_global.h>
#include <gperftools/profiler.h>
#include <sys/time.h>
#include <unistd.h>
#define PROFILER_ON

ABC_NAMESPACE_IMPL_START

int Abc_RealMain(int argc, char *argv[]);

ABC_NAMESPACE_IMPL_END

int main(int argc, char *argv[])
{
#ifdef PROFILER_ON
   ProfilerStart("/home/ymc/prof/abc.prof");
#endif
   struct timeval t1,t2;
   double timeUsed;
   gettimeofday(&t1, NULL);

   int ret = ABC_NAMESPACE_PREFIX Abc_RealMain(argc, argv);

   gettimeofday(&t2, NULL);
   timeUsed = t2.tv_sec - t1.tv_sec + (t2.tv_usec - t1.tv_usec)/1000000.0;
   printf("total time used: %f\n", timeUsed);

#ifdef PROFILER_ON
   ProfilerStop();  
#endif
   return ret;
}
