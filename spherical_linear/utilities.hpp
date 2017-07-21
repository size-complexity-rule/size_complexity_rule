#ifndef __UTILITIES_HPP__
#define __UTILITIES_HPP__

#include <ctime>
#include <string>

#include <unistd.h>

int get_process_id() {
  return int( getpid() );
}

std::string get_hostname() {
  char hostname[256];
  int  gethostname_sucess = (gethostname(hostname, 256) == 0)?true:false;

  return std::string(gethostname_sucess?hostname:"");
}

std::string get_time_and_date() {
  time_t rawtime = std::time(nullptr);
  auto ct = std::localtime( &rawtime );

  std::string time_str =
      std::to_string(ct->tm_hour)       + std::string(":")
    + std::to_string(ct->tm_min)        + std::string(":")
    + std::to_string(ct->tm_sec)        + std::string("\t")
    + std::to_string(ct->tm_mday)       + std::string("/")
    + std::to_string(ct->tm_mon+1)      + std::string("/")
    + std::to_string(ct->tm_year+1900);
  return time_str;
}

#endif /* end of include guard: __UTILITIES_HPP__ */
