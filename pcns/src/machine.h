// -*- C++ -*-
#pragma once

struct machine_t {
  static void initialize(int argc, char** argv);
  static void finalize();
  static int rank();
  static int size();
  static void all_gather(const int * sendbuf, int sendcount,
                         int * recvbuf, int recvcount);
};
