#ifndef grad_h
#define grad_h

struct grad_t {
  enum gd_mthd;
  enum lm_mthd;

  double*** gd;
}
