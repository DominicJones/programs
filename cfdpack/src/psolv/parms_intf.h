#ifndef parms_intf_h
#define parms_intf_h

#ifdef __cplusplus
extern "C" {
#endif

void parms_init_ (int* mp);
void parms_solv_ (int* mp, int* im, int* ia, int* ja,
                  double* aij, double* rhs, double* sol,
                  double* res0);

#ifdef __cplusplus
}
#endif

#endif
