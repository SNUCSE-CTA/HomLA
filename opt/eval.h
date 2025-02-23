#ifndef EVAL_H
#define EVAL_H

#include "core.h"

void score(LS ***H, LS ***Hx, LS ***Hy, LS **xa, LS **xb, int n, int m, int s_m,
           int s_s, int g_o, int g_e, int w_sco, int w_pos, int w_alp, CK *ck);

void find(LS *, LS *, LS *, LS *, LS *, LS ***, LS ***, LS ***, int, int, int,
          int, CK *);

#endif
