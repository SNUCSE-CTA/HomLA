#include "eval.h"
#include "par_circuit.h"
#include "seq_circuit.h"
#include <iostream>

using namespace par;

// Pre-computing similarity scores
void sim_score(LS ***S, LS **x, LS **y, int n, int m, int w_sco, int w_alp,
               int s_m, int s_s, CK *ck) {

    LS *es_m = new_gate_bootstrapping_ciphertext_array(w_sco, ck->params);
    LS *es_s = new_gate_bootstrapping_ciphertext_array(w_sco, ck->params);
    seq::enc_const(es_m, s_m, w_sco, ck);
    seq::enc_const(es_s, s_s, w_sco, ck);

#pragma omp parallel for collapse(2)
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            LS *eq = new_gate_bootstrapping_ciphertext(ck->params);
            equals(eq, x[i], y[j], w_alp, ck);
            select(S[i][j], eq, es_m, es_s, w_sco, ck);
            delete_gate_bootstrapping_ciphertext(eq);
        }
    }

    delete_gate_bootstrapping_ciphertext_array(w_sco, es_m);
    delete_gate_bootstrapping_ciphertext_array(w_sco, es_s);
}

// Computing scores and positions by diagonal-parallel way
void score(LS ***H, LS ***Hx, LS ***Hy, LS **cx, LS **cy, int n, int m, int s_m,
           int s_s, int g_o, int g_e, int w_sco, int w_pos, int w_alp, CK *ck) {

    // Encrypt params
    LS *eg_o = new_gate_bootstrapping_ciphertext_array(w_sco, ck->params);
    LS *eg_e = new_gate_bootstrapping_ciphertext_array(w_sco, ck->params);
    LS *zero = new_gate_bootstrapping_ciphertext_array(w_sco, ck->params);

    seq::enc_const(eg_o, g_o, w_sco, ck);
    seq::enc_const(eg_e, g_e, w_sco, ck);
    seq::fill_zero(zero, w_sco, ck);

    // INIT TABLES

    LS ***P = new LS **[n + 1];
    LS ***Q = new LS **[n + 1];

#ifdef POSITION
    LS ***Px = new LS **[n + 1];
    LS ***Py = new LS **[n + 1];
    LS ***Qx = new LS **[n + 1];
    LS ***Qy = new LS **[n + 1];
#endif

    for (int i = 0; i <= n; i++) {

        P[i] = new LS *[m + 1];
        Q[i] = new LS *[m + 1];
#ifdef POSITION
        Px[i] = new LS *[m + 1];
        Py[i] = new LS *[m + 1];
        Qx[i] = new LS *[m + 1];
        Qy[i] = new LS *[m + 1];
#endif
        for (int j = 0; j <= m; j++) {

            P[i][j] =
                new_gate_bootstrapping_ciphertext_array(w_sco, ck->params);
            Q[i][j] =
                new_gate_bootstrapping_ciphertext_array(w_sco, ck->params);
#ifdef POSITION
            Px[i][j] =
                new_gate_bootstrapping_ciphertext_array(w_pos, ck->params);
            Py[i][j] =
                new_gate_bootstrapping_ciphertext_array(w_pos, ck->params);
            Qx[i][j] =
                new_gate_bootstrapping_ciphertext_array(w_pos, ck->params);
            Qy[i][j] =
                new_gate_bootstrapping_ciphertext_array(w_pos, ck->params);
#endif
        }
    }

    for (int i = 0; i <= n; i++) {
        seq::copy(H[i][0], zero, w_sco, ck);
        seq::copy(P[i][0], zero, w_sco, ck);
        seq::copy(Q[i][0], zero, w_sco, ck);

        seq::enc_const(Hx[i][0], i, w_pos, ck);
        seq::enc_const(Hy[i][0], 0, w_pos, ck);
#ifdef POSITION
        seq::enc_const(Qx[i][0], i, w_pos, ck);
        seq::enc_const(Qy[i][0], 0, w_pos, ck);
#endif
    }

    for (int j = 0; j <= m; j++) {
        seq::copy(H[0][j], zero, w_sco, ck);
        seq::copy(P[0][j], zero, w_sco, ck);
        seq::copy(Q[0][j], zero, w_sco, ck);

        seq::enc_const(Hx[0][j], 0, w_pos, ck);
        seq::enc_const(Hy[0][j], j, w_pos, ck);
#ifdef POSITION
        seq::enc_const(Px[0][j], 0, w_pos, ck);
        seq::enc_const(Py[0][j], j, w_pos, ck);
#endif
    }

    int start, end;

    LS ***S = new LS **[n];
    for (int i = 0; i < n; i++) {
        S[i] = new LS *[m];
        for (int j = 0; j < m; j++) {
            S[i][j] =
                new_gate_bootstrapping_ciphertext_array(w_sco, ck->params);
        }
    }

    sim_score(S, cx, cy, n, m, w_sco, w_alp, s_m, s_s, ck);

    // SW reccurence computed by diagonal-parallel way
    for (int d = 1; d <= n + m - 1; d++) {
        end = std::min(d, n);
        start = std::max(1, d - m + 1);

#pragma omp parallel for
        for (int k = start; k <= end; k++) {

            int i = k;
            int j = d - k + 1;

            LS *cmp = new_gate_bootstrapping_ciphertext(ck->params);
            LS *tmp1 =
                new_gate_bootstrapping_ciphertext_array(w_sco, ck->params);
            LS *tmp2 =
                new_gate_bootstrapping_ciphertext_array(w_sco, ck->params);

            add(tmp1, P[i - 1][j], eg_e, w_sco, ck);
            add(tmp2, H[i - 1][j], eg_o, w_sco, ck);
            compare(cmp, tmp1, tmp2, w_sco, ck);
            select(P[i][j], cmp, tmp1, tmp2, w_sco, ck);
#ifdef POSITION
            select(Px[i][j], cmp, Px[i - 1][j], Hx[i - 1][j], w_pos, ck);
            select(Py[i][j], cmp, Py[i - 1][j], Hy[i - 1][j], w_pos, ck);
#endif
            add(tmp1, Q[i][j - 1], eg_e, w_sco, ck);
            add(tmp2, H[i][j - 1], eg_o, w_sco, ck);
            compare(cmp, tmp1, tmp2, w_sco, ck);
            select(Q[i][j], cmp, tmp1, tmp2, w_sco, ck);
#ifdef POSITION
            select(Qx[i][j], cmp, Qx[i][j - 1], Hx[i][j - 1], w_pos, ck);
            select(Qy[i][j], cmp, Qy[i][j - 1], Hy[i][j - 1], w_pos, ck);
#endif
            compare(cmp, P[i][j], Q[i][j], w_sco, ck);
            select(H[i][j], cmp, P[i][j], Q[i][j], w_sco, ck);
#ifdef POSITION
            select(Hx[i][j], cmp, Px[i][j], Qx[i][j], w_pos, ck);
            select(Hy[i][j], cmp, Py[i][j], Qy[i][j], w_pos, ck);
#endif
            add(tmp1, H[i - 1][j - 1], S[i - 1][j - 1], w_sco, ck);
            compare(cmp, tmp1, H[i][j], w_sco, ck);
            select(H[i][j], cmp, tmp1, H[i][j], w_sco, ck);
#ifdef POSITION
            select(Hx[i][j], cmp, Hx[i - 1][j - 1], Hx[i][j], w_pos, ck);
            select(Hy[i][j], cmp, Hy[i - 1][j - 1], Hy[i][j], w_pos, ck);
#endif
            compare(cmp, H[i][j], zero, w_sco, ck);
            select(H[i][j], cmp, H[i][j], zero, w_sco, ck);
#ifdef POSITION
            seq::enc_const(tmp1, i, w_pos, ck);
            seq::enc_const(tmp2, j, w_pos, ck);
            select(Hx[i][j], cmp, Hx[i][j], tmp1, w_pos, ck);
            select(Hy[i][j], cmp, Hy[i][j], tmp2, w_pos, ck);
#endif
            // FREE
            delete_gate_bootstrapping_ciphertext(cmp);
            delete_gate_bootstrapping_ciphertext_array(w_sco, tmp1);
            delete_gate_bootstrapping_ciphertext_array(w_sco, tmp2);
        }
    }

    delete_gate_bootstrapping_ciphertext_array(w_sco, eg_e);
    delete_gate_bootstrapping_ciphertext_array(w_sco, eg_o);
    delete_gate_bootstrapping_ciphertext_array(w_sco, zero);

    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= m; j++) {

            delete_gate_bootstrapping_ciphertext_array(w_sco, P[i][j]);
            delete_gate_bootstrapping_ciphertext_array(w_sco, Q[i][j]);
#ifdef POSITION
            delete_gate_bootstrapping_ciphertext_array(w_pos, Px[i][j]);
            delete_gate_bootstrapping_ciphertext_array(w_pos, Py[i][j]);
            delete_gate_bootstrapping_ciphertext_array(w_pos, Qx[i][j]);
            delete_gate_bootstrapping_ciphertext_array(w_pos, Qy[i][j]);
#endif
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            delete_gate_bootstrapping_ciphertext_array(w_sco, S[i][j]);
        }
    }
}

// Finding maximum score and its position via tree-based parallel computation
void find(LS *M, LS *start_x, LS *start_y, LS *end_x, LS *end_y, LS ***H,
          LS ***Hx, LS ***Hy, int n, int m, int w_sco, int w_pos, CK *ck) {

#ifdef POSITION
    LS ***end_i = new LS **[n + 1];
    LS ***end_j = new LS **[n + 1];

    for (int i = 0; i <= n; i++) {
        end_i[i] = new LS *[m + 1];
        end_j[i] = new LS *[m + 1];
        for (int j = 0; j <= m; j++) {
            end_i[i][j] =
                new_gate_bootstrapping_ciphertext_array(w_pos, ck->params);
            end_j[i][j] =
                new_gate_bootstrapping_ciphertext_array(w_pos, ck->params);

            seq::enc_const(end_i[i][j], i, w_pos, ck);
            seq::enc_const(end_j[i][j], j, w_pos, ck);
        }
    }
#endif

    for (int d = 1; d <= ceil(log2(n * m)); d++) {
#pragma omp parallel for
        for (int k = 1; k <= n * m; k = k + (1 << d)) {
            int i = (k - 1) / m + 1;
            int j = (k - 1) % m + 1;
            int ip = (k + (1 << (d - 1)) - 1) / m + 1;
            int jp = (k + (1 << (d - 1)) - 1) % m + 1;
            if (ip <= n) {
                LS *cmp = new_gate_bootstrapping_ciphertext(ck->params);

                compare(cmp, H[i][j], H[ip][jp], w_sco, ck);

                select(H[i][j], cmp, H[i][j], H[ip][jp], w_sco, ck);
#ifdef POSITION
                select(Hx[i][j], cmp, Hx[i][j], Hx[ip][jp], w_pos, ck);
                select(Hy[i][j], cmp, Hy[i][j], Hy[ip][jp], w_pos, ck);
                select(end_i[i][j], cmp, end_i[i][j], end_i[ip][jp], w_pos, ck);
                select(end_j[i][j], cmp, end_j[i][j], end_j[ip][jp], w_pos, ck);
#endif
                delete_gate_bootstrapping_ciphertext(cmp);
            }
        }
    }

    seq::copy(M, H[1][1], w_sco, ck);

#ifdef POSITION
    seq::copy(start_x, Hx[1][1], w_pos, ck);
    seq::copy(start_y, Hy[1][1], w_pos, ck);
    seq::copy(end_x, end_i[1][1], w_pos, ck);
    seq::copy(end_y, end_j[1][1], w_pos, ck);

    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= m; j++) {
            delete_gate_bootstrapping_ciphertext_array(w_pos, end_i[i][j]);
            delete_gate_bootstrapping_ciphertext_array(w_pos, end_j[i][j]);
        }
    }
#endif
}
