#include "eval.h"
#include "omp.h"
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;

const int w_alp = 2;

int A[w_alp] = {0, 0};
int C[w_alp] = {1, 0};
int G[w_alp] = {0, 1};
int T[w_alp] = {1, 1};

string read_str(string file_name) {
    string line;
    ifstream file(file_name);
    if (getline(file, line)) {
        return line;
    } else {
        return NULL;
    };
}

void str_to_bit_array(int **arr, string str) {

    for (int i = 0; i < str.length(); i++) {
        if (str[i] == 'A') {
            arr[i] = A;
        } else if (str[i] == 'C') {
            arr[i] = C;
        } else if (str[i] == 'G') {
            arr[i] = G;
        } else {
            arr[i] = T;
        }
    }
}

void encrypt(int **x, int n, int **y, int m,
             const TFheGateBootstrappingSecretKeySet *key,
             const TFheGateBootstrappingParameterSet *params) {

    LS **cx = new LS *[n];
    LS **cy = new LS *[m];

    FILE *cloud_data = fopen("cloud.data", "wb");

    for (int i = 0; i < n; i++) {
        cx[i] = new_gate_bootstrapping_ciphertext_array(w_alp, params);
        for (int j = 0; j < w_alp; j++) {
            bootsSymEncrypt(&cx[i][j], x[i][j], key);
            export_gate_bootstrapping_ciphertext_toFile(cloud_data, &cx[i][j],
                                                        params);
        }
    }

    for (int i = 0; i < m; i++) {
        cy[i] = new_gate_bootstrapping_ciphertext_array(w_alp, params);
        for (int j = 0; j < w_alp; j++) {
            bootsSymEncrypt(&cy[i][j], y[i][j], key);
            export_gate_bootstrapping_ciphertext_toFile(cloud_data, &cy[i][j],
                                                        params);
        }
    }

    fclose(cloud_data);

    for (int i = 0; i < n; i++) {
        delete_gate_bootstrapping_ciphertext_array(w_alp, cx[i]);
    }
    for (int i = 0; i < m; i++) {
        delete_gate_bootstrapping_ciphertext_array(w_alp, cy[i]);
    }
}

int decrypt_score(LS *a, int w,
                  const TFheGateBootstrappingSecretKeySet *seckey) {

    int mask = -1;
    int val = 0;
    int ai;
    for (int k = 0; k < w; k++) {
        ai = bootsSymDecrypt(&a[k], seckey);
        if (ai == 1) {
            val += 1 << k;
        }
    }

    if (ai == 1) {
        val = val | (mask << w);
    }

    return val;
}

int decrypt_pos(LS *a, int w, const TFheGateBootstrappingSecretKeySet *seckey) {

    int mask = -1;

    int val = 0;
    int ai;
    for (int k = 0; k < w; k++) {
        ai = bootsSymDecrypt(&a[k], seckey);
        if (ai == 1) {
            val += 1 << k;
        }
    }

    return val;
}

int main(int argc, char *argv[]) {
    // OpenMP - Nested Parallelism (Enabled 1, Disabled 0)
    omp_set_nested(1);

    int g_o = -9;
    int g_e = -1;
    int s_m = 5;
    int s_s = -3;

    string x_file = "x.fasta";
    string y_file = "y.fasta";

    if (argc == 3) {
        x_file = argv[1];
        y_file = argv[2];
    }

    if (argc == 7) {
        x_file = argv[1];
        y_file = argv[2];
        s_m = atoi(argv[3]);
        s_s = atoi(argv[4]);
        g_o = atoi(argv[5]);
        g_e = atoi(argv[6]);
    }

    // Reading two strings from files
    string str_x = read_str(x_file);
    string str_y = read_str(y_file);

    int n = str_x.length();
    int m = str_y.length();

    int **x = new int *[n];
    for (int i = 0; i < n; i++) {
        x[i] = new int[w_alp];
    }
    int **y = new int *[m];
    for (int i = 0; i < n; i++) {
        y[i] = new int[w_alp];
    }

    str_to_bit_array(x, str_x);
    str_to_bit_array(y, str_y);

    cout << "Scoring scheme: " << s_m << "/" << s_s << "," << g_o << "," << g_e
         << endl;
    cout << "Input: " << x_file << " " << y_file << endl;
    cout << "Lengths: " << n << " " << m << endl;

    int cell_max_val = min(n, m) * s_m;
    int comp_max_val = max(max(-g_o, -s_s), cell_max_val);
    int w_sco = ceil(log2(comp_max_val + 1)) + 1;
    int w_pos = ceil(log2(max(n, m) + 1));

    // Initialzing security params and generating keys
    const int minimum_lambda = 128;

    TFheGateBootstrappingParameterSet *params =
        new_default_gate_bootstrapping_parameters(minimum_lambda);

    uint32_t seed[] = {314, 1592, 657};
    tfhe_random_generator_setSeed(seed, 3);
    TFheGateBootstrappingSecretKeySet *key =
        new_random_gate_bootstrapping_secret_keyset(params);

    FILE *secret_key = fopen("secret.key", "wb");
    export_tfheGateBootstrappingSecretKeySet_toFile(secret_key, key);
    fclose(secret_key);

    FILE *cloud_key = fopen("cloud.key", "wb");
    export_tfheGateBootstrappingCloudKeySet_toFile(cloud_key, &key->cloud);
    fclose(cloud_key);

    FILE *cloud_key2 = fopen("cloud.key", "rb");
    CK *ck = new_tfheGateBootstrappingCloudKeySet_fromFile(cloud_key2);
    fclose(cloud_key2);

    // Encrypting two strings
    auto start = chrono::steady_clock::now();
    encrypt(x, n, y, m, key, params);
    auto end = chrono::steady_clock::now();

    auto diff = end - start;
    cout << setprecision(2);
    cout << "Encryption time: "
         << chrono::duration<double, milli>(diff).count() / 1000 << "s" << endl;

    FILE *cloud_data = fopen("cloud.data", "rb");

    // Reading cihertexts from file
    LS **cx = new LS *[n];
    LS **cy = new LS *[m];

    for (int i = 0; i < n; i++) {
        cx[i] = new_gate_bootstrapping_ciphertext_array(w_alp, params);
        for (int j = 0; j < w_alp; j++) {
            import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &cx[i][j],
                                                          params);
        }
    }

    for (int i = 0; i < m; i++) {
        cy[i] = new_gate_bootstrapping_ciphertext_array(w_alp, params);
        for (int j = 0; j < w_alp; j++) {
            import_gate_bootstrapping_ciphertext_fromFile(cloud_data, &cy[i][j],
                                                          params);
        }
    }

    fclose(cloud_data);

    // Initializing tables
    LS ***H = new LS **[n + 1];
    LS ***Hx = new LS **[n + 1];
    LS ***Hy = new LS **[n + 1];

    for (int i = 0; i <= n; i++) {
        H[i] = new LS *[m + 1];
        Hx[i] = new LS *[m + 1];
        Hy[i] = new LS *[m + 1];
        for (int j = 0; j <= m; j++) {
            H[i][j] =
                new_gate_bootstrapping_ciphertext_array(w_sco, ck->params);
            Hx[i][j] =
                new_gate_bootstrapping_ciphertext_array(w_pos, ck->params);
            Hy[i][j] =
                new_gate_bootstrapping_ciphertext_array(w_pos, ck->params);
        }
    }

    // Initializing score and positions
    LS *M = new_gate_bootstrapping_ciphertext_array(w_sco, ck->params);
    LS *start_x = new_gate_bootstrapping_ciphertext_array(w_pos, ck->params);
    LS *start_y = new_gate_bootstrapping_ciphertext_array(w_pos, ck->params);
    LS *end_x = new_gate_bootstrapping_ciphertext_array(w_pos, ck->params);
    LS *end_y = new_gate_bootstrapping_ciphertext_array(w_pos, ck->params);

    cout << "Homomorphically evaluating alignment..." << endl;

    start = chrono::steady_clock::now();
    score(H, Hx, Hy, cx, cy, n, m, s_m, s_s, g_o, g_e, w_sco, w_pos, w_alp, ck);
    find(M, start_x, start_y, end_x, end_y, H, Hx, Hy, n, m, w_sco, w_pos, ck);
    end = chrono::steady_clock::now();

    diff = end - start;
    cout << "Homomorphic computation time: "
         << chrono::duration<double, milli>(diff).count() / 1000 << "s" << endl;

    start = chrono::steady_clock::now();
    int score_decrypted = decrypt_score(M, w_sco, key);

#ifdef POSITION
    int start_x_decrypted = decrypt_pos(start_x, w_pos, key) + 1;
    int start_y_decrypted = decrypt_pos(start_y, w_pos, key) + 1;
    int end_x_decrypted = decrypt_pos(end_x, w_pos, key);
    int end_y_decrypted = decrypt_pos(end_y, w_pos, key);
#endif
    end = chrono::steady_clock::now();

    diff = end - start;
    cout << "Decryption time: "
         << chrono::duration<double, milli>(diff).count() / 1000 << "s" << endl;
    cout << "RESULT" << endl;
    cout << "------------------" << endl;
    cout << "Score: " << score_decrypted << endl;
#ifdef POSITION
    cout << "Starting pos: " << start_x_decrypted << " " << start_y_decrypted
         << endl;
    cout << "Ending pos: " << end_x_decrypted << " " << end_y_decrypted << endl;
#endif

    for (int i = 0; i < n; i++) {
        delete_gate_bootstrapping_ciphertext_array(w_alp, cx[i]);
    }
    for (int i = 0; i < m; i++) {
        delete_gate_bootstrapping_ciphertext_array(w_alp, cy[i]);
    }

    delete_gate_bootstrapping_ciphertext_array(w_sco, M);
    delete_gate_bootstrapping_ciphertext_array(w_pos, start_x);
    delete_gate_bootstrapping_ciphertext_array(w_pos, start_y);
    delete_gate_bootstrapping_ciphertext_array(w_pos, end_x);
    delete_gate_bootstrapping_ciphertext_array(w_pos, end_y);

    delete_gate_bootstrapping_secret_keyset(key);
    delete_gate_bootstrapping_cloud_keyset(ck);
    delete_gate_bootstrapping_parameters(params);

    return 0;
}
