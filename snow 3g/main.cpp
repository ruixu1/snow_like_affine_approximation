#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <set>
#include "gurobi_c++.h"
#include "equation.h"

using namespace std;
#define uint32 unsigned int

// Used to read files
vector<vector<uint32>> ReadData(string filename, int base) {
    //  Can read integers with sign with various length, base can take 10 or 16
    ifstream in(filename);
    if (!in) cout << "No such file!" << endl;
    string s;
    vector<vector<string>> res_str;
    int i, j;
    while (getline(in, s)) {
        string t;
        //  Filter data to store integers (include hexadecimal code), space and x (for 16)
        for (i = 0; i < s.size(); i++) {
            if ((s[i] >= 48 && s[i] <= 57) || s[i] == 32 || s[i] == 120 || s[i] >= 97 && s[i] <= 102 || s[i] == 45) t += s[i];
        }

        j = 0;
        vector<string> tmp;
        for (i = 0; i < t.size(); i++) {
            if (t[i] == 32) {   //  deal with space
                string tmp2 = t.substr(j, i - j + 1);
                tmp.push_back(tmp2);
                j = i + 1;
            }
        }

        if (t[t.size() - 1] != 32) {
            string tmp2 = t.substr(j);
            tmp.push_back(tmp2);
        }
        res_str.push_back(tmp);
    }

    vector<vector<uint32>> result;
    for (i = 0; i < res_str.size(); i++) {
        vector<uint32> tmp;
        for (j = 0; j < res_str[i].size(); j++) {
            const char* nptr = res_str[i][j].c_str();
            char* endptr = NULL;
            errno = 0;
            int val = strtol(nptr, &endptr, base);
            tmp.push_back(val);
        }
        result.push_back(tmp);

    }

    return result;
}


vector<vector<int>> ReadDataInt(string filename, int base) {
    //  Can read integers with sign with various length, base can take 10 or 16
    ifstream in(filename);
    if (!in) cout << "No such file!" << endl;
    string s;
    vector<vector<string>> res_str;
    int i, j;
    while (getline(in, s)) {
        string t;
        //  Filter data to store integers (include hexadecimal code), space and x (for 16)
        for (i = 0; i < s.size(); i++) {
            if ((s[i] >= 48 && s[i] <= 57) || s[i] == 32 || s[i] == 120 || s[i] >= 97 && s[i] <= 102 || s[i] == 45) t += s[i];
        }

        j = 0;
        vector<string> tmp;
        for (i = 0; i < t.size(); i++) {
            if (t[i] == 32) {   //  deal with space
                string tmp2 = t.substr(j, i - j + 1);
                tmp.push_back(tmp2);
                j = i + 1;
            }
        }

        if (t[t.size() - 1] != 32) {
            string tmp2 = t.substr(j);
            tmp.push_back(tmp2);
        }
        res_str.push_back(tmp);
    }

    vector<vector<int>> result;
    for (i = 0; i < res_str.size(); i++) {
        vector<int> tmp;
        for (j = 0; j < res_str[i].size(); j++) {
            const char* nptr = res_str[i][j].c_str();
            char* endptr = NULL;
            errno = 0;
            int val = strtol(nptr, &endptr, base);
            tmp.push_back(val);
        }
        result.push_back(tmp);

    }

    return result;
}


// The system of equations used to describe the linear approximation of 32-bit modulo addition
void Add32(vector<GRBVar>& A, vector<GRBVar>& B, vector<GRBVar>& C, vector<GRBVar>& S, GRBModel& model) {
    int i, j, t;

    for (i = 0; i < 8; i++) {
        model.addConstr(equa_add[i][1] * A[0] + equa_add[i][2] * B[0] + equa_add[i][3] * C[0] + equa_add[i][4] * S[0] + equa_add[i][5] >= 0);
    }
    for (j = 1; j < 32; j++) {
        for (i = 0; i < 8; i++) {
            model.addConstr(equa_add[i][0] * S[j - 1] + equa_add[i][1] * A[j] + equa_add[i][2] * B[j] + equa_add[i][3] * C[j] + equa_add[i][4] * S[j] + equa_add[i][5] >= 0);
        }
    }
}


// Used to describe the relationship between the input mask and output mask of the affine function part in the S-box of aes
void Sbox_AES_Affine(vector<GRBVar>& A, vector<GRBVar>& B, GRBModel& model) {
    int i, j;
    for (i = 0; i < 32; i++) {
        GRBLinExpr Sum;
        for (j = 0; j < 32; j++) {
            Sum += matrix_aes_tran[i][j] * A[j];
        }
        GRBVar s;
        s = model.addVar(0, 32, 0, GRB_INTEGER);
        model.addConstr(Sum - 2 * s - B[i] == 0);
    }
}


// Used to describe the relationship between the input mask and output mask of the affine function part in the inverse S-box of AES
void Sbox_AES_Inv_Affine(vector<GRBVar>& A, vector<GRBVar>& B, GRBModel& model) {
    int i, j;
    for (i = 0; i < 32; i++) {
        GRBLinExpr Sum;
        for (j = 0; j < 32; j++) {
            Sum += matrix_aes_inv_tran[i][j] * A[j];
        }
        GRBVar s;
        s = model.addVar(0, 32, 0, GRB_INTEGER);
        model.addConstr(Sum - 2 * s - B[i] == 0);
    }
}


// Used to describe the relationship between the input mask and output mask of the affine function part in the inverse SQ
void Sbox_SQ_Inv_Affine(vector<GRBVar>& A, vector<GRBVar>& B, GRBModel& model) {
    int i, j;
    for (i = 0; i < 32; i++) {
        GRBLinExpr Sum;
        for (j = 0; j < 32; j++) {
            Sum += matrix_sq_inv_tran[i][j] * A[j];
        }
        GRBVar s;
        s = model.addVar(0, 32, 0, GRB_INTEGER);
        model.addConstr(Sum - 2 * s - B[i] == 0);
    }
}


// Used to describe the nonlinear part of the S-box of AES, that is, the relationship between the corresponding deviations of the input mask and the output mask of the four S-boxes
void Linear_AES_Sbox(vector<vector<int>> equa_sbox_aes, vector<GRBVar>& A, vector<GRBVar>& B, vector<GRBVar>& C, GRBModel& model) {
    int i, j, t;
    for (t = 0; t < 4; ++t) {
        for (i = 0; i < equa_sbox_aes.size(); ++i) {
            GRBLinExpr Sum;
            for (j = 0; j < 8; j++) {
                Sum += equa_sbox_aes[i][19 - j] * A[j + 8 * t];
                Sum += equa_sbox_aes[i][11 - j] * B[j + 8 * t];
            }
            for (j = 0; j < 4; j++) {
                Sum += equa_sbox_aes[i][3 - j] * C[j + 4 * t];
            }
            model.addConstr(Sum + equa_sbox_aes[i][20] >= 0);
        }
    }
}


// Used to describe the nonlinear part of the S-box inverse transform of AES, that is, the relationship between the corresponding deviations of the four S-box input masks and output masks
void Linear_AES_Inv_Sbox(vector<vector<int>> equa_sbox_aes_inv, vector<GRBVar>& A, vector<GRBVar>& B, vector<GRBVar>& C, GRBModel& model) {
    int i, j, t;
    
    for (t = 0; t < 4; ++t) {
        for (i = 0; i < equa_sbox_aes_inv.size(); ++i) {
            GRBLinExpr Sum;
            for (j = 0; j < 8; j++) {
                Sum += equa_sbox_aes_inv[i][19 - j] * A[j + 8 * t];
                Sum += equa_sbox_aes_inv[i][11 - j] * B[j + 8 * t];
            }
            for (j = 0; j < 4; j++) {
                Sum += equa_sbox_aes_inv[i][3 - j] * C[j + 4 * t];
            }
            model.addConstr(Sum + equa_sbox_aes_inv[i][20] >= 0);
        }
    }
}


// Used to describe the nonlinear part of the inverse transform of the S-box SQ, that is, the relationship between the corresponding deviations of the four S-box input masks and output masks
void Linear_SQ_Inv_Sbox(vector<vector<int>> equa_sbox_sq_inv, vector<GRBVar>& A, vector<GRBVar>& B, vector<GRBVar>& C, GRBModel& model) {
    int i, j, t;
    for (t = 0; t < 4; ++t) {
        for (i = 0; i < equa_sbox_sq_inv.size(); ++i) {
            GRBLinExpr Sum;
            for (j = 0; j < 8; j++) {
                Sum += equa_sbox_sq_inv[i][18 - j] * A[j + 8 * t];
                Sum += equa_sbox_sq_inv[i][10 - j] * B[j + 8 * t];
            }
            for (j = 0; j < 3; j++) {
                Sum += equa_sbox_sq_inv[i][2 - j] * C[j + 3 * t];
            }
            model.addConstr(Sum + equa_sbox_sq_inv[i][19] >= 0);
        }
    }
}


// Used to describe the nonlinear part of the S-box and inverse S-box of AES, that is, the relationship between the deviation of the four S-boxes and the corresponding logarithm with base 2
void Sbox_AES_Bias_Log(vector<GRBVar>& A, vector<GRBVar>& B, GRBModel& model) {
    int i, j, t;
    for (i = 0; i < 4; i++) {
        for (t = 0; t < 4; t++) {
            GRBLinExpr Sum;
            for (j = 0; j < 4; j++) {
                Sum += equa_aes_bias_log[i][j] * A[j + 4 * t];
            }
            for (j = 0; j < 9; j++) {
                Sum += equa_aes_bias_log[i][j + 4] * B[j + 9 * t];
            }
            model.addConstr(Sum + equa_aes_bias_log[i][13] >= 0);
        }
    }
}


// Used to describe the nonlinear part of the inverse S-box of sq, that is, the relationship between the deviations of the four S-boxes and the corresponding logarithms with base 2
void Sbox_SQ_Bias_Log(vector<GRBVar>& A, vector<GRBVar>& B, GRBModel& model) {
    int i, j, t;
    for (i = 0; i < 2; i++) {
        for (t = 0; t < 4; t++) {
            GRBLinExpr Sum;
            for (j = 0; j < 3; j++) {
                Sum += equa_sq_inv_log[i][j] * A[j + 3 * t];
            }
            for (j = 0; j < 5; j++) {
                Sum += equa_sq_inv_log[i][j + 3] * B[j + 5 * t];
            }
            model.addConstr(Sum + equa_sq_inv_log[i][8] >= 0);
        }
    }
}


//  Building the MILP model
vector<vector<int>> MILP(vector<vector<int>> equa1, vector<vector<int>> equa2, vector<vector<int>> equa3) {
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 1);    //    print log or not
    GRBModel model = GRBModel(env);
    model.set(GRB_DoubleParam_TimeLimit, 100000.00);    //    set time    return 0;

    int i, j;

    //  Variables of all masks, each mask is a 32-bit variable
    vector<vector<GRBVar>> MaskVar(17, (vector<GRBVar>)32);

    for (i = 0; i < MaskVar.size(); i++) {
        ostringstream si;   //   This part is to name the variables, so that you can check them when you look at the LP file.
        si << i;
        for (j = 0; j < MaskVar[i].size(); j++) {
            ostringstream s;
            s << j;
            MaskVar[i][j] = model.addVar(0, 1, 0, GRB_BINARY, "M_" + si.str() + "[" + s.str() + "]");
        }
    }

    //  Describes the linearity of the modulo addition. Each mask is a 32-bit variable.
    vector<vector<GRBVar>> CorVar1(4, (vector<GRBVar>)32);

    for (i = 0; i < CorVar1.size(); i++) {
        ostringstream si;   //   This part is to name the variables, so that you can check them when you look at the LP file.
        si << i;
        for (j = 0; j < CorVar1[i].size(); j++) {
            ostringstream s;
            s << j;
            CorVar1[i][j] = model.addVar(0, 1, 0, GRB_BINARY, "S1_" + si.str() + "[" + s.str() + "]");
        }
    }

    //  Describes the linearity of the S-box and inverse S-box of AES. Each mask is a 16-bit variable.
    vector<vector<GRBVar>> CorVar2(2, (vector<GRBVar>)16);

    for (i = 0; i < CorVar2.size(); i++) {
        ostringstream si;   //   This part is to name the variables, so that you can check them when you look at the LP file.
        si << i;
        for (j = 0; j < CorVar2[i].size(); j++) {
            ostringstream s;
            s << j;
            CorVar2[i][j] = model.addVar(0, 1, 0, GRB_BINARY, "S2_" + si.str() + "[" + s.str() + "]");
        }
    }

    //  A variable describing the logarithm of the linearity corresponding to the S-box and the inverse S-box of AES (i.e., converting the variable of CorVar2 into the form of an exponential of 2). Each mask is a 36-bit variable.
    vector<vector<GRBVar>> CorVar3(2, (vector<GRBVar>)36);

    for (i = 0; i < CorVar3.size(); i++) {
        ostringstream si;   //   This part is to name the variables, so that you can check them when you look at the LP file.
        si << i;
        for (j = 0; j < CorVar3[i].size(); j++) {
            ostringstream s;
            s << j;
            CorVar3[i][j] = model.addVar(0, 1, 0, GRB_BINARY, "S3_" + si.str() + "[" + s.str() + "]");
        }
    }

    //  Variable describing the linearity corresponding to the inverse SQ. Each mask is a 16-bit variable.
    vector<vector<GRBVar>> CorVar4(1, (vector<GRBVar>)12);

    for (i = 0; i < CorVar4.size(); i++) {
        ostringstream si;   //   This part is to name the variables, so that you can check them when you look at the LP file.
        si << i;
        for (j = 0; j < CorVar4[i].size(); j++) {
            ostringstream s;
            s << j;
            CorVar4[i][j] = model.addVar(0, 1, 0, GRB_BINARY, "S4_" + si.str() + "[" + s.str() + "]");
        }
    }

    //  A variable describing the logarithm of the linearity corresponding to the inverse SQ (i.e., converting the CorVar2 variable to the exponential form of 2). Each mask is a 36-bit variable.
    vector<vector<GRBVar>> CorVar5(1, (vector<GRBVar>)20);

    for (i = 0; i < CorVar5.size(); i++) {
        ostringstream si;   //   This part is to name the variables, so that you can check them when you look at the LP file.
        si << i;
        for (j = 0; j < CorVar5[i].size(); j++) {
            ostringstream s;
            s << j;
            CorVar5[i][j] = model.addVar(0, 1, 0, GRB_BINARY, "S5_" + si.str() + "[" + s.str() + "]");
        }
    }

    cout << "Variable initialization completed" << endl;

    // Requires that the weight of some mask is not 0
    int L[3] = { 0,7,10 };
    GRBLinExpr Sum;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 32; j++) {
            Sum += MaskVar[L[i]][j];
        }
    }
    model.addConstr(Sum >= 1);

    //  The system of equations used to describe the linear approximation of 32-bit modulo addition
    Add32(MaskVar[0], MaskVar[1], MaskVar[2], CorVar1[0], model);
    Add32(MaskVar[7], MaskVar[8], MaskVar[9], CorVar1[1], model);
    Add32(MaskVar[10], MaskVar[11], MaskVar[12], CorVar1[2], model);
    Add32(MaskVar[11], MaskVar[14], MaskVar[15], CorVar1[3], model);

    // The equation used to describe the deviation between the input mask and the output mask of the inverse transform of SQ is described
    Linear_SQ_Inv_Sbox(equa1, MaskVar[3], MaskVar[0], CorVar4[0], model);
    // The equation used to describe the deviation between the input mask and the output mask of the inverse transform of AES
    Linear_AES_Inv_Sbox(equa2, MaskVar[4], MaskVar[1], CorVar2[0], model);
    // The equation used to describe the deviation between the input mask and output mask of the AES S-box
    Linear_AES_Sbox(equa3, MaskVar[16], MaskVar[13], CorVar2[1], model);

    // Used to describe the relationship between the input mask and the output mask of the affine function part in the inverse transformation of sq
    Sbox_SQ_Inv_Affine(MaskVar[3], MaskVar[5], model);
    // Used to describe the relationship between the input mask and output mask of the affine function part in the inverse transform of AES
    Sbox_AES_Inv_Affine(MaskVar[4], MaskVar[6], model);
    // Used to describe the relationship between the input mask and output mask of the affine function part in the S-box of AES
    Sbox_AES_Affine(MaskVar[10], MaskVar[13], model);

    // The equation used to describe the deviation of the AES S-box and the corresponding base 2 logarithm
    Sbox_AES_Bias_Log(CorVar2[0], CorVar3[0], model);
    // The equation used to describe the deviation of the inverse transformation of AES and the corresponding base 2 logarithm is described
    Sbox_AES_Bias_Log(CorVar2[1], CorVar3[1], model);
    // The equation used to describe the deviation of the inverse transformation of SQ and the corresponding base 2 logarithm is described
    Sbox_SQ_Bias_Log(CorVar4[0], CorVar5[0], model);

    // The coefficients of R1_t, R2_t, and R3_t are all required to be 0. The formula after adding the three moments is (M8+M16)R1_t + (M6+M7+M15)R2_t + (M5+M14)R3_t + ...
    for (i = 0; i < MaskVar[0].size(); ++i) {
        model.addConstr(MaskVar[8][i] - MaskVar[16][i] == 0);
        model.addConstr(MaskVar[5][i] - MaskVar[14][i] == 0);
        GRBVar s;
        s = model.addVar(0, 2, 0, GRB_INTEGER);
        model.addConstr(MaskVar[6][i] + MaskVar[7][i] + MaskVar[15][i] - 2 * s == 0);
    }


    // Constructing the objective function
    GRBLinExpr Obj;

    for (i = 0; i < CorVar3.size(); i++) {
        for (j = 0; j < 4; j++) {
            Obj += 100 * CorVar3[i][0 + 9 * j];
            Obj += 200 * CorVar3[i][1 + 9 * j];
            Obj += 258 * CorVar3[i][2 + 9 * j];
            Obj += 300 * CorVar3[i][3 + 9 * j];
            Obj += 332 * CorVar3[i][4 + 9 * j];
            Obj += 358 * CorVar3[i][5 + 9 * j];
            Obj += 381 * CorVar3[i][6 + 9 * j];
            Obj += 400 * CorVar3[i][7 + 9 * j];
            Obj += 700 * CorVar3[i][8 + 9 * j];
        }
    }

    for (i = 0; i < CorVar5.size(); i++) {
        for (j = 0; j < 4; j++) {
            Obj += 300 * CorVar5[i][0 + 5 * j];
            Obj += 400 * CorVar5[i][1 + 5 * j];
            Obj += 458 * CorVar5[i][2 + 5 * j];
            Obj += 500 * CorVar5[i][3 + 5 * j];
            Obj += 700 * CorVar5[i][4 + 5 * j];
        }
    }
    for (i = 0; i < CorVar1.size(); i++) {
        for (j = 0; j < 31; j++) {
            Obj -= 100 * CorVar1[i][j];
        }
    }

    model.setObjective(Obj, GRB_MAXIMIZE);

    model.update();
    model.write("snow3g.lp");    //    write lp file
    model.optimize();
    model.write("snow3g.sol");   //    write sol file
    
    vector<int> bias_mod(4, 0), bias_sbox(12, 0);
    vector<vector<int>> Sol(17, vector<int>(32, 0));    //  Collect the solutions of MaskVar
    if (model.get(GRB_IntAttr_Status) != GRB_INFEASIBLE) {
        for (i = 0; i < MaskVar.size(); i++) {
            for (j = 0; j < MaskVar[i].size(); j++) {
                if (MaskVar[i][j].get(GRB_DoubleAttr_Xn) == 1) {
                    Sol[i][j] = 1;
                    cout << "M" << i << "[" << j << "]" << " ";
                }
            }
        }
        cout << endl;
        for (i = 0; i < CorVar1.size(); i++) {
            for (j = 0; j < CorVar1[i].size(); j++) {
                if (CorVar1[i][j].get(GRB_DoubleAttr_Xn) == 1) {
                    cout << "C1" << i << "[" << j << "]" << " ";
                }
            }
            for (j = 0; j < 31; ++j) {
                if (CorVar1[i][j].get(GRB_DoubleAttr_Xn) == 1) {
                    bias_mod[i] += 1;
                }
            }
        }
        cout << endl;
        for (i = 0; i < CorVar2.size(); i++) {
            for (j = 0; j < CorVar2[i].size(); j++) {
                if (CorVar2[i][j].get(GRB_DoubleAttr_Xn) == 1) {
                    cout << "C2" << i << "[" << j << "]" << " ";
                }
            }
        }
        cout << endl;
        for (i = 0; i < CorVar4.size(); i++) {
            for (j = 0; j < CorVar4[i].size(); j++) {
                if (CorVar4[i][j].get(GRB_DoubleAttr_Xn) == 1) {
                    cout << "C4" << i << "[" << j << "]" << " ";
                }
            }
        }
        cout << endl;
    }

    else printf("No solution found\n");

    return Sol;
}


int maxvalue(const vector<int>& seq) {
    int result = seq[0];
    for (int i = 1; i < seq.size(); ++i) {
        if (result < seq[i]) result = seq[i];
    }
    return result;
}


int main() {
    vector<vector<int>> equa_sq_inv = ReadDataInt("ieq.txt", 10);
    vector<vector<int>> equa_aes_inv = ReadDataInt("ieq20.txt", 10);
    vector<vector<int>> equa_aes = ReadDataInt("finalineq.txt", 10);
    cout << "Inequality reading completed" << endl;

    vector<vector<int>> mask;
    mask = MILP(equa_sq_inv, equa_aes_inv, equa_aes);

    // Write the sol vector to a file
    ofstream file;
    file.open("mask.txt");

    for (auto& row : mask) {
        for (size_t i = 0; i < row.size(); ++i) {
            file << row[i];
            if (i != row.size() - 1) {
                file << " "; // Separate values with a space
            }
        }
        file << "\n"; // Newline for each row
    }

    return 0;
}

