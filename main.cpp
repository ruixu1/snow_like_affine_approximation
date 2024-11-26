#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
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


// Used to describe the relationship between the input mask and output mask of the affine function part in the S-box of snow2.0
void Sbox_Affine(vector<GRBVar>& A, vector<GRBVar>& B, GRBModel& model) {
    int i, j;
    for (i = 0; i < 32; i++) {
        GRBLinExpr Sum;
        for (j = 0; j < 32; j++) {
            Sum += matrix_affine_transpose[i][j] * A[j];
        }
        GRBVar s;
        s = model.addVar(0, 32, 0, GRB_INTEGER);
        model.addConstr(Sum - 2 * s - B[i] == 0);
    }
}


// Used to describe the nonlinear part of the S-box of Snow2.0, that is, the relationship between the corresponding bias of the input mask and output mask of the four S-boxes
void Linear_AES_Sbox(vector<vector<int>> equa_sbox, vector<GRBVar>& A, vector<GRBVar>& B, vector<GRBVar>& C, GRBModel& model) {
    int i, j, t;
    for (i = 0; i < equa_sbox.size(); i++) {
        for (t = 0; t < 4; t++) {
            GRBLinExpr Sum;
            for (j = 0; j < 8; j++) {
                Sum += equa_sbox[i][19 - j] * A[j + 8 * t];
                Sum += equa_sbox[i][11 - j] * B[j + 8 * t];
            }
            for (j = 0; j < 4; j++) {
                Sum += equa_sbox[i][3 - j] * C[j + 4 * t];
            }
            model.addConstr(Sum + equa_sbox[i][20] >= 0);
        }
    }
}



// Used to describe the nonlinear part of the S-box of snow2.0, that is, the relationship between the bias of the four S-boxes and the corresponding logarithm with base 2
void AES_Sbox_Bias_Log(vector<GRBVar>& A, vector<GRBVar>& B, GRBModel& model) {
    int i, j, t;
    for (i = 0; i < 4; i++) {
        for (t = 0; t < 4; t++) {
            GRBLinExpr Sum;
            for (j = 0; j < 4; j++) {
                Sum += equa_sbox_bias_log[i][j] * A[j + 4 * t];
            }
            for (j = 0; j < 9; j++) {
                Sum += equa_sbox_bias_log[i][j + 4] * B[j + 9 * t];
            }
            model.addConstr(Sum + equa_sbox_bias_log[i][13] >= 0);
        }
    }
}


//  Building the MILP model
vector<vector<int>> MILP(vector<vector<int>> ineq_equa) {
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 1);    //    print log or not
    GRBModel model = GRBModel(env);
    model.set(GRB_DoubleParam_TimeLimit, 100000.00);    //    set time    return 0;

    int i, j;

    //  Variables of all masks, each mask is a 32-bit variable
    vector<vector<GRBVar>> MaskVar(8, (vector<GRBVar>)32);

    for (i = 0; i < MaskVar.size(); i++) {
        ostringstream si;   //  This part is to name the variables, so that you can check them when you look at the LP file.
        si << i;
        for (j = 0; j < MaskVar[i].size(); j++) {
            ostringstream s;
            s << j;
            MaskVar[i][j] = model.addVar(0, 1, 0, GRB_BINARY, "M_" + si.str() + "[" + s.str() + "]");
        }
    }

    //  Describes the linearity of the modulo addition. Each mask is a 32-bit variable.
    vector<vector<GRBVar>> CorVar1(3, (vector<GRBVar>)32);

    for (i = 0; i < CorVar1.size(); i++) {
        ostringstream si;   // This part is to name the variables, so that you can check them when you look at the LP file.
        si << i;
        for (j = 0; j < CorVar1[i].size(); j++) {
            ostringstream s;
            s << j;
            CorVar1[i][j] = model.addVar(0, 1, 0, GRB_BINARY, "S1_" + si.str() + "[" + s.str() + "]");
        }
    }

    //  Variable describing the linearity corresponding to the S-box. Each mask is a 16-bit variable.
    vector<vector<GRBVar>> CorVar2(1, (vector<GRBVar>)16);

    for (i = 0; i < CorVar2.size(); i++) {
        ostringstream si;   //  This part is to name the variables, so that you can check them when you look at the LP file.
        si << i;
        for (j = 0; j < CorVar2[i].size(); j++) {
            ostringstream s;
            s << j;
            CorVar2[i][j] = model.addVar(0, 1, 0, GRB_BINARY, "S2_" + si.str() + "[" + s.str() + "]");
        }
    }

    //  A variable that describes the logarithm of the linearity corresponding to the S-box (i.e., converting the CorVar2 variable to the exponential form of 2). Each mask is a 36-bit variable.
    vector<vector<GRBVar>> CorVar3(1, (vector<GRBVar>)36);

    for (i = 0; i < CorVar3.size(); i++) {
        ostringstream si;   //  This part is to name the variables, so that you can check them when you look at the LP file.
        si << i;
        for (j = 0; j < CorVar3[i].size(); j++) {
            ostringstream s;
            s << j;
            CorVar3[i][j] = model.addVar(0, 1, 0, GRB_BINARY, "S3_" + si.str() + "[" + s.str() + "]");
        }
    }

    cout << "Variable initialization completed" << endl;

    // Requires that the weight of each mask is not 0
    for (i = 0; i < 7; i++) {
        GRBLinExpr Sum;
        for (j = 0; j < 32; j++) {
            Sum += MaskVar[i][j];
        }
        model.addConstr(Sum >= 1);
    }

    //  The system of equations used to describe the linear approximation of 32-bit modulo addition
    Add32(MaskVar[0], MaskVar[1], MaskVar[2], CorVar1[0], model);
    Add32(MaskVar[2], MaskVar[3], MaskVar[4], CorVar1[1], model);
    Add32(MaskVar[4], MaskVar[5], MaskVar[6], CorVar1[2], model);

    // The equation used to describe the bias between the input mask and output mask of the AES S-box
    Linear_AES_Sbox(ineq_equa, MaskVar[6], MaskVar[7], CorVar2[0], model);

    // Used to describe the relationship between the input mask and output mask of the affine function part in the S-box of snow2.0
    Sbox_Affine(MaskVar[0], MaskVar[7], model);

    // The equation used to describe the deviation of the AES S-box and the corresponding base 2 logarithm
    AES_Sbox_Bias_Log(CorVar2[0], CorVar3[0], model);

    // Constructing the objective function
    GRBLinExpr Obj;

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 31; j++) {
            Obj -= 100 * CorVar1[i][j];
        }
    }

    for (j = 0; j < 4; j++) {
        Obj += 100 * CorVar3[0][0 + 9 * j];
        Obj += 200 * CorVar3[0][1 + 9 * j];
        Obj += 258 * CorVar3[0][2 + 9 * j];
        Obj += 300 * CorVar3[0][3 + 9 * j];
        Obj += 332 * CorVar3[0][4 + 9 * j];
        Obj += 358 * CorVar3[0][5 + 9 * j];
        Obj += 381 * CorVar3[0][6 + 9 * j];
        Obj += 400 * CorVar3[0][7 + 9 * j];
        Obj += 700 * CorVar3[0][8 + 9 * j];
    }

    model.setObjective(Obj, GRB_MAXIMIZE);

    model.update();
    model.write("snow2.lp");    //    write lp file
    model.optimize();
    model.write("snow2.sol");   //    write sol file

    vector<vector<int>> Sol(8, vector<int>(32, 0));    //  Collect the solutions of MaskVar
    //    Get soultions
    if (model.get(GRB_IntAttr_Status) != GRB_INFEASIBLE) {
        for (i = 0; i < MaskVar.size(); i++) {
            for (j = 0; j < MaskVar[i].size(); j++) {
                if (MaskVar[i][j].get(GRB_DoubleAttr_Xn) == 1) {
                    Sol[i][j] = 1;
                    cout << i << "[" << j << "]" << " ";
                }
            }
        }
        cout << endl;
    }

    else printf("No solution found\n");

    return Sol;
}


int main() {
    vector<vector<int>> equa = ReadDataInt("finalineq.txt", 10);
    cout << "Inequality reading completed" << endl;

    vector<vector<int>> mask;

    mask = MILP(equa);

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

