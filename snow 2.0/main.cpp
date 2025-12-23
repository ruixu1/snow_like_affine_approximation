#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
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


vector<vector<vector<int>>> Read3DInt(string filename, int base) {
    ifstream in(filename);
    if (!in) {
        cerr << "No such file!" << endl;
        return {};
    }

    vector<vector<vector<int>>> result;  // ��ά���
    vector<vector<int>> current_block;   // ��ǰ��ά block
    string s;

    while (getline(in, s)) {
        // ȥ��ǰ��ո�
        s.erase(0, s.find_first_not_of(" \t\r\n"));
        s.erase(s.find_last_not_of(" \t\r\n") + 1);

        if (s.empty()) {
            // ���б�ʾһ�� block ����
            if (!current_block.empty()) {
                result.push_back(current_block);
                current_block.clear();
            }
            continue;
        }

        // ������ Solution ����
        if (s.find("Solution") != string::npos)
            continue;

        istringstream iss(s);
        vector<int> row;
        string token;

        while (iss >> token) {
            const char* nptr = token.c_str();
            char* endptr = nullptr;
            int val = strtol(nptr, &endptr, base);
            row.push_back(val);
        }

        if (!row.empty())
            current_block.push_back(row);
    }

    // ���һ�� block �ǵüӽ�ȥ
    if (!current_block.empty()) {
        result.push_back(current_block);
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
                Sum += equa_sbox[i][j] * A[j + 8 * t];
                Sum += equa_sbox[i][j + 8] * B[j + 8 * t];
            }
            for (j = 0; j < 9; j++) {
                Sum += equa_sbox[i][j + 16] * C[j + 9 * t];
            }
            model.addConstr(Sum + equa_sbox[i][25] >= 0);
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


void XOR(GRBVar& x, GRBVar& y, GRBVar& z, GRBModel& model) {
    model.addConstr(-x + y + z >= 0);
    model.addConstr(x - y + z >= 0);
    model.addConstr(x + y - z >= 0);
    model.addConstr(-x - y - z >= -2);
}


//  Building the MILP model
vector<vector<int>> MILP(vector<vector<int>> ineq_equa, vector<int> mask_text) {
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

    //  Variable describing the linearity corresponding to the S-box. Each mask is a 36-bit variable.
    vector<vector<GRBVar>> CorVar2(1, vector<GRBVar>(36));
    for (i = 0; i < CorVar2.size(); i++) {
        ostringstream si;
        si << i;
        for (j = 0; j < CorVar2[i].size(); j++) {
            ostringstream s;
            s << j;
            CorVar2[i][j] = model.addVar(0, 1, 0, GRB_BINARY, "S1_" + si.str() + "[" + s.str() + "]");
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
    Sbox_Affine(MaskVar[0], MaskVar[7], model);

    vector<int> index = { 0, 4 };
    for (i = 0; i < index.size(); ++i) {
        for (j = 0; j < 32; ++j) {
            model.addConstr(MaskVar[index[i]][j] - mask_text[j] == 0);
        }
    }

    // Constructing the objective function
    GRBLinExpr Obj;

    for (i = 0; i < 3; i++) {
        for (j = 0; j < 31; j++) {
            Obj -= 100 * CorVar1[i][j];
        }
    }

    for (j = 0; j < 4; j++) {
        Obj += 100 * CorVar2[0][0 + 9 * j];
        Obj += 200 * CorVar2[0][1 + 9 * j];
        Obj += 258 * CorVar2[0][2 + 9 * j];
        Obj += 300 * CorVar2[0][3 + 9 * j];
        Obj += 332 * CorVar2[0][4 + 9 * j];
        Obj += 358 * CorVar2[0][5 + 9 * j];
        Obj += 381 * CorVar2[0][6 + 9 * j];
        Obj += 400 * CorVar2[0][7 + 9 * j];
        Obj += 700 * CorVar2[0][8 + 9 * j];
    }

    model.setObjective(Obj, GRB_MAXIMIZE);

    model.update();
    //model.write("snow2_text.lp");    //    write lp file
    model.optimize();
    //model.write("snow2_text.sol");   //    write sol file

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


// �������е��������Լ����оۼ�
vector<vector<vector<int>>> MILP_max(vector<vector<int>> ineq_equa, int num_solutions) {

    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 1);

    GRBModel model = GRBModel(env);
    model.set(GRB_DoubleParam_TimeLimit, 100000.0);

    // ===== ������� =====
    int i, j;
    vector<vector<GRBVar>> MaskVar(8, vector<GRBVar>(32));
    for (i = 0; i < MaskVar.size(); i++) {
        ostringstream si;
        si << i;
        for (j = 0; j < MaskVar[i].size(); j++) {
            ostringstream s;
            s << j;
            MaskVar[i][j] = model.addVar(0, 1, 0, GRB_BINARY, "M_" + si.str() + "[" + s.str() + "]");
        }
    }

    vector<vector<GRBVar>> CorVar1(3, vector<GRBVar>(32));
    for (i = 0; i < CorVar1.size(); i++) {
        ostringstream si;
        si << i;
        for (j = 0; j < CorVar1[i].size(); j++) {
            ostringstream s;
            s << j;
            CorVar1[i][j] = model.addVar(0, 1, 0, GRB_BINARY, "S0_" + si.str() + "[" + s.str() + "]");
        }
    }

    vector<vector<GRBVar>> CorVar2(1, vector<GRBVar>(36));
    for (i = 0; i < CorVar2.size(); i++) {
        ostringstream si;
        si << i;
        for (j = 0; j < CorVar2[i].size(); j++) {
            ostringstream s;
            s << j;
            CorVar2[i][j] = model.addVar(0, 1, 0, GRB_BINARY, "S1_" + si.str() + "[" + s.str() + "]");
        }
    }

    cout << "Variable initialization completed" << endl;

    // ===== ����Լ�� =====
    Add32(MaskVar[0], MaskVar[1], MaskVar[2], CorVar1[0], model);
    Add32(MaskVar[2], MaskVar[3], MaskVar[4], CorVar1[1], model);
    Add32(MaskVar[4], MaskVar[5], MaskVar[6], CorVar1[2], model);

    Linear_AES_Sbox(ineq_equa, MaskVar[6], MaskVar[7], CorVar2[0], model);
    Sbox_Affine(MaskVar[0], MaskVar[7], model);

    // ===== ����Ŀ�꺯�� =====
    GRBLinExpr Obj;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 31; j++) {
            Obj -= 100 * CorVar1[i][j];
        }
    }
    for (j = 0; j < 4; j++) {
        Obj += 100 * CorVar2[0][0 + 9 * j];
        Obj += 200 * CorVar2[0][1 + 9 * j];
        Obj += 258 * CorVar2[0][2 + 9 * j];
        Obj += 300 * CorVar2[0][3 + 9 * j];
        Obj += 332 * CorVar2[0][4 + 9 * j];
        Obj += 358 * CorVar2[0][5 + 9 * j];
        Obj += 381 * CorVar2[0][6 + 9 * j];
        Obj += 400 * CorVar2[0][7 + 9 * j];
        Obj += 700 * CorVar2[0][8 + 9 * j];
    }
    model.addConstr(Obj == 1258);
    model.setObjective(Obj, GRB_MAXIMIZE);

    // ===== ���ý�ز��� =====
    model.set(GRB_IntParam_PoolSearchMode, 2);   // ���������
    model.set(GRB_IntParam_PoolSolutions, num_solutions); // ���� num_solutions ����
    model.set(GRB_DoubleParam_PoolGap, 1.0);    // ���� gap��������ܶ�Ľ�

    // ===== ��� =====
    model.optimize();

    int solCount = model.get(GRB_IntAttr_SolCount);
    cout << "Found " << solCount << " solutions in solution pool." << endl;

    vector<vector<vector<int>>> all_sol;


    for (int sol = 0; sol < solCount && sol < num_solutions; ++sol) {
        model.set(GRB_IntParam_SolutionNumber, sol);

        vector<vector<int>> curr_sol(MaskVar.size(), vector<int>(MaskVar[0].size(), 0));
        for (i = 0; i < MaskVar.size(); ++i) {
            for (j = 0; j < MaskVar[0].size(); ++j) {
                if (MaskVar[i][j].get(GRB_DoubleAttr_Xn) > 0.5) {
                    curr_sol[i][j] = 1;
                }
            }
        }
        all_sol.push_back(curr_sol);

        vector<int> Sol(MaskVar[6].size(), 0);
        for (int b = 0; b < MaskVar[6].size(); ++b) {
            if (MaskVar[6][b].get(GRB_DoubleAttr_Xn) > 0.5)
                Sol[b] = 1;
        }
    }

    return all_sol;
}


// �������е��������Լ����оۼ�
vector<vector<vector<int>>> MILP_gather(vector<vector<int>> ineq_equa, vector<vector<int>> mask, int num_solutions, vector<double>& AllObjectiveValues) {

    AllObjectiveValues.clear();
    vector<vector<vector<int>>> AllSolutions;

    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 1);

    GRBModel model = GRBModel(env);
    model.set(GRB_DoubleParam_TimeLimit, 100000.0);

    // ===== ������� =====
    int i, j;
    vector<vector<GRBVar>> MaskVar(8, vector<GRBVar>(32));
    for (i = 0; i < MaskVar.size(); i++) {
        ostringstream si;
        si << i;
        for (j = 0; j < MaskVar[i].size(); j++) {
            ostringstream s;
            s << j;
            MaskVar[i][j] = model.addVar(0, 1, 0, GRB_BINARY, "M_" + si.str() + "[" + s.str() + "]");
        }
    }

    vector<vector<GRBVar>> CorVar1(3, vector<GRBVar>(32));
    for (i = 0; i < CorVar1.size(); i++) {
        ostringstream si;
        si << i;
        for (j = 0; j < CorVar1[i].size(); j++) {
            ostringstream s;
            s << j;
            CorVar1[i][j] = model.addVar(0, 1, 0, GRB_BINARY, "S0_" + si.str() + "[" + s.str() + "]");
        }
    }

    vector<vector<GRBVar>> CorVar2(1, vector<GRBVar>(36));
    for (i = 0; i < CorVar2.size(); i++) {
        ostringstream si;
        si << i;
        for (j = 0; j < CorVar2[i].size(); j++) {
            ostringstream s;
            s << j;
            CorVar2[i][j] = model.addVar(0, 1, 0, GRB_BINARY, "S1_" + si.str() + "[" + s.str() + "]");
        }
    }

    cout << "Variable initialization completed" << endl;

    // ===== ����Լ�� =====
    Add32(MaskVar[0], MaskVar[1], MaskVar[2], CorVar1[0], model);
    Add32(MaskVar[2], MaskVar[3], MaskVar[4], CorVar1[1], model);
    Add32(MaskVar[4], MaskVar[5], MaskVar[6], CorVar1[2], model);

    Linear_AES_Sbox(ineq_equa, MaskVar[6], MaskVar[7], CorVar2[0], model);
    Sbox_Affine(MaskVar[0], MaskVar[7], model);

    // �̶���������
    vector<int> index = { 0, 1, 3, 4, 5, 7 };
    for (i = 0; i < index.size(); ++i) {
        for (j = 0; j < 32; ++j) {
            model.addConstr(MaskVar[index[i]][j] - mask[index[i]][j] == 0);
        }
    }

    // ===== ����Ŀ�꺯�� =====
    GRBLinExpr Obj;
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 31; j++) {
            Obj -= 100 * CorVar1[i][j];
        }
    }
    for (j = 0; j < 4; j++) {
        Obj += 100 * CorVar2[0][0 + 9 * j];
        Obj += 200 * CorVar2[0][1 + 9 * j];
        Obj += 258 * CorVar2[0][2 + 9 * j];
        Obj += 300 * CorVar2[0][3 + 9 * j];
        Obj += 332 * CorVar2[0][4 + 9 * j];
        Obj += 358 * CorVar2[0][5 + 9 * j];
        Obj += 381 * CorVar2[0][6 + 9 * j];
        Obj += 400 * CorVar2[0][7 + 9 * j];
        Obj += 700 * CorVar2[0][8 + 9 * j];
    }
    model.setObjective(Obj, GRB_MAXIMIZE);

    // ===== ���ý�ز��� =====
    model.set(GRB_IntParam_PoolSearchMode, 2);   // ���������
    model.set(GRB_IntParam_PoolSolutions, num_solutions); // ���� num_solutions ����
    model.set(GRB_DoubleParam_PoolGap, 1.0);    // ���� gap��������ܶ�Ľ�

    // ===== ��� =====
    model.optimize();

    int solCount = model.get(GRB_IntAttr_SolCount);
    cout << "Found " << solCount << " solutions in solution pool." << endl;

    for (int sol = 0; sol < solCount && sol < num_solutions; ++sol) {
        model.set(GRB_IntParam_SolutionNumber, sol);

        vector<vector<int>> Sol(2, vector<int>(MaskVar[6].size(), 0));
        for (int b = 0; b < MaskVar[6].size(); ++b) {
            if (MaskVar[2][b].get(GRB_DoubleAttr_Xn) > 0.5)
                Sol[0][b] = 1;
            if (MaskVar[6][b].get(GRB_DoubleAttr_Xn) > 0.5)
                Sol[1][b] = 1;
        }

        AllSolutions.push_back(Sol);
        AllObjectiveValues.push_back(model.get(GRB_DoubleAttr_PoolObjVal));

        cout << "Solution " << sol << " objective value: " << model.get(GRB_DoubleAttr_PoolObjVal) << endl;
    }

    return AllSolutions;
}



int main() {
    vector<vector<int>> ineq = ReadDataInt("ineq_aes_25.txt", 10);
    vector<vector<vector<int>>> all_sol = Read3DInt("all_sol1.txt", 10);
    cout << "Inequality reading completed" << endl;

    int i, j;

    vector<vector<vector<int>>> gather_sol;
    mask = all_sol[0];
    int n = 100;
    vector<double> objvalue(n, 0);
    gather_sol = MILP_gather(ineq, mask, n, objvalue);

    for (const auto& elem : objvalue) {
        cout << elem << " ";
    }
    cout << endl;

    const char HEXMAP[] = "0123456789ABCDEF";

    ofstream file;
    file.open("sol_gather1.txt");

    for (auto& sol : gather_sol) {
        for (auto& row : sol) {
            file << "0x";
            for (i = 0; i < 8; ++i) {
                unsigned char byte = 0;
                for (j = 0; j < 4; ++j) {
                    // ��λ��ǰ�����Ե�һ��bit�����λ��bit7��
                    byte = (byte << 1) | (row[i * 4 + j] & 1);
                }
                file << HEXMAP[byte];
                //if (i != 15) file << " ";
            }
            file << " & "; // Newline for each row
        }
        file << "\n"; // Newline for each row
    }

    file.close();

    return 0;
}