#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include "gurobi_c++.h"
#include "equation.h"

using namespace std;
#define uint32 unsigned int

// ������ȡ�ļ�
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


// �����̻�32λ���ӵ����Աƽ��ķ�����
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


// ��������aes��S���еķ��亯����������������������֮��Ĺ�ϵ
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


// ��������aes��S���еķ��亯����������������������֮��Ĺ�ϵ
void Sbox_SQ_Affine(vector<GRBVar>& A, vector<GRBVar>& B, GRBModel& model) {
    int i, j;
    for (i = 0; i < 32; i++) {
        GRBLinExpr Sum;
        for (j = 0; j < 32; j++) {
            Sum += matrix_sq_tran[i][j] * A[j];
        }
        GRBVar s;
        s = model.addVar(0, 32, 0, GRB_INTEGER);
        model.addConstr(Sum - 2 * s - B[i] == 0);
    }
}


// ��������aes��S���еķ����Բ���Ҳ����4��S�������������������Ӧƫ��֮��Ĺ�ϵ
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


// ��������aes��S���еķ����Բ���Ҳ����4��S�������������������Ӧƫ��֮��Ĺ�ϵ
void Linear_SQ_Sbox(vector<vector<int>> equa_sbox, vector<GRBVar>& A, vector<GRBVar>& B, vector<GRBVar>& C, GRBModel& model) {
    int i, j, t;
    for (i = 0; i < equa_sbox.size(); i++) {
        for (t = 0; t < 4; t++) {
            GRBLinExpr Sum;
            for (j = 0; j < 8; j++) {
                Sum += equa_sbox[i][j] * A[j + 8 * t];
                Sum += equa_sbox[i][j + 8] * B[j + 8 * t];
            }
            for (j = 0; j < 5; j++) {
                Sum += equa_sbox[i][j + 16] * C[j + 5 * t];
            }
            model.addConstr(Sum + equa_sbox[i][21] >= 0);
        }
    }
}


// ��������aes��S�к���S���еķ����Բ���Ҳ����4��S�е�ƫ��Ͷ�Ӧ����2Ϊ�׵Ķ����Ĺ�ϵ
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


// ��������sq����S���еķ����Բ���Ҳ����4��S�е�ƫ��Ͷ�Ӧ����2Ϊ�׵Ķ����Ĺ�ϵ
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


void XOR(GRBVar& x, GRBVar& y, GRBVar& z, GRBModel& model) {
    model.addConstr(-x + y + z >= 0);
    model.addConstr(x - y + z >= 0);
    model.addConstr(x + y - z >= 0);
    model.addConstr(-x - y - z >= -2);
}


//  ����MILPģ��
vector<vector<int>> MILP(vector<vector<int>> ineq_aes, vector<vector<int>> ineq_sq) {
    GRBEnv env = GRBEnv();
    env.set(GRB_IntParam_LogToConsole, 1);    //    print log or not
    GRBModel model = GRBModel(env);
    model.set(GRB_DoubleParam_TimeLimit, 1000000.00);    //    set time    return 0;

    int i, j;

    //  ȫ������ı�����ÿ�����붼��һ��32���صı���
    vector<vector<GRBVar>> MaskVar(20, (vector<GRBVar>)32);

    for (i = 0; i < MaskVar.size(); i++) {
        ostringstream si;   //  ��һ������Ϊ�˸�������������LP�ļ���ʱ��ÿ��ü��
        si << i;
        for (j = 0; j < MaskVar[i].size(); j++) {
            ostringstream s;
            s << j;
            MaskVar[i][j] = model.addVar(0, 1, 0, GRB_BINARY, "M_" + si.str() + "[" + s.str() + "]");
        }
    }

    //  �������Ӷ�Ӧ�����Զȵı�����ÿ�����붼��һ��32���صı���
    vector<vector<GRBVar>> CorVar1(5, (vector<GRBVar>)32);

    for (i = 0; i < CorVar1.size(); i++) {
        ostringstream si;   //  ��һ������Ϊ�˸�������������LP�ļ���ʱ��ÿ��ü��
        si << i;
        for (j = 0; j < CorVar1[i].size(); j++) {
            ostringstream s;
            s << j;
            CorVar1[i][j] = model.addVar(0, 1, 0, GRB_BINARY, "S1_" + si.str() + "[" + s.str() + "]");
        }
    }

    //  ����aes��S�к���S�ж�Ӧ�����Զȵı�����ÿ�����붼��һ��16���صı���
    vector<vector<GRBVar>> CorVar2(2, (vector<GRBVar>)36);

    for (i = 0; i < CorVar2.size(); i++) {
        ostringstream si;   //  ��һ������Ϊ�˸�������������LP�ļ���ʱ��ÿ��ü��
        si << i;
        for (j = 0; j < CorVar2[i].size(); j++) {
            ostringstream s;
            s << j;
            CorVar2[i][j] = model.addVar(0, 1, 0, GRB_BINARY, "S2_" + si.str() + "[" + s.str() + "]");
        }
    }

    //  ����aes��S�к���S�ж�Ӧ�����ԶȵĶ���������CorVar2�ı���ת��Ϊ2��ָ������ʽ���ı�����ÿ�����붼��һ��36���صı���
    vector<vector<GRBVar>> CorVar3(1, (vector<GRBVar>)20);

    for (i = 0; i < CorVar3.size(); i++) {
        ostringstream si;   //  ��һ������Ϊ�˸�������������LP�ļ���ʱ��ÿ��ü��
        si << i;
        for (j = 0; j < CorVar3[i].size(); j++) {
            ostringstream s;
            s << j;
            CorVar3[i][j] = model.addVar(0, 1, 0, GRB_BINARY, "S3_" + si.str() + "[" + s.str() + "]");
        }
    }

    cout << "��Ԫ��ʼ�����" << endl;

    // Ҫ��ÿ���������������Ϊ0
    GRBLinExpr Sum, Sum1, Sum2;
    for (j = 0; j < 32; j++) {
        Sum += MaskVar[0][j];
    }
    model.addConstr(Sum >= 1);
    for (j = 0; j < 32; j++) {
        Sum1 += MaskVar[7][j];
    }
    model.addConstr(Sum1 >= 1);
    for (j = 0; j < 32; j++) {
        Sum2 += MaskVar[18][j];
    }
    model.addConstr(Sum2 >= 1);

    // ���ȶ�z_{t+2}�������Աƽ���z_{t+2}������ΪMaskVar[0]
    Add32(MaskVar[0], MaskVar[1], MaskVar[2], CorVar1[0], model);
    Add32(MaskVar[1], MaskVar[3], MaskVar[4], CorVar1[1], model);
    Sbox_AES_Affine(MaskVar[0], MaskVar[5], model);
    Linear_AES_Sbox(ineq_aes, MaskVar[6], MaskVar[5], CorVar2[0], model);

	// ���Ŷ�z_{t+1}�������Աƽ���z_{t+1}������ΪMaskVar[7]
    Add32(MaskVar[7], MaskVar[8], MaskVar[9], CorVar1[2], model);

	// Ȼ�����M[0]z_{t+2}��M[1]z_{t+1}�Ľ�������ҽ��н�һ�����Աƽ�
    for (i = 0; i < MaskVar[0].size(); ++i) {
        XOR(MaskVar[6][i], MaskVar[8][i], MaskVar[10][i], model);
    }
    for (i = 0; i < MaskVar[0].size(); ++i) {
        XOR(MaskVar[4][i], MaskVar[7][i], MaskVar[11][i], model);
    }
    Add32(MaskVar[10], MaskVar[12], MaskVar[13], CorVar1[3], model);
    Sbox_AES_Affine(MaskVar[11], MaskVar[14], model);
    Linear_AES_Sbox(ineq_aes, MaskVar[15], MaskVar[14], CorVar2[1], model);
    Sbox_SQ_Affine(MaskVar[3], MaskVar[16], model);
    Linear_SQ_Sbox(ineq_sq, MaskVar[17], MaskVar[16], CorVar3[0], model);
    // ��Ҫ����R3_t������Ϊ0
    for (i = 0; i < MaskVar[0].size(); ++i) {
        model.addConstr(MaskVar[12][i] == 0);
    }
    for (i = 0; i < MaskVar[0].size(); ++i) {
        XOR(MaskVar[13][i], MaskVar[17][i], MaskVar[18][i], model);
    }

	// ���Ŷ�z_{t}�������Աƽ���z_{t}������ΪMaskVar[18]
    Add32(MaskVar[18], MaskVar[15], MaskVar[19], CorVar1[4], model);

    // ����Ŀ�꺯��
    GRBLinExpr Obj;

    for (i = 0; i < CorVar1.size(); i++) {
        for (j = 0; j < 31; j++) {
            Obj -= 100 * CorVar1[i][j];
        }
    }

    for (i = 0; i < CorVar2.size(); i++) {
        for (j = 0; j < 4; j++) {
            Obj += 100 * CorVar2[i][0 + 9 * j];
            Obj += 200 * CorVar2[i][1 + 9 * j];
            Obj += 258 * CorVar2[i][2 + 9 * j];
            Obj += 300 * CorVar2[i][3 + 9 * j];
            Obj += 332 * CorVar2[i][4 + 9 * j];
            Obj += 358 * CorVar2[i][5 + 9 * j];
            Obj += 381 * CorVar2[i][6 + 9 * j];
            Obj += 400 * CorVar2[i][7 + 9 * j];
            Obj += 700 * CorVar2[i][8 + 9 * j];
        }
    }

    for (i = 0; i < CorVar3.size(); i++) {
        for (j = 0; j < 4; j++) {
            Obj += 300 * CorVar3[i][0 + 5 * j];
            Obj += 400 * CorVar3[i][1 + 5 * j];
            Obj += 458 * CorVar3[i][2 + 5 * j];
            Obj += 500 * CorVar3[i][3 + 5 * j];
            Obj += 700 * CorVar3[i][4 + 5 * j];
        }
    }

    model.setObjective(Obj, GRB_MAXIMIZE);

    model.update();
    model.write("snow3g.lp");    //    write lp file
    model.optimize();
    model.write("snow3g.sol");   //    write sol file

    // ��ȡĿ�꺯����Ŀ��ֵ
    double obj_temp = model.get(GRB_DoubleAttr_ObjVal);

    // ��Ŀ��ֵת��Ϊ�������Ͳ���ֵ�� obj
    int obj = static_cast<int>(obj_temp);

    vector<vector<int>> Sol(20, vector<int>(32, 0));    //  �ռ��⣬��ı�ʾ��Sol[i]����һ������Ԫ�أ���Sol[i] = {j,k}�������j���ֵĵ�k���������Ծ
    //    Get soultions
    if (model.get(GRB_IntAttr_Status) != GRB_INFEASIBLE) {
        for (i = 0; i < MaskVar.size(); i++) {
            for (j = 0; j < MaskVar[i].size(); j++) {
                if (MaskVar[i][j].get(GRB_DoubleAttr_Xn) == 1) {
                    //Sol.push_back({ i,j });
                    Sol[i][j] = 1;
                    cout << "M" << i << "[" << j << "]" << " ";
                }
            }
        }
        cout << endl;
        for (i = 0; i < CorVar1.size(); i++) {
            for (j = 0; j < CorVar1[i].size(); j++) {
                if (CorVar1[i][j].get(GRB_DoubleAttr_Xn) == 1) {
                    //Sol.push_back({ i+17,j });
                    cout << "C1" << i << "[" << j << "]" << " ";
                }
            }
        }
        cout << endl;
        for (i = 0; i < CorVar2.size(); i++) {
            for (j = 0; j < CorVar2[i].size(); j++) {
                if (CorVar2[i][j].get(GRB_DoubleAttr_Xn) == 1) {
                    //Sol.push_back({ i+21,j });
                    cout << "C2" << i << "[" << j << "]" << " ";
                }
            }
        }
        cout << endl;
        for (i = 0; i < CorVar3.size(); i++) {
            for (j = 0; j < CorVar3[i].size(); j++) {
                if (CorVar3[i][j].get(GRB_DoubleAttr_Xn) == 1) {
                    //Sol.push_back({ i+25,j });
                    cout << "C3" << i << "[" << j << "]" << " ";
                }
            }
        }
        cout << endl;
    }

    else printf("No solution found\n");

    return Sol;
}


int main() {
    vector<vector<int>> ineq_aes = ReadDataInt("ineq_aes_25.txt", 10);
    vector<vector<int>> ineq_sq = ReadDataInt("ineq_sq_21.txt", 10);
    cout << "����ʽ��ȡ���" << endl;

    vector<vector<int>> mask;

    mask = MILP(ineq_aes, ineq_sq);

    // ��sol����д���ļ�
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
