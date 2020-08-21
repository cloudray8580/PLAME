//
// Created by Cloud on 2020/8/21.
//

#ifndef PLAME_LP_SOLVELP_H
#define PLAME_LP_SOLVELP_H

#include "ilcplex/ilocplex.h" // from additional search path
#include <vector>
#include <cmath>
#include <iostream>

using namespace std;

#define Kh(z) 2.0 / (exp(-z/2) + exp(z/2))

/*
double Khf(double z){
    return 2.0 / (exp(-z/2) + exp(z/2));
}

void TestMacro(){
    double Z = 2.5;
    cout << "Kh calculated by macro: " << Kh(Z) << endl;
    cout << "Kh calculated by function: " << Khf(Z) << endl;
}
 */

/*
 * param z: z
 * param omega: omega
 * param gamma: gamma
 * param results: the linear programming results, x_0 ~ x_N-1, epsilon, objective value
 */
void SolveLP(const vector<double> &z, const vector<double> &omega, const double gamma, vector<double> &results) {
    IloEnv env;
    IloModel model(env);
    IloCplex cplex(model);
    IloObjective obj(env);
    IloNumVarArray vars(env);
    IloRangeArray ranges(env);

    //env.setOut(env.getNullStream());
    //cplex.setOut(env.getNullStream());

    //cplex.setParam(IloCplex::NumericalEmphasis, CPX_ON);
    //cplex.setParam(IloCplex::Param::Advance, 0); // turnoff advanced start
    //cplex.setParam(IloCplex::Param::Preprocessing::Presolve, false); // turnoff presolve

    //cplex.setParam(IloCplex::RootAlg, IloCplex::AutoAlg);
    //cplex.setParam(IloCplex::RootAlg, IloCplex::Primal); // using simplex
    //cplex.setParam(IloCplex::RootAlg, IloCplex::Dual); // using dual simplex
    //cplex.setParam(IloCplex::RootAlg, IloCplex::Network);
    //cplex.setParam(IloCplex::RootAlg, IloCplex::Barrier); // set optimizer used interior point method
    //cplex.setParam(IloCplex::RootAlg, IloCplex::Sifting); // set optimizer used sifting method
    //cplex.setParam(IloCplex::RootAlg, IloCplex::Concurrent);

    int M = z.size();
    int N = omega.size(); // corresponding to each x

    // set variables
    for (int i = 0; i < N; i++) {
        vars.add(IloNumVar(env, 0, INFINITY, ILOFLOAT));  // add variable x
    }
    IloNumVar epsilon(env, 0, INFINITY, ILOFLOAT); // variable epsilon

    // set target expression
    IloExpr target_term(env);
    target_term = vars[0]; // D_0 = 1
    for (int i = 1; i < N; i++) { // D_i = 2 (i != 0)
        target_term += 2 * vars[i];
    }
    target_term += gamma * epsilon; // the last variable

    // set objective as target expression
    obj.setExpr(target_term);
    obj.setSense(IloObjective::Minimize);
    model.add(obj);

    // add constraints
    for (int m = 0; m < M; m++) {
        IloExpr model_term(env);
        for (int i = 0; i < N; i++){
            model_term += vars[i] * cos(omega[i]*z[m]);
        }
        model.add(model_term <= Kh(z[m]) + exp(z[m]/2) * epsilon);
        model.add(model_term >= Kh(z[m]) - exp(z[m]/2) * epsilon);
    }

    // solve
    IloNum starttime = cplex.getTime();
    cplex.solve();
//    try {
//        cplex.solve();
//    }
//    catch (IloException &e) {
//        std::cerr << "iloexception: " << e << endl;
//    }
//    catch (std::exception &e) {
//        std::cerr << "standard exception: " << e.what() << endl;
//    }
//    catch (...) {
//        std::cerr << "some other exception: " << endl;
//    }
    IloNum endtime = cplex.getTime();

    // save model
    cplex.exportModel("/mnt/d/PLAME_LP/model.sav");
    cplex.exportModel("/mnt/d/PLAME_LP/model.lp");

    // retrieve results
    results.clear();
    for (int i = 0; i < N; i++) {
        results.push_back(cplex.getValue(vars[i]));
    }
    results.push_back(cplex.getValue(epsilon));
    double target = cplex.getObjValue();
    results.push_back(target);

    cout << "cplex solve time: " << endtime - starttime << endl;

    env.end();
}

void enumerateGamma(double start_value, double stop_value, double step){

    vector<double> Z = {0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375 ,0.5, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875}; // z
    vector<double> W = {0.0, 1.0, 2.0, 3.0, 4.0}; // omega

    double gamma = start_value;
    vector<double> results;
    while(gamma <= stop_value + 0.001){
        SolveLP(Z, W, gamma, results);

        cout << " = = = = = gamma: " << gamma << " = = = = = " << endl;
        for(int i = 0; i < W.size(); i ++){
            cout << "X[" << i << "]: " << results[i] << endl;
        }
        cout << "epsilon: " << results[W.size()] << endl;
        cout << "obj value: " << results[W.size()+1] << endl;

        gamma += step;
    }
}

#endif //PLAME_LP_SOLVELP_H
