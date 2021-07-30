#ifndef __NON_LINEAR_OPTIMAIZATION_HH__
#define __NON_LINEAR_OPTIMAIZATION_HH__
//#define _CRT_SECURE_NO_WARNINGS  
#define _SCL_SECURE_NO_WARNINGS

#define LAMBDA1 1.0
#define LAMBDA2 1.0
#define LAMBDA3 1.0
#define KERNEL_KIND 2
#define Neighbor_NUM  6
#pragma once
#include "lbfgs.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "VectorOperations.h"
#include <minmax.h>
#include <numeric>
#include <algorithm>
#include <iomanip>
using namespace std;
static int large_time = 0;
extern lbfgsfloatval_t *previous_x;
class Non_Linear_Optimization
{
protected:
    vector<lbfgsfloatval_t> m_variables;
	vector<lbfgsfloatval_t> mp_variables;
	lbfgs_parameter_t m_parameters;
	vector<double> m_energysequence;
	vector<double> m_normsequence;
	bool m_printprogress;
	bool first;
public:
    Non_Linear_Optimization() : m_parameters(_defparam)
    {	
		first = true;
		m_printprogress = true;
    }
	bool& GetPrintProgress()
	{
		return m_printprogress;
	}
	const vector<lbfgsfloatval_t> & GetVariables() const
	{
		return m_variables;
	}
	const vector<lbfgsfloatval_t> & GetPreVariables() const
	{
		return mp_variables;
	}
	//return the number of iterations
	lbfgs_parameter_t& GetParameters()
	{
		return m_parameters;
	}
	virtual double GetTerminationAccuracy() const { return m_parameters.epsilon; }
	const vector<double>& GetEnergySequence() const
	{
		return m_energysequence;
	}
	const vector<double>& GetNormSequence() const
	{
		return m_normsequence;
	}
	vector<double>& GetEnergySequence()
	{
		return m_energysequence;
	}
	vector<double>& GetNormSequence()
	{
		return m_normsequence;
	}
	
    virtual int  run()
    { 
		 /*
            Start the L-BFGS optimization; this will invoke the callback functions
            evaluate() and progress() when necessary.
         */
		if (first)
		{
			m_energysequence.clear();
			m_normsequence.clear();
			first = false;
		}
		lbfgsfloatval_t fx = 0;
		int numOfIterations = 0;
        int ret = lbfgs(m_variables.size(), &m_variables[0], &fx, _evaluate, _progress, this, &m_parameters, &numOfIterations);
		if (ret < 0)
		{
			for (int mpi = 0; mpi < m_variables.size(); mpi++)
			{
				mp_variables.push_back(previous_x[mpi]);
			}

		}
        /* Report the result. */
		if (m_printprogress)
			cerr << "L-BFGS optimization terminated with status code = " << ret << endl;
        
		return ret;
	}
	static int SaveSequenceIntoM(const vector<double>& sequence, const char* style, const char* filename, int base = 0)
	{
		if (sequence.empty())
			return base;
		ofstream out(filename);
		out.setf(ios::fixed); 
		out << "clc;clear;close;\n";
		out << "figure(1);\n";
		out << "X = 1:" << sequence.size() << ";\n";
		out << "Y = [";
		for (int i = 0; i < (int)sequence.size(); ++i)
		{
			out << setprecision(10) << sequence[i] << " ";
		}
		out << "];\n";
		out << "plot(X, Y, \'" << style << "\');\n"; 
		out << "hold on" << endl;
		out << "legend(\'" << "???" << "\');\n";
		double xMin = 1;
		double xMax = sequence.size();
		double yMin = *min_element(sequence.begin(), sequence.end());
		double yMax = *max_element(sequence.begin(), sequence.end());
		out << "xlim([" << xMin << " " << xMax << "]);\n";
		out << "ylim([" << yMin << " " << yMax << "]);\n";
		out << "set(gca,\'FontSize\',15);\n";
		out << "set(gca,\'FontName\',\'TimesNewRoman\');\n";
		out << "width = 500;\n";
		out << "height= 400;\n";
		out << "x = 400;\n";
		out << "y = 150;\n";
		out << "set(figure(1), \'Position\', [x y width height]);\n";
		out.close();
		return sequence.size() + base;
	}	
	
	int SaveEnergySequence(const char* style, const char* filename)
	{
		return Non_Linear_Optimization::SaveSequenceIntoM(m_energysequence, style, filename);
	}
	
	int SaveNormSequence(const char* style, const char* filename)
	{
		return Non_Linear_Optimization::SaveSequenceIntoM(m_normsequence, style, filename);
	}
	
protected:
	virtual int progress(
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
    {
		m_energysequence.push_back(fx);
		m_normsequence.push_back(Length(g, n));
		if (m_printprogress)
		{
			cerr << " Iter#" << k << ": ";
			cerr << " Cost = " << setprecision(10) << m_energysequence.back() << ", ";
			cerr << " GraErr = " << setprecision(10) << m_normsequence.back() << ", ";
			cerr << " Stp = " << step << "\n";
			cerr << "-----------------------------------------\n";

			char findex[10];
			_itoa_s(large_time, findex, 9);
			char before[256] = "result/large_result/optimized_";
			strcat_s(before, findex);
			strcat_s(before, ".obj");
			std::cout<<before;
			ofstream out(before, ios::out | ios::trunc);
			if (out.fail())
				throw "fail to read file";
			for (int i = 0; i < m_variables.size()/3; i++)
			{
				out << "v  " << m_variables[i * 3] << ' ';
				out << m_variables[i * 3 + 1] << ' ';
				out << m_variables[i * 3 + 2] << "\n";
			}
			cout << "one_offset..." << endl;
			out.close();
			large_time++;
		}
		
		return 0;
    }
	virtual lbfgsfloatval_t evaluate(
		const lbfgsfloatval_t *x,
		lbfgsfloatval_t *g,
		const int n,
		const lbfgsfloatval_t step
		) = 0;

    static lbfgsfloatval_t _evaluate(
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
		//if (reinterpret_cast<Non_Linear_Optimization*>(instance)->m_printprogress)
		//	cerr << "\tCall your evaluation function once..........\n";
        return reinterpret_cast<Non_Linear_Optimization*>(instance)->evaluate(x, g, n, step);
    }

	static int _progress(
        void *instance,
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
    {
        return reinterpret_cast<Non_Linear_Optimization*>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
    }   
};

#endif // !__NON_LINEAR_OPTIMAIZATION_HH__
