using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TransferEquation;

class UnidirectionalTransferEquationSolver
{
    private double _h_t; // внешний шаг по времени
    private double _h_tau; // внутренний шаг по времни
    private double _h_x; // шаг по пространственной координате
    private double[] _y; // дискретное множество текущих значений системы

    public double H_t { 
        get { return _h_t; }
        set
        {
            _h_t = value;
        }
    } // Тут условие Куранта
    public double H_x { get; }
    public double InputValue { get; set; }
    public double OutputValue 
    {
        get
        {
            return _y[^1];
        }
    }
    
    public UnidirectionalTransferEquationSolver(double[] initialValues, double h_t, double h_x)
    {
        _y = new double[initialValues.Length];
        Array.Copy(initialValues, _y, initialValues.Length);
        _h_t = h_t;
        _h_x = h_x;
    }

    //public double Tick(double inputValue, double a)
    //{
    //    _h_tau = _h_x / Math.Abs(a);
    //    int N = 1;
    //    double lambda;
    //    if (_h_tau >= _h_t)
    //    {
    //        _h_tau = _h_t;
    //    } else
    //    {
    //        N = (int)Math.Ceiling(_h_t / _h_tau);
    //        _h_tau = _h_t / N;
    //    }
    //    for (int t = 0; t < N; t++)
    //    {
    //        lambda = a * _h_tau / _h_x;
    //        for (int i = _y.Length - 1; i > 0; i--)
    //        {
    //            _y[i] = (1.0 - lambda) * _y[i] + lambda * _y[i - 1];
    //        }
    //        _y[0] = (1.0 - lambda) * _y[0] + lambda * inputValue;
    //    }
    //}
}
