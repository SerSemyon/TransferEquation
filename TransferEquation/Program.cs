using System.Numerics;

namespace TransferEquation
{
    internal class Program
    {
        static double U_0 (double x, double X)
        {
            return Math.Sqrt(1 - x / X) * Math.Abs(Math.Sin(3 * x));
        }

        static double Mu(double t)
        {
            return 1.0 / (1.0 + t / 100) * Math.Pow(Math.Sin(t / 8), 2);
        }

        static double[] U(double a, double t, double[] X, double L, int n_x)
        {
            double[] Uex = new double[n_x];
            for (int i = 0; i < n_x; i++)
            {
                if (a >= 0)
                {
                    if (X[i] >= a * t)
                    {
                        Uex[i] = U_0(X[i] - a * t, L);
                    } 
                    else
                    {
                        Uex[i] = Mu(t - X[i] / a);
                    }
                }
                else
                {
                    if (L - X[i] >= Math.Abs(a) * t)
                    {
                        Uex[i] = U_0(L - X[i] - Math.Abs(a) * t, L);
                    }
                    else
                    {
                        Uex[i] = Mu(t - (L - X[i]) / Math.Abs(a));
                    }
                }
            }
            return Uex;
        }

        static void Main(string[] args)
        {
            int Nx = 41; // Количество точек
            double L = 10;
            double a = -1;
            double Hx = L / (Nx - 1);
            double Ht = Hx / Math.Abs(a);
            double gm = a * Ht / Hx;
            double t = 0;
            double[] X = new double[Nx];
            for (int i = 0; i < Nx; i++)
            {
                X[i] = i * Hx;
            }
            double[] Uj = new double[Nx];
            if (a > 0)
            {
                for (int i = 0; i < Nx; i++)
                {
                    Uj[i] = U_0(X[i], L);
                }
            }
            else
            {
                for (int i = 0; i < Nx; i++)
                {
                    Uj[i] = U_0(X[Nx - i - 1], L);
                }
            }

            double[] Uex = new double[1];
            double[] Uj1 = new double[Nx];
            while (t <= L)
            {
                Uj1[0] = Mu(t);
                for (int i = 1; i < Nx; i++)
                {
                    if (a >= 0)
                    {
                        Uj1[i] = (1 - gm) * Uj[i] + gm * Uj[i - 1];
                    }
                    else
                    {
                        Uj1[i] = (1 + gm) * Uj[Nx - i - 1] - gm * Uj[Nx - i];
                    }
                }
                for (int i = 0; i < Nx; i++)
                {
                    Uj[i] = Uj1[Nx - i - 1];
                }
                t += Ht;
            }

            Console.WriteLine("Приближенное решение");
            for (int i = 0; i < Uj.Length; i++)
            {
                Console.Write(Uj[i] + ", ");
            }
            Console.WriteLine();
            Console.WriteLine();
            Uex = U(a, t - Ht, X, L, Nx);

            Console.WriteLine("Точное решение");
            for (int i = 0; i < Uex.Length; i++)
            {
                Console.Write(Uex[i] + ", ");
            }
        }
    }
}
