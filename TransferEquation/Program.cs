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

        static double[] U(double[] a, double t, double[] X, double L, int n_x)
        {
            double[] Uex = new double[n_x];
            for (int i = 0; i < n_x; i++)
            {
                if (a[i] >= 0)
                {
                    if (X[i] >= a[i] * t)
                    {
                        Uex[i] = U_0(X[i] - a[i] * t, L);
                    }
                    else
                    {
                        Uex[i] = Mu(t - X[i] / a[i]);
                    }
                }
                else
                {
                    if (L - X[i] >= Math.Abs(a[i]) * t)
                    {
                        Uex[i] = U_0(L - X[i] - Math.Abs(a[i]) * t, L);
                    }
                    else
                    {
                        Uex[i] = Mu(t - (L - X[i]) / Math.Abs(a[i]));
                    }
                }
            }
            return Uex;
        }

        static double A(double t, double x)
        {
            return Math.Sin(t / 4.5 - x / 2);
        }

        static double[] A(double t, double h_x, int n_x)
        {
            double[] a = new double[n_x]; 
            for (int i = 0; i < n_x; i++)
            {
                a[i] = A(t, h_x * i);
            }
            return a;
        }

        static double Ht(double[] a, double h_x)
        {
            double max_a = 0;
            double abs_a;
            for (int i = 0; i < a.Length; i++)
            {
                abs_a = Math.Abs(a[i]);
                if (abs_a > max_a)
                {
                    max_a = abs_a;
                }
            }
            if (max_a != 0)
            {
                return h_x / max_a;
            }
            else
            {
                return h_x;
            }
        }

        static void Main(string[] args)
        {
            int Nx = 41; // Количество точек
            double L = 10;
            double Hx = L / (Nx - 1);
            double t = 0;
            double[] X = new double[Nx];
            for (int i = 0; i < Nx; i++)
            {
                X[i] = i * Hx;
            }
            double[] Uj = new double[Nx];

            for (int i = 0; i < Nx; i++)
            {
                Uj[i] = U_0(X[i], L);
            }

            double[] Uex = new double[1];
            double[] Uj1 = new double[Nx];
            double[] a;
            double h_t;
            double gm;
            while (t < L)
            {
                a = A(t, Hx, Nx);
                h_t = Ht(a, Hx);
                t += h_t;
                if (a[0] > 0)
                {
                    Uj1[0] = Mu(t);
                }
                else
                {
                    gm = a[0] * h_t / Hx;
                    Uj1[0] = (1 + gm) * Uj[0] - gm * Uj[1];
                }
                for (int i = 1; i < Nx - 1; i++)
                {
                    gm = a[i] * h_t / Hx;
                    if (a[i] >= 0)
                    {
                        Uj1[i] = (1 - gm) * Uj[i] + gm * Uj[i - 1];
                    }
                    else
                    {
                        Uj1[i] = (1 + gm) * Uj[i] - gm * Uj[i + 1];
                    }
                }
                if (a[^1] < 0)
                {
                    Uj1[^1] = Uj[^1];
                }
                else
                {
                    gm = a[^1] * h_t / Hx;
                    Uj1[^1] = (1 - gm) * Uj[^1] + gm * Uj[^2];
                }
                for (int i = 0; i < Nx; i++)
                {
                    Uj[i] = Uj1[i];
                    Console.Write(Uj[i] + ", ");
                }
                Console.WriteLine();
                Console.WriteLine();
            }

            Console.WriteLine("Приближенное решение");
            for (int i = 0; i < Uj.Length; i++)
            {
                Console.Write(Uj[i] + ", ");
            }
            Console.WriteLine();
            Console.WriteLine();

            //Uex = U(a, t - h_t, X, L, Nx);

            //Console.WriteLine("Точное решение");
            //for (int i = 0; i < Uex.Length; i++)
            //{
            //    Console.Write(Uex[i] + ", ");
            //}
        }
    }
}
