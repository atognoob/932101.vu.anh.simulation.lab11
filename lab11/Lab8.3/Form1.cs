using MathNet.Numerics.Distributions;
using MathNet.Numerics.Statistics;

namespace Lab8._3
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        private Random rand = new Random();
        private const double alpha = 0.05;
        private void textBox_TextChanged(object sender, EventArgs e)
        {
            button1.Enabled = true;
            
        }

        private void button1_Click(object sender, EventArgs e)
        {
            double a1;
            double a2;
            int N = (int)numericUpDown1.Value;
            int K = (int)intervalUpDawn.Value;

            if(K > N)
            {
                throw new Exception("");
            }
            double[] normalDistribution = new double[N];
            double mean =   Double.Parse(Meaninput.Text);
            double variance =  double.Parse(Varnput.Text);
            for (int i = 0; i < N; i++)
            {
                a1 = rand.NextDouble();
                a2 = rand.NextDouble();
                normalDistribution[i] = ( mean + variance * Math.Sqrt(-2 * Math.Log(a1)) * Math.Sin(2 * Math.PI * a2));
            }

            double min = normalDistribution.Min();
            double max = normalDistribution.Max();

            int[] statistic = new int[K];

            double step = (max - min) / K;

            for(int i = 0; i < N; i++)
            {
                int j = 0;
                while (j < K && normalDistribution[i] >= min + step * j)
                {
                    j++;
                }
                statistic[j - 1]++;
            }

            chart.Series[0].Points.Clear();
            for (int j = 0; j < K; j++)
            {
                chart.Series[0].Points.AddXY(j, statistic[j]);
            }

            double[] freqs = statistic.Select(interval => (double)interval / (double)N).ToArray(); 

            double eMean = getMean(freqs, step, min);
            double eVar = getVar(freqs, mean, step, min);
            Mean.Text = $"{eMean} ({(Math.Abs(eMean - mean) / mean) * 100}%)";
            Dispersion.Text = $"{eVar} ({(Math.Abs(eVar - variance) / variance) * 100}%)";

            double[] probs = new double[K];

            for (int i = 0; i < K; i++)
            {
                double a = min + i * step;
                double b = min + (i + 1) * step;
                probs[i] = (b - a) * NormalPDF((b + a) / 2, mean, variance);
            }

            double chiSquare = getChiSquare(statistic, probs, N);
            ChiSquare.Text = chiSquare.ToString();

            int df = N - 1;
            double tableValue = ChiSquared.InvCDF(df, 1 - alpha);
            bool isAccepted = chiSquare <= tableValue;
            HypRes.Text = isAccepted ? "TRUE" : "FALSE";
        }

        static double NormalPDF(double x, double mean, double stdDev)
        {
            double coefficient = 1.0 / (Math.Sqrt(2 * Math.PI) * stdDev);
            double exponent = Math.Exp(-(Math.Pow(x - mean, 2)) / (2 * Math.Pow(stdDev, 2)));

            return coefficient * exponent;
        }

        private double getChiSquare(int[] observed, double[] expected, int N)
        {
            double chiSquare = 0;

            for (int i = 0; i < observed.Length; i++)
            {
                if (expected[i] > 0)
                    chiSquare += (observed[i] - (expected[i] * N)) * (observed[i] - (expected[i] * N)) / (expected[i] * N);
            }

            return chiSquare;
        }

        private double getVar(IList<double> freqs, double mean, double step, double min = 0)
        {
            double variance = 0;
            for (int i = 0; i < freqs.Count(); i++)
            {
                variance += freqs[i] * Math.Pow(min + step * i, 2);
            }
            return variance - Math.Pow(mean, 2);
        }

        private double getMean(IList<double> freqs, double step, double min = 0)
        {
            double sum = 0;
            for (int i = 0; i < freqs.Count(); i++)
            {
                sum += (double)freqs[i] * (min + i * step);
            }

            return sum;
        }

        private void chart_Click(object sender, EventArgs e)
        {

        }

        private void Form1_Load(object sender, EventArgs e)
        {

        }

        private void chart_Click_1(object sender, EventArgs e)
        {

        }

        private void button1_Click_1(object sender, EventArgs e)
        {
            button1_Click(sender, e);
        }

        private void prob1input_TextChanged(object sender, EventArgs e)
        {

        }

        private void panel1_Paint(object sender, PaintEventArgs e)
        {

        }
    }
}