using System.Windows.Forms.DataVisualization.Charting;

namespace Lab8._3
{
    partial class Form1
    {
        /// <summary>
        ///  Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        ///  Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        ///  Required method for Designer support - do not modify
        ///  the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            ChartArea chartArea1 = new ChartArea();
            Series series1 = new Series();
            panel1 = new Panel();
            intervalUpDawn = new NumericUpDown();
            label3 = new Label();
            Varnput = new TextBox();
            Meaninput = new TextBox();
            numericUpDown1 = new NumericUpDown();
            button1 = new Button();
            label6 = new Label();
            label2 = new Label();
            label1 = new Label();
            chart = new Chart();
            label7 = new Label();
            label8 = new Label();
            Mean = new Label();
            Dispersion = new Label();
            label10 = new Label();
            label11 = new Label();
            ChiSquare = new Label();
            HypRes = new Label();
            panel1.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)intervalUpDawn).BeginInit();
            ((System.ComponentModel.ISupportInitialize)numericUpDown1).BeginInit();
            ((System.ComponentModel.ISupportInitialize)chart).BeginInit();
            SuspendLayout();
            // 
            // panel1
            // 
            panel1.BackColor = Color.LightBlue;
            panel1.Controls.Add(intervalUpDawn);
            panel1.Controls.Add(label3);
            panel1.Controls.Add(Varnput);
            panel1.Controls.Add(Meaninput);
            panel1.Controls.Add(numericUpDown1);
            panel1.Controls.Add(button1);
            panel1.Controls.Add(label6);
            panel1.Controls.Add(label2);
            panel1.Controls.Add(label1);
            panel1.Location = new Point(35, 28);
            panel1.Margin = new Padding(3, 2, 3, 2);
            panel1.Name = "panel1";
            panel1.Size = new Size(204, 159);
            panel1.TabIndex = 0;
            panel1.Paint += panel1_Paint;
            // 
            // intervalUpDawn
            // 
            intervalUpDawn.Location = new Point(118, 91);
            intervalUpDawn.Margin = new Padding(3, 2, 3, 2);
            intervalUpDawn.Maximum = new decimal(new int[] { 1000000, 0, 0, 0 });
            intervalUpDawn.Name = "intervalUpDawn";
            intervalUpDawn.Size = new Size(83, 23);
            intervalUpDawn.TabIndex = 15;
            intervalUpDawn.Value = new decimal(new int[] { 5, 0, 0, 0 });
            // 
            // label3
            // 
            label3.AutoSize = true;
            label3.Location = new Point(10, 91);
            label3.Name = "label3";
            label3.Size = new Size(51, 15);
            label3.TabIndex = 14;
            label3.Text = "Intervals";
            // 
            // Varnput
            // 
            Varnput.Location = new Point(67, 38);
            Varnput.Margin = new Padding(3, 2, 3, 2);
            Varnput.Name = "Varnput";
            Varnput.Size = new Size(56, 23);
            Varnput.TabIndex = 9;
            Varnput.Text = "1";
            Varnput.TextChanged += textBox_TextChanged;
            // 
            // Meaninput
            // 
            Meaninput.Location = new Point(67, 14);
            Meaninput.Margin = new Padding(3, 2, 3, 2);
            Meaninput.Name = "Meaninput";
            Meaninput.Size = new Size(56, 23);
            Meaninput.TabIndex = 8;
            Meaninput.Text = "0";
            Meaninput.TextChanged += textBox_TextChanged;
            // 
            // numericUpDown1
            // 
            numericUpDown1.Location = new Point(118, 65);
            numericUpDown1.Margin = new Padding(3, 2, 3, 2);
            numericUpDown1.Maximum = new decimal(new int[] { 1000000, 0, 0, 0 });
            numericUpDown1.Name = "numericUpDown1";
            numericUpDown1.Size = new Size(83, 23);
            numericUpDown1.TabIndex = 7;
            numericUpDown1.Value = new decimal(new int[] { 50, 0, 0, 0 });
            // 
            // button1
            // 
            button1.BackColor = Color.DeepSkyBlue;
            button1.ForeColor = Color.White;
            button1.Location = new Point(67, 125);
            button1.Margin = new Padding(3, 2, 3, 2);
            button1.Name = "button1";
            button1.Size = new Size(82, 22);
            button1.TabIndex = 6;
            button1.Text = "Start";
            button1.UseVisualStyleBackColor = false;
            button1.Click += button1_Click_1;
            // 
            // label6
            // 
            label6.AutoSize = true;
            label6.Location = new Point(10, 65);
            label6.Name = "label6";
            label6.Size = new Size(93, 15);
            label6.TabIndex = 5;
            label6.Text = "Number of trials";
            // 
            // label2
            // 
            label2.AutoSize = true;
            label2.Location = new Point(10, 38);
            label2.Name = "label2";
            label2.Size = new Size(51, 15);
            label2.TabIndex = 1;
            label2.Text = "Variance";
            // 
            // label1
            // 
            label1.AutoSize = true;
            label1.Location = new Point(10, 14);
            label1.Name = "label1";
            label1.Size = new Size(37, 15);
            label1.TabIndex = 0;
            label1.Text = "Mean";
            // 
            // chart
            // 
            chartArea1.Name = "ChartArea1";
            chart.ChartAreas.Add(chartArea1);
            chart.Location = new Point(295, 28);
            chart.Margin = new Padding(3, 2, 3, 2);
            chart.Name = "chart";
            series1.ChartArea = "ChartArea1";
            series1.Name = "Series1";
            chart.Series.Add(series1);
            chart.Size = new Size(380, 290);
            chart.TabIndex = 14;
            chart.Click += chart_Click_1;
            // 
            // label7
            // 
            label7.AutoSize = true;
            label7.Location = new Point(48, 189);
            label7.Name = "label7";
            label7.Size = new Size(50, 15);
            label7.TabIndex = 14;
            label7.Text = "Average";
            // 
            // label8
            // 
            label8.AutoSize = true;
            label8.Location = new Point(48, 221);
            label8.Name = "label8";
            label8.Size = new Size(51, 15);
            label8.TabIndex = 16;
            label8.Text = "Variance";
            // 
            // Mean
            // 
            Mean.AutoSize = true;
            Mean.Location = new Point(107, 189);
            Mean.Name = "Mean";
            Mean.Size = new Size(13, 15);
            Mean.TabIndex = 17;
            Mean.Text = "0";
            // 
            // Dispersion
            // 
            Dispersion.AutoSize = true;
            Dispersion.Location = new Point(107, 221);
            Dispersion.Name = "Dispersion";
            Dispersion.Size = new Size(13, 15);
            Dispersion.TabIndex = 18;
            Dispersion.Text = "0";
            // 
            // label10
            // 
            label10.AutoSize = true;
            label10.Location = new Point(48, 283);
            label10.Name = "label10";
            label10.Size = new Size(42, 15);
            label10.TabIndex = 20;
            label10.Text = "Вывод";
            // 
            // label11
            // 
            label11.AutoSize = true;
            label11.Location = new Point(48, 254);
            label11.Name = "label11";
            label11.Size = new Size(72, 15);
            label11.TabIndex = 21;
            label11.Text = "Кси квадрат";
            // 
            // ChiSquare
            // 
            ChiSquare.AutoSize = true;
            ChiSquare.Location = new Point(126, 254);
            ChiSquare.Name = "ChiSquare";
            ChiSquare.Size = new Size(72, 15);
            ChiSquare.TabIndex = 22;
            ChiSquare.Text = "Кси квадрат";
            // 
            // HypRes
            // 
            HypRes.AutoSize = true;
            HypRes.Location = new Point(96, 283);
            HypRes.Name = "HypRes";
            HypRes.Size = new Size(42, 15);
            HypRes.TabIndex = 27;
            HypRes.Text = "Вывод";
            // 
            // Form1
            // 
            AutoScaleDimensions = new SizeF(7F, 15F);
            AutoScaleMode = AutoScaleMode.Font;
            ClientSize = new Size(700, 331);
            Controls.Add(HypRes);
            Controls.Add(ChiSquare);
            Controls.Add(label11);
            Controls.Add(label10);
            Controls.Add(Dispersion);
            Controls.Add(Mean);
            Controls.Add(label8);
            Controls.Add(label7);
            Controls.Add(panel1);
            Controls.Add(chart);
            Margin = new Padding(3, 2, 3, 2);
            Name = "Form1";
            Text = "Form1";
            panel1.ResumeLayout(false);
            panel1.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)intervalUpDawn).EndInit();
            ((System.ComponentModel.ISupportInitialize)numericUpDown1).EndInit();
            ((System.ComponentModel.ISupportInitialize)chart).EndInit();
            ResumeLayout(false);
            PerformLayout();
        }

        #endregion

        private Panel panel1;
        private Label label2;
        private Label label1;
        private TextBox Varnput;
        private TextBox Meaninput;
        private NumericUpDown numericUpDown1;
        private Button button1;
        private Label label6;
        private Chart chart;
        private Label label7;
        private Label label8;
        private Label Mean;
        private Label Dispersion;
        private Label label10;
        private Label label11;
        private Label ChiSquare;
        private Label HypRes;
        private NumericUpDown intervalUpDawn;
        private Label label3;
    }
}