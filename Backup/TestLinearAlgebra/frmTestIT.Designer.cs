namespace TestLinearAlgebra
{
	partial class frmTestIt
	{
		/// <summary>
		/// Required designer variable.
		/// </summary>
		private System.ComponentModel.IContainer components = null;

		/// <summary>
		/// Clean up any resources being used.
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
		/// Required method for Designer support - do not modify
		/// the contents of this method with the code editor.
		/// </summary>
		private void InitializeComponent()
		{
      this.components = new System.ComponentModel.Container();
      this.lstMpos = new System.Windows.Forms.ListBox();
      this.lstTP = new System.Windows.Forms.ListBox();
      this.lblMAxisNames = new System.Windows.Forms.Label();
      this.txtMachinePos = new System.Windows.Forms.TextBox();
      this.btnSaveIt = new System.Windows.Forms.Button();
      this.myTips = new System.Windows.Forms.ToolTip(this.components);
      this.groupBox1 = new System.Windows.Forms.GroupBox();
      this.button1 = new System.Windows.Forms.Button();
      this.textBox1 = new System.Windows.Forms.TextBox();
      this.listBox1 = new System.Windows.Forms.ListBox();
      this.groupBox2 = new System.Windows.Forms.GroupBox();
      this.button2 = new System.Windows.Forms.Button();
      this.textBox2 = new System.Windows.Forms.TextBox();
      this.listBox2 = new System.Windows.Forms.ListBox();
      this.chkConnectToCNC = new System.Windows.Forms.CheckBox();
      this.tmrLookAtCNC = new System.Windows.Forms.Timer(this.components);
      this.btnVarTogglePluss = new System.Windows.Forms.Button();
      this.btnVarToggleMinus = new System.Windows.Forms.Button();
      this.btnCompToggleMinus = new System.Windows.Forms.Button();
      this.btnCompTogglePlus = new System.Windows.Forms.Button();
      this.groupBox1.SuspendLayout();
      this.groupBox2.SuspendLayout();
      this.SuspendLayout();
      // 
      // lstMpos
      // 
      this.lstMpos.Font = new System.Drawing.Font("Microsoft Sans Serif", 15.75F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
      this.lstMpos.FormattingEnabled = true;
      this.lstMpos.ItemHeight = 25;
      this.lstMpos.Location = new System.Drawing.Point(1, 71);
      this.lstMpos.Name = "lstMpos";
      this.lstMpos.Size = new System.Drawing.Size(201, 154);
      this.lstMpos.TabIndex = 10;
      // 
      // lstTP
      // 
      this.lstTP.Font = new System.Drawing.Font("Microsoft Sans Serif", 15.75F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
      this.lstTP.FormattingEnabled = true;
      this.lstTP.ItemHeight = 25;
      this.lstTP.Location = new System.Drawing.Point(208, 71);
      this.lstTP.Name = "lstTP";
      this.lstTP.Size = new System.Drawing.Size(201, 154);
      this.lstTP.TabIndex = 11;
      // 
      // lblMAxisNames
      // 
      this.lblMAxisNames.AutoSize = true;
      this.lblMAxisNames.Font = new System.Drawing.Font("Microsoft Sans Serif", 15.75F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
      this.lblMAxisNames.Location = new System.Drawing.Point(-5, 8);
      this.lblMAxisNames.Name = "lblMAxisNames";
      this.lblMAxisNames.Size = new System.Drawing.Size(315, 25);
      this.lblMAxisNames.TabIndex = 12;
      this.lblMAxisNames.Text = "Machine Position (X, Y, Z, A, B)";
      // 
      // txtMachinePos
      // 
      this.txtMachinePos.Font = new System.Drawing.Font("Microsoft Sans Serif", 15.75F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
      this.txtMachinePos.Location = new System.Drawing.Point(0, 36);
      this.txtMachinePos.Name = "txtMachinePos";
      this.txtMachinePos.Size = new System.Drawing.Size(391, 31);
      this.txtMachinePos.TabIndex = 13;
      this.txtMachinePos.Text = "0, 0, 0, 0, 0, 0";
      // 
      // btnSaveIt
      // 
      this.btnSaveIt.Font = new System.Drawing.Font("Microsoft Sans Serif", 12F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(0)));
      this.btnSaveIt.Location = new System.Drawing.Point(156, 231);
      this.btnSaveIt.Name = "btnSaveIt";
      this.btnSaveIt.Size = new System.Drawing.Size(109, 40);
      this.btnSaveIt.TabIndex = 16;
      this.btnSaveIt.Text = "Save Me";
      this.btnSaveIt.UseVisualStyleBackColor = true;
      this.btnSaveIt.Click += new System.EventHandler(this.btnSaveIt_Click);
      // 
      // groupBox1
      // 
      this.groupBox1.Controls.Add(this.btnVarToggleMinus);
      this.groupBox1.Controls.Add(this.btnVarTogglePluss);
      this.groupBox1.Controls.Add(this.button1);
      this.groupBox1.Controls.Add(this.textBox1);
      this.groupBox1.Controls.Add(this.listBox1);
      this.groupBox1.Location = new System.Drawing.Point(463, 23);
      this.groupBox1.Name = "groupBox1";
      this.groupBox1.Size = new System.Drawing.Size(123, 522);
      this.groupBox1.TabIndex = 21;
      this.groupBox1.TabStop = false;
      this.groupBox1.Text = "groupBox1";
      // 
      // button1
      // 
      this.button1.Location = new System.Drawing.Point(18, 409);
      this.button1.Name = "button1";
      this.button1.Size = new System.Drawing.Size(88, 39);
      this.button1.TabIndex = 21;
      this.button1.Text = "Change Value";
      this.button1.UseVisualStyleBackColor = true;
      this.button1.Click += new System.EventHandler(this.button1_Click);
      // 
      // textBox1
      // 
      this.textBox1.Location = new System.Drawing.Point(6, 383);
      this.textBox1.Name = "textBox1";
      this.textBox1.Size = new System.Drawing.Size(100, 20);
      this.textBox1.TabIndex = 20;
      // 
      // listBox1
      // 
      this.listBox1.FormattingEnabled = true;
      this.listBox1.Location = new System.Drawing.Point(21, 21);
      this.listBox1.Name = "listBox1";
      this.listBox1.Size = new System.Drawing.Size(75, 342);
      this.listBox1.TabIndex = 18;
      this.listBox1.SelectedIndexChanged += new System.EventHandler(this.listBox1_SelectedIndexChanged);
      // 
      // groupBox2
      // 
      this.groupBox2.Controls.Add(this.btnCompToggleMinus);
      this.groupBox2.Controls.Add(this.btnCompTogglePlus);
      this.groupBox2.Controls.Add(this.button2);
      this.groupBox2.Controls.Add(this.textBox2);
      this.groupBox2.Controls.Add(this.listBox2);
      this.groupBox2.Location = new System.Drawing.Point(593, 23);
      this.groupBox2.Name = "groupBox2";
      this.groupBox2.Size = new System.Drawing.Size(120, 522);
      this.groupBox2.TabIndex = 22;
      this.groupBox2.TabStop = false;
      this.groupBox2.Text = "groupBox2";
      // 
      // button2
      // 
      this.button2.Location = new System.Drawing.Point(16, 409);
      this.button2.Name = "button2";
      this.button2.Size = new System.Drawing.Size(75, 39);
      this.button2.TabIndex = 23;
      this.button2.Text = "Change Value";
      this.button2.UseVisualStyleBackColor = true;
      this.button2.Click += new System.EventHandler(this.button2_Click);
      // 
      // textBox2
      // 
      this.textBox2.Location = new System.Drawing.Point(6, 383);
      this.textBox2.Name = "textBox2";
      this.textBox2.Size = new System.Drawing.Size(100, 20);
      this.textBox2.TabIndex = 22;
      // 
      // listBox2
      // 
      this.listBox2.FormattingEnabled = true;
      this.listBox2.Location = new System.Drawing.Point(16, 19);
      this.listBox2.Name = "listBox2";
      this.listBox2.Size = new System.Drawing.Size(75, 342);
      this.listBox2.TabIndex = 21;
      this.listBox2.SelectedIndexChanged += new System.EventHandler(this.listBox2_SelectedIndexChanged);
      // 
      // chkConnectToCNC
      // 
      this.chkConnectToCNC.AutoSize = true;
      this.chkConnectToCNC.Location = new System.Drawing.Point(12, 254);
      this.chkConnectToCNC.Name = "chkConnectToCNC";
      this.chkConnectToCNC.Size = new System.Drawing.Size(107, 17);
      this.chkConnectToCNC.TabIndex = 23;
      this.chkConnectToCNC.Text = "Connect To CNC";
      this.chkConnectToCNC.UseVisualStyleBackColor = true;
      this.chkConnectToCNC.CheckedChanged += new System.EventHandler(this.chkConnectToCNC_CheckedChanged);
      // 
      // tmrLookAtCNC
      // 
      this.tmrLookAtCNC.Interval = 500;
      this.tmrLookAtCNC.Tick += new System.EventHandler(this.tmrLookAtCNC_Tick);
      // 
      // btnVarTogglePluss
      // 
      this.btnVarTogglePluss.Location = new System.Drawing.Point(21, 456);
      this.btnVarTogglePluss.Name = "btnVarTogglePluss";
      this.btnVarTogglePluss.Size = new System.Drawing.Size(85, 23);
      this.btnVarTogglePluss.TabIndex = 22;
      this.btnVarTogglePluss.Text = "Toggle Plus";
      this.btnVarTogglePluss.UseVisualStyleBackColor = true;
      this.btnVarTogglePluss.Click += new System.EventHandler(this.btnVarTogglePluss_Click);
      // 
      // btnVarToggleMinus
      // 
      this.btnVarToggleMinus.Location = new System.Drawing.Point(21, 485);
      this.btnVarToggleMinus.Name = "btnVarToggleMinus";
      this.btnVarToggleMinus.Size = new System.Drawing.Size(85, 23);
      this.btnVarToggleMinus.TabIndex = 23;
      this.btnVarToggleMinus.Text = "Toggle Minus";
      this.btnVarToggleMinus.UseVisualStyleBackColor = true;
      this.btnVarToggleMinus.Click += new System.EventHandler(this.btnVarToggleMinus_Click);
      // 
      // btnCompToggleMinus
      // 
      this.btnCompToggleMinus.Location = new System.Drawing.Point(16, 485);
      this.btnCompToggleMinus.Name = "btnCompToggleMinus";
      this.btnCompToggleMinus.Size = new System.Drawing.Size(85, 23);
      this.btnCompToggleMinus.TabIndex = 25;
      this.btnCompToggleMinus.Text = "Toggle Minus";
      this.btnCompToggleMinus.UseVisualStyleBackColor = true;
      this.btnCompToggleMinus.Click += new System.EventHandler(this.btnCompToggleMinus_Click);
      // 
      // btnCompTogglePlus
      // 
      this.btnCompTogglePlus.Location = new System.Drawing.Point(16, 456);
      this.btnCompTogglePlus.Name = "btnCompTogglePlus";
      this.btnCompTogglePlus.Size = new System.Drawing.Size(85, 23);
      this.btnCompTogglePlus.TabIndex = 24;
      this.btnCompTogglePlus.Text = "Toggle Plus";
      this.btnCompTogglePlus.UseVisualStyleBackColor = true;
      this.btnCompTogglePlus.Click += new System.EventHandler(this.btnCompTogglePlus_Click);
      // 
      // frmTestIt
      // 
      this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
      this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
      this.ClientSize = new System.Drawing.Size(858, 557);
      this.Controls.Add(this.chkConnectToCNC);
      this.Controls.Add(this.groupBox2);
      this.Controls.Add(this.groupBox1);
      this.Controls.Add(this.btnSaveIt);
      this.Controls.Add(this.txtMachinePos);
      this.Controls.Add(this.lblMAxisNames);
      this.Controls.Add(this.lstTP);
      this.Controls.Add(this.lstMpos);
      this.Name = "frmTestIt";
      this.Text = "Test Linear Algebra";
      this.Load += new System.EventHandler(this.frmTestIt_Load);
      this.groupBox1.ResumeLayout(false);
      this.groupBox1.PerformLayout();
      this.groupBox2.ResumeLayout(false);
      this.groupBox2.PerformLayout();
      this.ResumeLayout(false);
      this.PerformLayout();

		}

		#endregion

    private System.Windows.Forms.ListBox lstMpos;
    private System.Windows.Forms.ListBox lstTP;
    private System.Windows.Forms.Label lblMAxisNames;
    private System.Windows.Forms.TextBox txtMachinePos;
    private System.Windows.Forms.Button btnSaveIt;
    private System.Windows.Forms.ToolTip myTips;
    private System.Windows.Forms.GroupBox groupBox1;
    private System.Windows.Forms.Button button1;
    private System.Windows.Forms.TextBox textBox1;
    private System.Windows.Forms.ListBox listBox1;
    private System.Windows.Forms.GroupBox groupBox2;
    private System.Windows.Forms.TextBox textBox2;
    private System.Windows.Forms.ListBox listBox2;
    private System.Windows.Forms.Button button2;
    private System.Windows.Forms.CheckBox chkConnectToCNC;
    private System.Windows.Forms.Timer tmrLookAtCNC;
    private System.Windows.Forms.Button btnVarToggleMinus;
    private System.Windows.Forms.Button btnVarTogglePluss;
    private System.Windows.Forms.Button btnCompToggleMinus;
    private System.Windows.Forms.Button btnCompTogglePlus;


  }
}

