/*
备注：1、使用FFT库时，申明fftw_complex* 变量时，需要在执行完成后进行释放内存操作
*/
///

#include "fftw3.h"
#include <complex>
#include "Parameter_Definition.h"
#include "HOS_usercode.h"

//FFTW库工具========================================================================================================
//FFT real to complex
void FFT_1D(double *in, fftw_complex *out, int N)
{
	fftw_plan p;
	p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);
}

//IFFT complex to real
void IFFT_1D(fftw_complex *in, double *out, int N)
{
	fftw_plan p;
	p = fftw_plan_dft_c2r_1d(N, in, out, FFTW_ESTIMATE);
	fftw_execute(p);
	for (int i = 0; i < N; i++)
	{
		out[i] = out[i] / N;
	}
	fftw_destroy_plan(p);
}

//eta与phi_s求导
void Derivation_x(double *in, double *out, double *k)
{
	fftw_complex *FFT_out, *IFFT_in;
	FFT_out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * Num_x);
	IFFT_in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * Num_x);
	FFT_1D(in, FFT_out, Num_x); //执行FFT

	std::complex<double> a(0, 0);
	std::complex<double> b(0, 0);
	std::complex<double> c(0, 0);

	for (int i = 0; i < Num_x; i++)
	{
		a.real(0.0);
		a.imag(k[i]);
		b.real(FFT_out[i][0]);
		b.imag(FFT_out[i][1]);
		c = b * a;
		IFFT_in[i][0] = c.real();
		IFFT_in[i][1] = c.imag();
	}
	IFFT_1D(IFFT_in, out, Num_x); //执行IFFT

	fftw_free(FFT_out); //free
	fftw_free(IFFT_in); //free
}

//对速度势进行z向求导，默认为深水基函数，N为求导阶数
void Derivation_z(double *in, double *out, double *k, int N)
{
	fftw_complex *FFT_out, *IFFT_in;
	FFT_out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * Num_x);
	IFFT_in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * Num_x);
	FFT_1D(in, FFT_out, Num_x); //执行FFT

	std::complex<double> a(0, 0);
	std::complex<double> b(0, 0);
	std::complex<double> c(0, 0);

	for (int i = 0; i < Num_x; i++)
	{
		a.real(pow(abs(k[i]), N)); //abs(kx)
		a.imag(0.0);
		b.real(FFT_out[i][0]);
		b.imag(FFT_out[i][1]);
		c = b * a;
		IFFT_in[i][0] = c.real();
		IFFT_in[i][1] = c.imag();
	}
	IFFT_1D(IFFT_in, out, Num_x); //执行IFFT

	fftw_free(FFT_out); //free
	fftw_free(IFFT_in); //free
}
//=============================================================================================================================
//网格插值函数

//线性插值函数
double linear(double x1, double f1, double x2, double f2, double x)
{
	return (f2 - f1) / (x2 - x1) * (x - x1) + f1;
}

//插值得到eta数值

double Interpolation_Eta(double x, struct HOS_data_reconstruct *HOS_data_rc)
{
	int Index_x = floor((x - X1) / dx_re);
	double x1_eta = HOS_data_rc->X_re[Index_x];
	double f1_eta = HOS_data_rc->Eta[Index_x];
	double x2_eta = HOS_data_rc->X_re[Index_x + 1];
	double f2_eta = HOS_data_rc->Eta[Index_x + 1];

	double eta = linear(x1_eta, f1_eta, x2_eta, f2_eta, x);
	return eta;
}

//插值得到自由面速度
void Interpolation_Fs(double x, double *HOS_Vel_Fs, struct HOS_data_reconstruct *HOS_data_rc)
{
	int Index_x = floor((x - X1) / dx_re);
	double x1_phi_s_x = HOS_data_rc->X_re[Index_x];
	double f1_phi_s_x = HOS_data_rc->Phi_s_x[Index_x];
	double x2_phi_s_x = HOS_data_rc->X_re[Index_x + 1];
	double f2_phi_s_x = HOS_data_rc->Phi_s_x[Index_x + 1];

	HOS_Vel_Fs[0] = linear(x1_phi_s_x, f1_phi_s_x, x2_phi_s_x, f2_phi_s_x, x); //x方向速度

	double x1_phi_s_z = HOS_data_rc->X_re[Index_x];
	double f1_phi_s_z = HOS_data_rc->Phi_s_z[Index_x];
	double x2_phi_s_z = HOS_data_rc->X_re[Index_x + 1];
	double f2_phi_s_z = HOS_data_rc->Phi_s_z[Index_x + 1];

	HOS_Vel_Fs[1] = 0.0;
	HOS_Vel_Fs[2] = linear(x1_phi_s_z, f1_phi_s_z, x2_phi_s_z, f2_phi_s_z, x); //z方向速度
}

//二次线性插值
double Bilinear(double x1, double x2, double y1, double y2, double f1, double f2, double f3, double f4, double x, double y)
{
	double f_R1 = linear(x1, f1, x2, f4, x);
	double f_R2 = linear(x1, f2, x2, f3, x);

	return linear(y1, f_R1, y2, f_R2, y);
}

//快速定位网格位置找到左下角位置
void Find_Index(int *Index, double x, double y, double z)
{

	Index[0] = floor((x - X1) / dx_re);
	if (z >= Z0 && z < Z1)
		Index[1] = floor((z - Z0) / dz1);
	else if (z >= Z1 && z < Z2)
		Index[1] = floor((z - Z1) / dz2) + n1;
	else if (z >= Z2 && z <= Z3)
		Index[1] = floor((z - Z2) / dz3) + n1 + n2;
}

//速度场插值
void Interpolation_Vel(double x, double y, double z, double *HOS_Vel, struct Grid_position *Grid, struct HOS_data_reconstruct *HOS_data_rc)
{
	int Index_x, Index_z;
	Index_x = floor((x - X1) / dx_re);
	if (z >= Z0 && z < Z1)
		Index_z = floor((z - Z0) / dz1);
	else if (z >= Z1 && z < Z2)
		Index_z = floor((z - Z1) / dz2) + n1;
	else if (z >= Z2 && z <= Z3)
		Index_z = floor((z - Z2) / dz3) + n1 + n2;

	double x1 = HOS_data_rc->X_re[Index_x];
	double x2 = HOS_data_rc->X_re[Index_x + 1];
	double y1 = Grid->z[Index_z];
	double y2 = Grid->z[Index_z + 1];

	double U1 = HOS_data_rc->U[Index_x][Index_z];
	double U2 = HOS_data_rc->U[Index_x][Index_z + 1];
	double U3 = HOS_data_rc->U[Index_x + 1][Index_z + 1];
	double U4 = HOS_data_rc->U[Index_x + 1][Index_z];

	double W1 = HOS_data_rc->W[Index_x][Index_z];
	double W2 = HOS_data_rc->W[Index_x][Index_z + 1];
	double W3 = HOS_data_rc->W[Index_x + 1][Index_z + 1];
	double W4 = HOS_data_rc->W[Index_x + 1][Index_z];

	HOS_Vel[0] = Bilinear(x1, x2, y1, y2, U1, U2, U3, U4, x, z);
	HOS_Vel[1] = 0.;
	HOS_Vel[2] = Bilinear(x1, x2, y1, y2, W1, W2, W3, W4, x, z);
}

//=============================================================================================================================
//初始化函数HOS计算网格
void Init_Grid(struct Grid_position *Grid, struct HOS_data_reconstruct *Grid_re)
{
	for (int i = 0; i < Num_x; i++)
	{
		Grid->x[i] = X1 + i * dx;
	}
	for (int i = 0; i < Num_x_re; i++)
	{
		Grid_re->X_re[i] = X1 + i * dx_re;
		//Shift变换Kx_re
		if (i <= Num_x_re / 2 - 1)
			Grid_re->Kx_re[i] = 2 * PI / Lx * (0. + i);
		else
			Grid_re->Kx_re[i] = 2 * PI / Lx * (i - (Num_x_re));
	}
	for (int i = 0; i < Num_z; i++)
	{
		int N12 = n1 + n2;
		int N13 = n1 + n2 + n3;
		if (i < n1)
			Grid->z[i] = Z0 + i * dz1;
		else if ((i >= n1) && (i < N12))
			Grid->z[i] = Z1 + (i - n1) * dz2;
		else if ((i >= N12) && (i <= N13))
			Grid->z[i] = Z2 + (i - (N12)) * dz3;
	}
}

void Initi_HOS_data(struct HOS_parameters *HOS_para, struct Grid_position *Grid,
					struct Wave_parameters *wave_para, double t)
{
	double x = 0;
	double z = 0;
	HOS_para->M = Num_M;
	HOS_para->istep = 0;
	HOS_para->dt = Timestep;
	////读取速度势，读取速度势与波面文件
	//Read_array1D("phi0.txt",HOS_para->Phi_s,Num_x);
	//Read_array1D("eta0.txt", HOS_para->Eta, Num_x);

	for (int i = 0; i < Num_x; i++)
	{
		//初始化波面和速度为五阶stokes
		x = Grid->x[i];
		z = Eta_fifth(x, 0);
		HOS_para->Eta[i] = z;
		HOS_para->Phi_s[i] = Phi_fifth(x, z, 0);

		//Shift变换Kx
		if (i <= Num_x / 2 - 1)
			HOS_para->Kx[i] = 2 * PI / Lx * (0. + i);
		else
			HOS_para->Kx[i] = 2 * PI / Lx * (i - (Num_x));

		//filter_x滤波函数
		if (abs(HOS_para->Kx[i]) > Kmax)
			HOS_para->Filter_x[i] = 0.;
		else
		{
			HOS_para->Filter_x[i] = 1.;
		}

		//Kx预先滤波
		HOS_para->Kx[i] = HOS_para->Kx[i] * HOS_para->Filter_x[i];
	}
}
//初始化中间变量为 0
void Initial_phi_m_z(struct HOS_parameters *HOS_para)
{
	for (int m = 0; m < Num_M; m++)
	{
		for (int i = 0; i < Num_x; i++)
		{
			HOS_para->Phi_m[i][m] = 0.0;
		}
	}
	for (int i = 0; i < Num_x; i++)
	{
		HOS_para->Phi_s_z[i] = 0.;
	}
}
//HOS求解函数=============================================================================================================================
//基函数
double Eigen_function(double k, double z)
{
	double results;
	//if (z > 0) results = 1 + abs(k) * z + 0.5 * pow(abs(k) * z, 2);
	//else results = exp(abs(k) * z);
	results = exp(abs(k) * z);
	return results;
}

//阶乘函数
int Prod(int n)
{
	int count = 1;
	if (n > 0)
	{
		for (int i = 1; i <= n; i++)
			count = count * i;
	}
	else
	{
	}
	return count;
}

//由phi_s计算得到phi_m
void Get_Phi_m(struct HOS_parameters *HOS_pa)
{
	double temp_in[Num_x]{};
	double temp_out[Num_x]{};
	double Phi_m_z_list[Num_x][Num_M][Num_M]{}; //储存速度势各阶的导数

	//清零phi_m与phi_z，避免下次累加错误！！！！！

	Initial_phi_m_z(HOS_pa);

	for (int m = 0; m < Num_M; m++)
	{
		//phi_1=phi_s
		if (m == 0)
		{
			for (int i = 0; i < Num_x; i++)
			{
				double val = HOS_pa->Phi_s[i];
				HOS_pa->Phi_m[i][m] = val;
				temp_in[i] = val;
			}
		}
		else
		{
			for (int k = 0; k < m; k++)
			{
				int prod = Prod(m - k);
				for (int i = 0; i < Num_x; i++)
				{
					HOS_pa->Phi_m[i][m] = HOS_pa->Phi_m[i][m] +
										  (-1. * pow(HOS_pa->Eta[i], m - k) / prod * Phi_m_z_list[i][k][m - k - 1]);
				}
			}
			for (int i = 0; i < Num_x; i++)
				temp_in[i] = HOS_pa->Phi_m[i][m]; //FFT输入数据
		}
		//求摄动展开各项各阶导数至M阶
		for (int n = 0; n < Num_M - m; n++)
		{
			Derivation_z(temp_in, temp_out, HOS_pa->Kx, n + 1); //z方向求各阶导数
			for (int i = 0; i < Num_x; i++)
				Phi_m_z_list[i][m][n] = temp_out[i];
		}
	}
	//累加求Phi_z
	for (int m = 0; m < Num_M; m++)
	{
		for (int n = 0; n < Num_M - m; n++)
		{
			int prod = Prod(n);
			for (int i = 0; i < Num_x; i++)
			{
				HOS_pa->Phi_s_z[i] = HOS_pa->Phi_s_z[i] +
									 Phi_m_z_list[i][m][n] * pow(HOS_pa->Eta[i], n) / prod;
			}
		}
	}
	//求解x方向导数
	Derivation_x(HOS_pa->Phi_s, HOS_pa->Phi_s_dx, HOS_pa->Kx);
}

//计算phiz，从而计算时间导数eta_t与phi_t
void GetDiff(struct HOS_parameters *HOS_pa)
{
	double temp_in[Num_x]{};
	double temp_out[Num_x]{};
	double Phi_m_z_list[Num_x][Num_M][Num_M]{}; //储存速度势各阶的导数

	//清零phi_m与phi_z，避免下次累加错误！！！！！

	Initial_phi_m_z(HOS_pa);

	for (int m = 0; m < Num_M; m++)
	{
		//phi_1=phi_s
		if (m == 0)
		{
			for (int i = 0; i < Num_x; i++)
			{
				double val = HOS_pa->Phi_s[i];
				HOS_pa->Phi_m[i][m] = val;
				temp_in[i] = val;
			}
		}
		else
		{
			for (int k = 0; k < m; k++)
			{
				int prod = Prod(m - k);
				for (int i = 0; i < Num_x; i++)
				{
					HOS_pa->Phi_m[i][m] = HOS_pa->Phi_m[i][m] +
										  (-1. * pow(HOS_pa->Eta[i], m - k) / prod * Phi_m_z_list[i][k][m - k - 1]);
				}
			}
			for (int i = 0; i < Num_x; i++)
				temp_in[i] = HOS_pa->Phi_m[i][m]; //FFT输入数据
		}
		//求摄动展开各项各阶导数至M阶
		for (int n = 0; n < Num_M - m; n++)
		{
			Derivation_z(temp_in, temp_out, HOS_pa->Kx, n + 1); //z方向求各阶导数
			for (int i = 0; i < Num_x; i++)
				Phi_m_z_list[i][m][n] = temp_out[i];
		}
	}

	//累加求Phi_z
	for (int m = 0; m < Num_M; m++)
	{
		for (int n = 0; n < Num_M - m; n++)
		{
			int prod = Prod(n);
			for (int i = 0; i < Num_x; i++)
			{
				HOS_pa->Phi_s_z[i] = HOS_pa->Phi_s_z[i] +
									 Phi_m_z_list[i][m][n] * pow(HOS_pa->Eta[i], n) / prod;
			}
		}
	}
	//求解x方向导数
	Derivation_x(HOS_pa->Eta, HOS_pa->Eta_dx, HOS_pa->Kx);
	Derivation_x(HOS_pa->Phi_s, HOS_pa->Phi_s_dx, HOS_pa->Kx);

	double a, b, c, d;
	for (int i = 0; i < Num_x; i++)
	{
		a = HOS_pa->Phi_s_dx[i] * HOS_pa->Eta_dx[i];
		b = (1 + HOS_pa->Eta_dx[i] * HOS_pa->Eta_dx[i]);
		c = HOS_pa->Phi_s_dx[i] * HOS_pa->Phi_s_dx[i];
		d = HOS_pa->Phi_s_z[i];

		HOS_pa->Eta_dt[i] = -1. * (a - b * d);
		HOS_pa->Phi_s_dt[i] = -1. * (HOS_pa->Eta[i] * g + 0.5 * c - 0.5 * b * pow(d, 2));
	}
}

//四阶龙格库塔方法随着时间推进，需要调用四次GetDiff()函数
void Solver_RK4(struct HOS_parameters *HOS_pa)
{

	double K1_eta[Num_x], K1_phi[Num_x]; //存储中间时刻参数；
	double K2_eta[Num_x], K2_phi[Num_x];
	double K3_eta[Num_x], K3_phi[Num_x];
	double K4_eta[Num_x], K4_phi[Num_x];
	double eta[Num_x], phi_s[Num_x]; //暂时存储Eta，Phi_s；

	double dt = HOS_pa->dt; //时间步dt
	for (int i = 0; i < Num_x; i++)
	{
		eta[i] = HOS_pa->Eta[i];
		phi_s[i] = HOS_pa->Phi_s[i];
	}

	GetDiff(HOS_pa); //1

	for (int i = 0; i < Num_x; i++)
	{
		K1_eta[i] = HOS_pa->Eta_dt[i];
		K1_phi[i] = HOS_pa->Phi_s_dt[i];

		HOS_pa->Eta[i] = eta[i] + 0.5 * dt * K1_eta[i];
		HOS_pa->Phi_s[i] = phi_s[i] + 0.5 * dt * K1_phi[i];
	}

	GetDiff(HOS_pa); //2

	for (int i = 0; i < Num_x; i++)
	{
		K2_eta[i] = HOS_pa->Eta_dt[i];
		K2_phi[i] = HOS_pa->Phi_s_dt[i];

		HOS_pa->Eta[i] = eta[i] + 0.5 * dt * K2_eta[i];
		HOS_pa->Phi_s[i] = phi_s[i] + 0.5 * dt * K2_phi[i];
	}

	GetDiff(HOS_pa); //3

	for (int i = 0; i < Num_x; i++)
	{
		K3_eta[i] = HOS_pa->Eta_dt[i];
		K3_phi[i] = HOS_pa->Phi_s_dt[i];

		HOS_pa->Eta[i] = eta[i] + dt * K3_eta[i];
		HOS_pa->Phi_s[i] = phi_s[i] + dt * K3_phi[i];
	}

	GetDiff(HOS_pa); //4

	for (int i = 0; i < Num_x; i++)
	{
		K4_eta[i] = HOS_pa->Eta_dt[i];
		K4_phi[i] = HOS_pa->Phi_s_dt[i];

		HOS_pa->Eta[i] = eta[i] + dt / 6 *
									  (K1_eta[i] + 2.0 * K2_eta[i] + 2.0 * K3_eta[i] + K4_eta[i]);
		HOS_pa->Phi_s[i] = phi_s[i] + dt / 6 *
										  (K1_phi[i] + 2.0 * K2_phi[i] + 2.0 * K3_phi[i] + K4_phi[i]);
	}

	//滤波
	Derivation_z(HOS_pa->Eta, HOS_pa->Eta, HOS_pa->Filter_x, 1);
	Derivation_z(HOS_pa->Phi_s, HOS_pa->Phi_s, HOS_pa->Filter_x, 1);
	//更新时间步
	HOS_pa->istep += 1;
	/*if (HOS_pa->istep % 50 == 0) {
		printf("RK4求解次数：%d   ", HOS_pa->istep);
		printf("当前求解时间：%f\n", HOS_pa->istep * HOS_pa->dt);
	}*/
}

//=============================================================================================================================
//HOS波浪场重构函数
//FFT补插值函数：in=输入Num_x长度数组 out=输出Num_x_re长度数组
void FFT_interpolation(double *in, double *out)
{
	fftw_complex *FFT_out, *IFFT_in;
	FFT_out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * Num_x);
	IFFT_in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * Num_x_re);
	FFT_1D(in, FFT_out, Num_x); //执行FFT
	//FFT变换结果补零操作
	//IFFT_in存储补零后的数据
	for (int i = 0; i < Num_x_re; i++)
	{
		if (i <= Num_x / 2 - 1)
		{
			IFFT_in[i][0] = FFT_out[i][0];
			IFFT_in[i][1] = FFT_out[i][1];
		}
		else if (i >= (Num_x_re - Num_x / 2)) //shift操作使得对应
		{
			IFFT_in[i][0] = FFT_out[i - Num_x_re + Num_x][0];
			IFFT_in[i][1] = FFT_out[i - Num_x_re + Num_x][1];
		}
		else
		{
			IFFT_in[i][0] = 0;
			IFFT_in[i][1] = 0;
		}
	}
	//执行逆变换
	IFFT_1D(IFFT_in, out, Num_x_re); //执行IFFT
	//乘以padding_rate
	for (int i = 0; i < Num_x_re; i++)
	{
		out[i] = out[i] * Padding_rate;
	}
	fftw_free(FFT_out); //free
	fftw_free(IFFT_in); //free
}

void HOS_reconstruct(struct HOS_parameters *HOS_pa, struct Grid_position *Grid,
					 struct HOS_data_reconstruct *HOS_data_rc)
{
	double Phi_sum[Num_x]{};
	double temp_out[Num_x]{};  //U
	double temp_out1[Num_x]{}; //W

	for (int m = 0; m < Num_M; m++)
	{
		for (int i = 0; i < Num_x; i++)
		{
			Phi_sum[i] += HOS_pa->Phi_m[i][m];
		}
	}
	//FFT补零插值得到重构的Eta
	double temp_Eta_re[Num_x_re];
	FFT_interpolation(HOS_pa->Eta, HOS_data_rc->Eta);
	//FFT补零插值得到重构的Phi_s_x;phi_s_z
	double temp_Phi_s_x_re[Num_x_re];
	double temp_Phi_s_z_re[Num_x_re];
	FFT_interpolation(HOS_pa->Phi_s_dx, HOS_data_rc->Phi_s_x);
	FFT_interpolation(HOS_pa->Phi_s_z, HOS_data_rc->Phi_s_z);

	//FFT补零插值得到重构的速度场
	for (int zz = 0; zz < Num_z; zz++)
	{
		double z = Grid->z[zz]; //网格z向高度；
		//基函数系数，z>eta:U,W=0;0<z<=eta:exp=1+z+0.5*z^2;

		double kx[Num_x]; //kz*eigenfunction

		//执行fft
		for (int i = 0; i < Num_x; i++)
		{
			kx[i] = HOS_pa->Kx[i] * Eigen_function(abs(HOS_pa->Kx[i]), z);
		}
		Derivation_x(Phi_sum, temp_out, kx);	 //求解x方向速度
		Derivation_z(Phi_sum, temp_out1, kx, 1); //求解z方向速度
		/*
		FFT插值函数
		*/
		double temp_re[Num_x_re];  //暂存速度场U
		double temp_re1[Num_x_re]; //暂存速度场W
		FFT_interpolation(temp_out, temp_re);
		FFT_interpolation(temp_out1, temp_re1);

		for (int i = 0; i < Num_x_re; i++)
		{
			HOS_data_rc->U[i][zz] = temp_re[i];
			HOS_data_rc->W[i][zz] = temp_re1[i];
		}
	}
}

//五阶波浪函数
double Phi_fifth(double x, double z, double t)
{
	double a = 0.0;
	double b = 0.0;
	double c = 0.0;

	a = exp(k * z) * sin(k * x) * (E - 0.5 * pow(E, 3) - 34. / 27. * pow(E, 5));
	b = exp(2. * k * z) * sin(2. * k * x) * 0.5 * pow(E, 4);
	c = exp(3. * k * z) * sin(3. * k * x) * 1. / 12. * pow(E, 5);

	return (a + b + c) / (sqrt(pow(k, 3) / g));
}
double Eta_fifth(double x, double t)
{
	double a = 0.0;
	double b = 0.0;
	double c = 0.0;
	double d = 0.0;
	double e = 0.0;
	a = cos(k * x - w * t) * (E - 3.0 / 8.0 * pow(E, 3) - 422.0 / 384.0 * pow(E, 5));
	b = cos(2 * (k * x - w * t)) * (1.0 / 2.0 * pow(E, 2) + 1.0 / 3.0 * pow(E, 4));
	c = cos(3 * (k * x - w * t)) * (3.0 / 8.0 * pow(E, 3) + 297.0 / 384.0 * pow(E, 5));
	d = cos(4 * (k * x - w * t)) * (1.0 / 3.0 * pow(E, 4));
	e = cos(5 * (k * x - w * t)) * (125.0 / 384.0 * pow(E, 5));
	return (a + b + c + d + e) / k;
}

double mu_x1(double x)
{
	double temp = 0.0;
	double R = 0.0;
	double Rin = 0.0;
	double Rout = 0.0;
	double x_center = 0.0;
	x_center = (x_2 + x_3) / 2.0;
	Rin = abs(x_2 - x_center);
	Rout = abs(x_1 - x_center);
	R = abs(x - x_center);
	if (R <= Rin)
	{
		temp = 1;
	}
	else if (R >= Rout)
	{
		temp = 0;
	}
	else
	{
		temp = pow((cos((pi * (R - Rin)) / (2 * (Rout - Rin)))), 2);
	}
	return temp;
}

double mu_y1(double y)
{
	double temp = 0.0;
	double R = 0.0;
	double Rin = 0.0;
	double Rout = 0.0;
	double y_center = 0.0;
	y_center = (y_2 + y_3) / 2.0;
	Rin = abs(y_2 - y_center);
	Rout = abs(y_1 - y_center);
	R = abs(y - y_center);
	if (R <= Rin)
	{
		temp = 1;
	}
	else if (R >= Rout)
	{
		temp = 0;
	}
	else
	{
		temp = pow((cos((pi * (R - Rin)) / (2 * (Rout - Rin)))), 2);
	}
	return temp;
}
