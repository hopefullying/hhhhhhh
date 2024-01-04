
#include<stdio.h>
#include<windows.h>
#include <math.h>
#include <graphics.h>
#pragma once
using namespace std;
class test
{
public:

	//param
	FILE* DATA;
	FILE* PRINT;

	//从数据中读取的数据
	double δ;                //电极厚度
	int n;					  //电极总数	
	double* dz;				  //电极间距
	int* N;					  //相邻电极间要划分的步长数
	double* V;				  //电极电位
	double r1;				  //电极内孔径半径
	int M1;					  //r1范围内等步长划分网格数
	double r2;				  //电极内孔径到封闭边界处的径向距离
	int M2;					  //r2范围内等步长划分网格数
	double e;				  //迭代控制精度e
	int NST;				  //输出打印空间电位时网格点间隔数
	int INS;                  //轴上电位作等距插值时的步长数


	double dv;				  //要求扫描搜索等电位线的电位间隔值
	double* Ve;               //要求扫描搜索等电位线的电位值
	int m = 0;                //要求扫描搜索等电位线的电位值个数
	int m0;                   //要求扫描搜索等电位线的电位值个数最大值
	int row = 0, col = 0;     //row为行数，col为列数
	double** E;               //电场分布



	double** φd;             //格点电位残差
	double φds;              //格点电位残差和
	double φdm;              //格点电位最大残差
	double φda;              //格点电位平均残差
	double φda1;             //格点电位平均残差
	double φda2;             //格点电位平均残差
	double* zz;			      //每个小格子长度
	double* ZZ;				  //格点z坐标
	double* rr;		          //每个小格子宽度
	double* RR;				  //格点r坐标
	double c1, c2, c3, c4, c0;	 //五点差分公式中系数
	double ω;				  //SOR中超张弛迭代因子
	double λ;                //残差平均值之比λ
	double ωλ;              //加速因子ωλ
	double μλ;              //ωλ中的中间变量
	double ωm;               //修正的最佳迭代因子
	double ωm1;              //修正ω因子
	double φ0;				  //中心点电位
	int T;                    //迭代轮次
	int Time;                 //迭代次数
	int Time_total;           //迭代总次数
	int* dj;				//电极右边界格点位置（列数）
	double* ze;				//等电位点的z坐标
	double* re;             //等电位点的r坐标
	double* zeb;			//图中放大的等电位点坐标z坐标
	double* reb;            //图中放大的等电位点坐标r坐标



	// function 
	void Myset(FILE* data, FILE* out) {

		DATA = data;
		PRINT = out;

		int i, j;

		//判断是否打开成功
		if (data == NULL) {
			printf("CANNOT OPEN DATA FILE! \n");
			exit(0);

		}
		if (out == NULL) {
			printf("CANNOT OPEN PRINT FILE! \n");
			exit(0);

		}

		fscanf_s(data, "电极厚度δ=%lfmm;\n", &δ, sizeof(double));//从数据文件里读取数据
		fscanf_s(data, "电极总数n=%d;\n", &n, sizeof(int));
		fscanf_s(data, "相邻电极间距离dz=", sizeof(double));


		//dz = (double*)calloc(n * sizeof(double), sizeof(double
		dz = new double[n];

		for (i = 0; i < n; i++) {
			fscanf_s(data, "%lfmm;", &dz[i], sizeof(double));
		}
		fscanf_s(data, "\n相邻电极间要划分的步长数N=", sizeof(int));

		N = new int[n];
		//N = (int*)calloc(n * sizeof(int), sizeof(int));
		for (i = 0; i < n; i++) {
			fscanf_s(data, "%d;", &N[i], sizeof(int));
		}
		fscanf_s(data, "\n电极电位V=", sizeof(double));

		//V = (double*)calloc(n * sizeof(double), sizeof(double
		V = new	 double[n];

		for (i = 0; i < n; i++) {
			fscanf_s(data, "%lfV;", &V[i], sizeof(double));
		}
		fscanf_s(data, "\n电极内孔径半径r1=%lfmm;\n", &r1, sizeof(double));
		fscanf_s(data, "r1范围内等步长划分的网格数M1=%d;\n", &M1, sizeof(int));
		fscanf_s(data, "电极内孔径沿到封闭边界处的径向距离r2=%lfmm;\n", &r2, sizeof(double));
		fscanf_s(data, "r2范围内等步长划分的网格数M2=%d;\n", &M2, sizeof(int));
		fscanf_s(data, "迭代控制精度e=%lfV;\n", &e, sizeof(double));
		fscanf_s(data, "输出打印空间电位时网格点间隔数NST=%d;\n", &NST, sizeof(int));
		fscanf_s(data, "轴上电位作等距插值时的步长数INS=%d;\n", &INS, sizeof(int));
		fscanf_s(data, "等电位线的电位间隔或电位值=");
		//扫描输入多个电位值需手动向dat文件中填入
		//Ve = (double*)calloc(100 * sizeof(double), sizeof(double));//给要求扫描搜索等电位线的电位值分配空间
		Ve = new double[100];

		while (fscanf_s(data, "%lfV;", &Ve[m], sizeof(double)) != EOF)
			m++;//要求扫描的电位值个数
		if (m == 1)//若该项只有一个数则输入的是dv
		{
			dv = Ve[0];
			m0 = 100 / dv;
			for (int j = 0; j < m0; j++)
			{
				Ve[j] = (j + 1) * dv;
			}
		}


		// 赋值
		test::δ = δ;
		test::n = n;
		test::dz = dz;
		test::N = N;
		test::V = V;
		test::r1 = r1;
		test::M1 = M1;
		test::r2 = r2;
		test::M2 = M2;
		test::e = e;
		test::NST = NST;
		test::INS = INS;
		test::Ve = Ve;




		// print
		fprintf(out, "\t\t\t\t【原始数据】\n");
		fprintf(out, "电极厚度δ=%lfmm;\n", δ);
		fprintf(out, "电极总数n=%ld;\n", n);
		fprintf(out, "相邻电极间距离dz=");
		for (int i = 0; i < n; i++) {
			fprintf(out, "%lfmm;", dz[i]);
		}
		fprintf(out, "\n相邻电极间要划分的步长数N=");
		for (int i = 0; i < n; i++) {
			fprintf(out, "%ld;", N[i]);
		}
		fprintf(out, "\n电极电位V=");
		for (int i = 0; i < n; i++) {
			fprintf(out, "%lfV;", V[i]);
		}
		fprintf(out, "\n电极内孔径半径r1=%lfmm;\n", r1);
		fprintf(out, "r1范围内等步长划分网格数M1=%ld;\n", M1);
		fprintf(out, "电极内孔径沿到封闭边界处的径向距离r2=%lfmm;\n", r2);
		fprintf(out, "r2范围内等步长划分网格数M2=%ld;\n", M2);
		fprintf(out, "迭代控制精度e=%lfV;\n", e);
		fprintf(out, "输出打印空间电位时网格点间隔数NST=%ld;\n", NST);
		fprintf(out, "轴上电位作等距插值时的步长数INS=%ld;\n", INS);
		fprintf(out, "等电位线的电位间隔或电位值=");
		for (int i = 0; i < m; i++) {
			fprintf(out, "%lfV;", Ve[i]);
		}

	}
	void init(void) {
		int i, j;
		int a;
		int k = 0;
		for (i = 0; i < n; i++)  //纵向网格初始化
		{
			col = col + N[i] + 1;
		}
		row = M1 + M2 + 1;        //横向网格初始化
		printf("\n步长划分完毕。共 %d 行 ， %d 列；\n", row, col);
		fprintf(PRINT, "\n步长划分完毕。共 %d 行 ， %d 列；\n", row, col);
		E = (double**)calloc(row * sizeof(double*), sizeof(double));  //为格点的电场分配内存（二维）
		for (i = 0; i < row; i++)
			E[i] = (double*)calloc(col * sizeof(double), sizeof(double));

		for (i = 0; i < row; i++)//赋E初值全部为零
		{
			for (j = 0; j < col; j++)
			{
				E[i][j] = 0;
			}
		}


		//电极为等势体，电极两侧电压应相等且为已知电极电压
		for (i = 0; i <= M2; i++)
		{
			for (j = 0; j < n; j++)
			{
				E[i][k + N[j]] = V[j];
				E[i][k + N[j] + 1] = V[j];
				k = k + N[j] + 1;
			}
		}



		for (i = 0; i < row; i++)//荧光屏（整列）全部赋值V[n-1]
		{
			E[i][col - 1] = V[n - 1];
		}
		k = 0;

		//为电极之间电压赋值（插值）
		for (i = 0; i < M2; i++)//有电极的范围内电场赋值
		{
			k = 0;
			a = 0;
			for (j = k + 1; j <= k + N[a] - 1; j++)
			{
				E[i][j] = E[i][j - 1] + V[a] / N[a];   //相邻电极间线性插值	          
			}
			k = k + N[a] + 1;
			for (a = 1; a < n; a++)
			{
				for (j = k + 1; j <= k + N[a] - 1; j++)
				{
					E[i][j] = E[i][j - 1] + (V[a] - V[a - 1]) / N[a];   //相邻电极间线性插值	          
				}
				k = k + N[a] + 1;
			}
			for (j = 0; j <= col; j++)
			{
				E[i + 1][j] = E[i][j];
			}
			//第一组相邻电极间电场赋值结束
		}


		//测试是否成功赋值
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < col; j++)
			{
				printf(" %lf ", E[i][j]);
			}
			printf("\n");
		}

	};
	void SOR(void) {
		int i, j, k;
		double h1, h2, h3, h4, r0;

		φd = (double**)calloc(row * sizeof(double*), sizeof(double));  //为每个格子的电位残差分配内存（二维）
		for (i = 0; i < row; i++)
			φd[i] = (double*)calloc(col * sizeof(double), sizeof(double));

		for (i = 0; i < row; i++)
		{
			for (j = 0; j < col; j++)
			{
				int flag = 1;             //边界点判断参量：边界1，非边界0
				int K = 0;                //计数器
				if ((i != 0) && (j != 0) && (j != col - 1))        //上、左、右边界点判断
				{
					flag = 0;
					for (k = 0; k < n; k++)                       //电极位置边界点也不参与计算
					{
						if (i <= M2)
						{
							if (j >= (K + N[k]) && j <= (K + N[k] + 1))
							{
								flag = 1;
								break;
							}
						}
						K = K + N[k] + 1;
					}
				}
				if (flag == 0)//其余的非边界点参与计算
				{
					c1 = 2 / (zz[j - 1] * (zz[j - 1] + zz[j]));//SOR中c系数计算 c1=2/h1(h1+h2)
					c2 = 2 / (zz[j] * (zz[j - 1] + zz[j]));//c2=2/h2(h1+h2)

					if (i == (row - 1))//对于c3,c4需要分类讨论
					{
						c3 = 0;//最下面一行c3=0
						c4 = 4 / (rr[i - 1] * rr[i - 1]);//最下面一行c4=4/(h4*h4)
					}
					else
					{
						c3 = (2 * RR[i] - rr[i - 1]) / (RR[i] * rr[i] * (rr[i - 1] + rr[i]));//其余情况的c3
						c4 = (2 * RR[i] + rr[i]) / (RR[i] * rr[i - 1] * (rr[i - 1] + rr[i]));//其余情况的c4
					}
					c0 = c1 + c2 + c3 + c4;

					if (i == (row - 1))//对于电位φ0需要分类讨论
					{
						φ0 = (1 - ω) * E[i][j] + ω * (c1 * E[i][j - 1] + c2 * E[i][j + 1] + c4 * E[i - 1][j]) / c0;
					}
					else
					{
						φ0 = (1 - ω) * E[i][j] + ω * (c1 * E[i][j - 1] + c2 * E[i][j + 1] + c3 * E[i + 1][j] + c4 * E[i - 1][j]) / c0;
					}
					φd[i][j] = fabs(φ0 - E[i][j]);//计算每个点的电位残差（求绝对值）
					E[i][j] = φ0;
				}
			}
		}
		φds = 0;    //初始化残差和
		φdm = 0;    //初始化最大残差
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < col; j++)
			{
				φds = φds + φd[i][j];
				if (φd[i][j] > φdm)
					φdm = φd[i][j];
			}
		}
		φda = φds / (row * col);     //计算平均残差	
	};
	void calculate() {

		//计算网格尺寸
		int i, j, k;
		//dj = (int*)calloc((n + 1) * sizeof(int), sizeof(int));   //为电极右侧点位置（对应的列数）分配内存
		dj = new int[n + 1];

		dj[0] = 0;                                 //令dj[0]=0
		for (i = 1; i < n + 1; i++)//为dj的每个元素赋值
		{
			if (i == n)
				dj[i] = dj[i - 1] + N[i - 1];//荧光屏处（相当于第n个电极右边界）格点位置
			else
				dj[i] = dj[i - 1] + N[i - 1] + 1;//一般电极右边界格点位置
		}
		//zz = (double*)calloc((col - 1) * sizeof(double), sizeof(double));//为每个格子的（横向）长度分配内存
		zz = new double[col - 1];

		k = 0;
		for (i = 0; i < n; i++)                                          //计算每个格子的（横向）长度
		{
			for (j = k; j < (k + N[i]); j++)                             //计算电极之间的格子的（横向）长度
			{
				zz[j] = dz[i] / N[i];
			}
			k = k + N[i];
			if (k >= col - 1)                                             //计算到达屏上即跳出
				break;
			zz[j] = δ;
			k = k + 1;
		}

		//ZZ = (double*)calloc(col * sizeof(double), sizeof(double));//为每个位置的（横向）横坐标分配内存
		ZZ = new double[col];


		ZZ[0] = 0;
		for (i = 1; i < col; i++)
		{
			ZZ[i] = ZZ[i - 1] + zz[i - 1];
			fprintf(PRINT, "%lf\t\t\t", ZZ[i]);
		}
		fprintf(PRINT, "\n");

		//rr = (double*)calloc((row - 1) * sizeof(double), sizeof(double));//为每个格子的（纵向）宽度分配内存
		rr = new double[row - 1];

		for (i = 0; i < row - 1; i++)                                 //计算每个格子的（纵向）宽度
		{
			if (i < M2)
				rr[i] = r2 / M2;
			else
				rr[i] = r1 / M1;

		}

		//RR = (double*)calloc(row * sizeof(double), sizeof(double));//为每个位置的（纵向）纵坐标分配内存	
		RR = new double[row];
		RR[0] = 0;
		for (i = 1; i < row; i++)
		{
			RR[i] = RR[i - 1] + rr[i];
			fprintf(PRINT, "%lf\t\t\t", RR[i]);
		}


		//计算电场分布
		//int i, j, k;             //计数使用字母
		int sign;                //迭代标志
		ω = 1;                  //超张弛迭代因子
		T = 1;                   //迭代轮次
		Time = 1;                //迭代次数为1
		Time_total = 0;          //总迭代次数为0

		//第一轮迭代
		printf("\n第%d轮迭代\n", T);
		fprintf(PRINT, "\n第%d轮迭代\n", T);
		SOR();           //开始第一轮迭代
		Time_total = Time_total + 1;    //总迭代次数+1
		printf("\t\t迭代次数：%d\t\t迭代因子值:%lf\t\t", Time_total, ω);
		printf("平均残差：%lf\t\t最大残差：%.7lf\n", φda, φdm);
		fprintf(PRINT, "\t\t迭代次数：%d\t\t迭代因子值:%lf\t\t", Time_total, ω);
		fprintf(PRINT, "平均残差：%lf\t\t最大残差：%.7lf\n", φda, φdm);

		//第二轮迭代
		ω = 1.375;       //输入ω为“卡瑞建议值”
		T = T + 1;
		printf("\n第%d轮迭代\n", T);
		fprintf(PRINT, "\n第%d轮迭代\n", T);
		Time_total = 0;
		for (Time = 1; Time < 13; Time++)
		{
			SOR();                            //12次迭代
			Time_total = Time_total + 1;
			if (Time == 11)
				φda1 = φda;
			else if (Time == 12)
				φda2 = φda;

		}
		printf("\t\t迭代次数：%d\t\t迭代因子值:%lf\t\t", Time_total, ω);
		printf("平均残差：%lf\t\t最大残差：%.7lf\n", φda, φdm);
		fprintf(PRINT, "\t\t迭代次数：%d\t\t迭代因子值:%lf\t\t", Time_total, ω);
		fprintf(PRINT, "平均残差：%lf\t\t最大残差：%.7lf\n", φda, φdm);
		λ = φda2 / φda1;                        //求残差平均值之比λ
		μλ = (λ + ω - 1) / (sqrt(λ) * ω);      //求中间变量μλ
		ωλ = 2 / (1 + sqrt(1 - μλ * μλ));    //求加速因子ωλ
		ωm = 1.25 * ωλ - 0.5;                     //求得最佳迭代因子

		//第三轮迭代
		do
		{
			ω = ωm;                             //输入ω为上一轮求得的最佳迭代因子
			T = T + 1;
			printf("\n第%d轮迭代\n", T);
			fprintf(PRINT, "\n第%d轮迭代\n", T);
			Time_total = 0;
			for (Time = 1; Time < 13; Time++)
			{
				SOR();                            //12次迭代
				Time_total = Time_total + 1;
				if (Time == 11)
					φda1 = φda;
				else if (Time == 12)
					φda2 = φda;

			}
			printf("\t\t迭代次数：%d\t\t迭代因子值:%lf\t\t", Time_total, ω);
			printf("平均残差：%lf\t\t最大残差：%.7lf\n", φda, φdm);
			fprintf(PRINT, "\t\t迭代次数：%d\t\t迭代因子值:%lf\t\t", Time_total, ω);
			fprintf(PRINT, "平均残差：%lf\t\t最大残差：%.7lf\n", φda, φdm);
			ωm1 = ωm;
			λ = φda2 / φda1;
			μλ = (λ + ω - 1) / (sqrt(λ) * ω);
			ωλ = 2 / (1 + sqrt(1 - μλ * μλ));
			ωm = 1.25 * ωλ - 0.5;
		} while (fabs((ωm - ωm1) / (2 - ωm1)) >= 0.05);               //若连续两轮求得ωm1和ωm满足判定式，则固定ωm1
		ω = ωm;

		do
		{
			T = T + 1;
			printf("第%d轮迭代", T);
			fprintf(PRINT, "\n第%d轮迭代\n", T);
			SOR();//最佳迭代因子确定后继续迭代
			Time_total = Time_total + 1;
			printf("\t\t迭代次数：%d\t\t迭代因子值:%lf\t\t", Time_total, ω);
			printf("平均残差：%lf\t\t最大残差：%.7lf\n", φda, φdm);
			fprintf(PRINT, "\t\t迭代次数：%d\t\t迭代因子值:%lf\t\t", Time_total, ω);
			fprintf(PRINT, "平均残差：%lf\t\t最大残差：%.7lf\n", φda, φdm);
			sign = 0;
			for (i = 0; i < row; i++)
				for (j = 0; j < col; j++)
				{
					if (φd[i][j] >= e)//迭代达到控制精度ε则停止
						sign = sign + 1;
				}
		} while (sign > 0);
		//网格点电位计算完毕

		//输出计算得到的空间电场分布
		printf("\n\t\t\t\t【网格点上的电位值】\n");
		fprintf(PRINT, "\n\t\t\t\t【网格点上的电位值】\n");
		for (i = 0; i <= col; i++)
		{
			printf("%d\t\t", i);
			fprintf(PRINT, "%d\t\t", i);
			i = i + NST - 1;
		}
		printf("\n");
		fprintf(PRINT, "\n");
		for (i = 0; i < row; i++)
		{
			printf("%d\t", i);
			fprintf(PRINT, "%d\t", i);
			for (j = 0; j <= col; j++)
			{
				fprintf(PRINT, "%lf\t", E[i][j]);
				printf("%lf\t", E[i][j]);
				j = j + NST - 1;
			}
			printf("\n");
			fprintf(PRINT, "\n");
			i = i + NST - 1;
		}


	};
	void painting() {
		int i, i1, j, j1, k, l, c, b, pan, f, f1, f2, f3, f4, f5, f6, * hi;
		double a = 0, zlm, rlm, zrlm;
		if (m == 1)//如果输入的是间隔值（只有一个数）则计算等电位值
		{
			for (m = 1; a < (V[n - 1] - dv); m++)//计算电位数量
				a = a + dv;
			Ve = (double*)calloc((m - 1) * sizeof(double), sizeof(double));
			for (i = 0; Ve[i] < (V[n - 1] - dv); i++)
				Ve[i + 1] = Ve[i] + dv;//计算各等位线的电位值
		}

		//以下定义图形并计算等位线坐标
		initgraph(100 + int(10 * ZZ[col - 1]), 100 + int(10 * RR[row - 1]));//画面距窗口边缘四周均50，坐标值均扩大10倍
		setbkcolor(BLACK);//设定黑色背景
		cleardevice();
		rectangle(50, 50, 50 + int(10 * ZZ[col - 1]), 50 + int(10 * RR[row - 1]));//距窗口边缘四周均50画矩形

		dj = (int*)calloc((n + 1) * sizeof(int), sizeof(int));//为电极右侧点位置分配内存
		dj[0] = 0;//令dj[0]=0
		for (i = 1; i < n + 1; i++)//为dj的每个元素赋值
		{
			if (i == n)
				dj[i] = dj[i - 1] + N[i - 1];//荧光屏处（相当于第n个电极右边界）格点位置
			else
				dj[i] = dj[i - 1] + N[i - 1] + 1;//一般电极右边界格点位置
		}

		for (i = 1; i < n; i++)
			bar(50 + int(10 * ZZ[dj[i] - 1]), 50, 50 + int(10 * (ZZ[dj[i]])), 50 + int(10 * r2));   //画电极
		ze = (double*)calloc(row * col * sizeof(double), sizeof(double));   //等电位点z坐标值分配内存
		re = (double*)calloc(row * col * sizeof(double), sizeof(double));   //等电位点r坐标值分配内存
		zeb = (double*)calloc(row * col * sizeof(double), sizeof(double));  //扩大10倍的等电位点z坐标值分配内存
		reb = (double*)calloc(row * col * sizeof(double), sizeof(double));  //扩大10倍的等电位点r坐标值分配内存
		zlm = zz[0];
		for (i = 0; i < col - 1; i++)
		{
			if (zlm < zz[i])
				zlm = zz[i];        //计算z方向最大格点间距
		}
		rlm = rr[0];
		if (rlm < rr[row - 2])
			rlm = rr[row - 2];        //计算r方向最大格点间距
		zrlm = sqrt(zlm * zlm + rlm * rlm);	 //计算对角线方向最大间距
		for (k = 0; k < m; k++)
		{
			f1 = 0;    //记录顶行等电位点数量初始化
			l = 0;     //等电位点数量初始化
			for (i1 = 0; i1 < row; i1++)
			{
				for (j1 = 0; j1 < col - 1; j1++)   //z轴方向线性插值计算等电位点坐标
				{
					if (((E[i1][j1] < Ve[k]) & (Ve[k] < E[i1][j1 + 1])) || ((E[i1][j1] > Ve[k]) & (Ve[k] > E[i1][j1 + 1])) || (E[i1][j1] == Ve[k]))
					{
						re[l] = RR[i1];
						if (re[l] == RR[0])
							f1++;          //顶行标记值记录顶行等电位点数量用于控制等位线绘制
						if (E[i1][j1] == Ve[k])
							ze[l] = ZZ[j1];
						else
							ze[l] = ZZ[j1] + (Ve[k] - E[i1][j1]) * zz[j1] / (E[i1][j1 + 1] - E[i1][j1]);
						reb[l] = 50 + 10 * re[l];             //扩大10倍后的图形坐标 
						zeb[l] = 50 + 10 * ze[l];             //扩大10倍后的图形坐标 
						l++;
					}
				}

				for (j = 0; j < col; j++)        //r轴方向线性插值计算等电位点坐标
				{
					if (i1 < row - 1)
					{
						i = i1;
						if (((E[i][j] < Ve[k]) & (Ve[k] < E[i + 1][j])) || ((E[i][j] > Ve[k]) & (Ve[k] > E[i + 1][j])))
						{
							ze[l] = ZZ[j];
							re[l] = RR[i] + (Ve[k] - E[i][j]) * rr[i] / (E[i + 1][j] - E[i][j]);
							reb[l] = 50 + 10 * re[l];              //扩大10倍后的图形坐标
							zeb[l] = 50 + 10 * ze[l];              //扩大10倍后的图形坐标
							l++;
						}
					}
				}


			}


			//等位线绘制
			f = 0;   //连接线段起点位置变量初始化
			c = 0;   //储存连接过的等电位点位置循环变量初始化
			f6 = 0;  //用于记录电位等于电极电位时的电极位置初始化
			hi = (int*)calloc(l * sizeof(int), sizeof(int));   //储存连接过的等电位点位置
			do
			{

				fprintf(PRINT, "\n电位值为%lf的等位线的各点坐标值：\n", Ve[k]);
				setcolor(20 * k * RED + 10 * k * GREEN + BLUE);    //确定等位线颜色
				for (i = 0; i < l; i++)
				{
					pan = 1;
					for (b = 0; b < c; b++)    //判断是否为已连接过的点
					{
						if (i == hi[b])
						{
							pan = 0;
							break;
						}
					}
					if ((pan == 1) & (re[i] == RR[0]))    //用于有鞍点时找出下一条等位线的起点位置
					{
						f = i;
						break;
					}
				}
				for (i = f; i < l; i++)
				{
					for (j = 0; j < l; j++)
					{
						pan = 1;
						for (b = 0; b < c; b++)    //判断是否为已连接过的点
						{
							if (j == hi[b])
							{
								pan = 0;
								break;
							}
						}
						if ((j != i) & (pan == 1) & (sqrt((zeb[j] - zeb[i]) * (zeb[j] - zeb[i]) + (reb[j] - reb[i]) * (reb[j] - reb[i])) < 12 * zrlm) & (fabs(zeb[j] - zeb[i]) < 12 * zlm) & (fabs(reb[j] - reb[i]) < 12 * rlm))
						{    //如果第二点不是第一点 且 不是已连接过的点 且 两点之间距离满足要求
							f3 = 0;
							for (f2 = 0; f2 < n; f2++)    //如果等于电极电位则记录电极位置并做标记
							{
								if (Ve[k] == V[f2])
								{
									f3 = 1;    //做标记
									f6 = f2;   //记录电极位置
									break;
								}
							}
							if ((f3 == 1) & ((ZZ[dj[f6 + 1]] == ze[j]) || (ZZ[dj[f6 + 1] - 1] == ze[j])) & (RR[0] == re[j]))
							{    //如果是电极电位且第二点为电极顶部所在位置则从上往下“Z”字型连接电极电位点
								for (f2 = 0; f2 < M2; f2++)
								{
									fprintf(PRINT, "(%lf\t,\t%lf\t)\n", re[i], ze[i]);   //输出第一点坐标
									if (re[i] == RR[0])
										f1--;     //如果第一点在顶行则标记值减1
									line(int(zeb[i]), int(reb[i]), int(zeb[j]), int(reb[j]));
									hi[c] = i;      //记录第一点位置
									c++;
									fprintf(PRINT, "(%lf\t,\t%lf\t)\n", re[j], ze[j]);   //输出第二点坐标
									if (re[j] == RR[0])
										f1--;	  //如果第二点在顶行则标记值减1			
									for (f4 = 0; f4 < l; f4++)
									{
										if ((ze[f4] == ze[i]) & (re[f4] == re[i] + r2 / M2))
										{
											i = f4;      //寻找下一个第一点
											break;
										}
									}
									line(int(zeb[i]), int(reb[i]), int(zeb[j]), int(reb[j]));
									hi[c] = j;
									c++;
									for (f5 = 0; f5 < l; f5++)
									{
										if ((ze[f5] == ze[j]) & (re[f5] == re[j] + r2 / M2))
										{
											j = f5;      //寻找下一个第二点
											break;
										}
									}
								}
								fprintf(PRINT, "(%lf\t,\t%lf\t)\n", re[i], ze[i]);  //输出电极底端左侧电位点坐标
								if (re[i] == RR[0])
									f1--;			//如果此点在顶行则标记值减1		
								line(int(zeb[i]), int(reb[i]), int(zeb[j]), int(reb[j]));   //连接电极底端水平线
								hi[c] = i;            //记录电极底端左侧电位点位置
								c++;
								break;
							}
							else if ((f3 == 1) & ((zz[dj[f6 + 1]] == ze[j]) || (ZZ[dj[f6 + 1] - 1] == ze[j])) & (RR[0] != re[j]))
							{   //如果是电极电位且第二点不是电极顶部所在位置则从下往上倒“Z”字型连接电极电位点
								fprintf(PRINT, "(%lf\t,\t%lf\t)\n", re[i], ze[i]);
								if (re[i] == RR[0])
									f1--;
								line(int(zeb[i]), int(reb[i]), int(zeb[j]), int(reb[j]));
								hi[c] = i;
								c++;
								fprintf(PRINT, "(%lf\t,\t%lf\t)\n", re[j], ze[j]);
								if (re[j] == RR[0])
									f1--;
								for (f4 = 0; f4 < l; f4++)
								{
									if ((ze[f4] == ze[dj[f6 + 1]]) & (re[f4] == re[M2]))
									{
										i = f4;
										j = f4 - 1;
										break;
									}
								}
								for (f2 = 0; f2 < M2; f2++)
								{
									fprintf(PRINT, "(%lf\t,\t%lf\t)\n", re[i], ze[i]);
									if (re[i] == RR[0])
										f1--;
									line(int(zeb[i]), int(reb[i]), int(zeb[j]), int(reb[j]));
									hi[c] = i;
									c++;
									fprintf(PRINT, "(%lf\t,\t%lf\t)\n", re[j], ze[j]);
									if (re[j] == RR[0])
										f1--;
									for (f4 = 0; f4 < l; f4++)
									{
										if ((ze[f4] == ze[i]) & (re[f4] == re[i] - r2 / M2))
										{
											i = f4;
											break;
										}
									}
									line(int(zeb[j]), int(reb[j]), int(zeb[i]), int(reb[i]));
									hi[c] = j;
									c++;
									for (f5 = 0; f5 < l; f5++)
									{
										if ((ze[f5] == ze[j]) & (re[f5] == re[j] - r2 / M2))
										{
											j = f5;
											break;
										}
									}
								}
								fprintf(PRINT, "(%lf\t,\t%lf\t)\n", re[i], ze[i]);
								if (re[i] == RR[0])
									f1--;
								line(int(zeb[i]), int(reb[i]), int(zeb[j]), int(reb[j]));
								hi[c] = i;
								c++;
								break;
							}
							else
							{   //如果第二点不是电极所在位置则输出第一点坐标并连接第一二点
								fprintf(PRINT, "(%lf\t,\t%lf\t)\n", re[i], ze[i]);
								if (re[i] == RR[0])
									f1--;
								line(int(zeb[i]), int(reb[i]), int(zeb[j]), int(reb[j]));
								hi[c] = i;
								c++;
								break;
							}
						}
					}
					if ((re[j] == RR[0]) || (re[j] == RR[row - 1]))
					{   //如果第二点位于顶行或底行则输出第二点坐标并跳出循环判断是否有下一条等位线
						fprintf(PRINT, "(%lf\t,\t%lf\t)\n", re[j], ze[j]);
						if (re[j] == RR[0])
							f1--;
						hi[c] = j;
						c++;
						break;
					}
					else   //如果第二点不位于顶行或底行则将第二点作为下一次循环的第一点
						i = j - 1;
				}
			} while (f1 > 0);      //如果顶行标记值不为0则继续搜索连接下一条等位线
		}




	};



};