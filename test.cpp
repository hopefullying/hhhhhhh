
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

	//�������ж�ȡ������
	double ��;                //�缫���
	int n;					  //�缫����	
	double* dz;				  //�缫���
	int* N;					  //���ڵ缫��Ҫ���ֵĲ�����
	double* V;				  //�缫��λ
	double r1;				  //�缫�ڿ׾��뾶
	int M1;					  //r1��Χ�ڵȲ�������������
	double r2;				  //�缫�ڿ׾�����ձ߽紦�ľ������
	int M2;					  //r2��Χ�ڵȲ�������������
	double e;				  //�������ƾ���e
	int NST;				  //�����ӡ�ռ��λʱ���������
	int INS;                  //���ϵ�λ���Ⱦ��ֵʱ�Ĳ�����


	double dv;				  //Ҫ��ɨ�������ȵ�λ�ߵĵ�λ���ֵ
	double* Ve;               //Ҫ��ɨ�������ȵ�λ�ߵĵ�λֵ
	int m = 0;                //Ҫ��ɨ�������ȵ�λ�ߵĵ�λֵ����
	int m0;                   //Ҫ��ɨ�������ȵ�λ�ߵĵ�λֵ�������ֵ
	int row = 0, col = 0;     //rowΪ������colΪ����
	double** E;               //�糡�ֲ�



	double** ��d;             //����λ�в�
	double ��ds;              //����λ�в��
	double ��dm;              //����λ���в�
	double ��da;              //����λƽ���в�
	double ��da1;             //����λƽ���в�
	double ��da2;             //����λƽ���в�
	double* zz;			      //ÿ��С���ӳ���
	double* ZZ;				  //���z����
	double* rr;		          //ÿ��С���ӿ��
	double* RR;				  //���r����
	double c1, c2, c3, c4, c0;	 //����ֹ�ʽ��ϵ��
	double ��;				  //SOR�г��ųڵ�������
	double ��;                //�в�ƽ��ֵ֮�Ȧ�
	double �ئ�;              //�������Ӧئ�
	double �̦�;              //�ئ��е��м����
	double ��m;               //��������ѵ�������
	double ��m1;              //����������
	double ��0;				  //���ĵ��λ
	int T;                    //�����ִ�
	int Time;                 //��������
	int Time_total;           //�����ܴ���
	int* dj;				//�缫�ұ߽���λ�ã�������
	double* ze;				//�ȵ�λ���z����
	double* re;             //�ȵ�λ���r����
	double* zeb;			//ͼ�зŴ�ĵȵ�λ������z����
	double* reb;            //ͼ�зŴ�ĵȵ�λ������r����



	// function 
	void Myset(FILE* data, FILE* out) {

		DATA = data;
		PRINT = out;

		int i, j;

		//�ж��Ƿ�򿪳ɹ�
		if (data == NULL) {
			printf("CANNOT OPEN DATA FILE! \n");
			exit(0);

		}
		if (out == NULL) {
			printf("CANNOT OPEN PRINT FILE! \n");
			exit(0);

		}

		fscanf_s(data, "�缫��Ȧ�=%lfmm;\n", &��, sizeof(double));//�������ļ����ȡ����
		fscanf_s(data, "�缫����n=%d;\n", &n, sizeof(int));
		fscanf_s(data, "���ڵ缫�����dz=", sizeof(double));


		//dz = (double*)calloc(n * sizeof(double), sizeof(double
		dz = new double[n];

		for (i = 0; i < n; i++) {
			fscanf_s(data, "%lfmm;", &dz[i], sizeof(double));
		}
		fscanf_s(data, "\n���ڵ缫��Ҫ���ֵĲ�����N=", sizeof(int));

		N = new int[n];
		//N = (int*)calloc(n * sizeof(int), sizeof(int));
		for (i = 0; i < n; i++) {
			fscanf_s(data, "%d;", &N[i], sizeof(int));
		}
		fscanf_s(data, "\n�缫��λV=", sizeof(double));

		//V = (double*)calloc(n * sizeof(double), sizeof(double
		V = new	 double[n];

		for (i = 0; i < n; i++) {
			fscanf_s(data, "%lfV;", &V[i], sizeof(double));
		}
		fscanf_s(data, "\n�缫�ڿ׾��뾶r1=%lfmm;\n", &r1, sizeof(double));
		fscanf_s(data, "r1��Χ�ڵȲ������ֵ�������M1=%d;\n", &M1, sizeof(int));
		fscanf_s(data, "�缫�ڿ׾��ص���ձ߽紦�ľ������r2=%lfmm;\n", &r2, sizeof(double));
		fscanf_s(data, "r2��Χ�ڵȲ������ֵ�������M2=%d;\n", &M2, sizeof(int));
		fscanf_s(data, "�������ƾ���e=%lfV;\n", &e, sizeof(double));
		fscanf_s(data, "�����ӡ�ռ��λʱ���������NST=%d;\n", &NST, sizeof(int));
		fscanf_s(data, "���ϵ�λ���Ⱦ��ֵʱ�Ĳ�����INS=%d;\n", &INS, sizeof(int));
		fscanf_s(data, "�ȵ�λ�ߵĵ�λ������λֵ=");
		//ɨ����������λֵ���ֶ���dat�ļ�������
		//Ve = (double*)calloc(100 * sizeof(double), sizeof(double));//��Ҫ��ɨ�������ȵ�λ�ߵĵ�λֵ����ռ�
		Ve = new double[100];

		while (fscanf_s(data, "%lfV;", &Ve[m], sizeof(double)) != EOF)
			m++;//Ҫ��ɨ��ĵ�λֵ����
		if (m == 1)//������ֻ��һ�������������dv
		{
			dv = Ve[0];
			m0 = 100 / dv;
			for (int j = 0; j < m0; j++)
			{
				Ve[j] = (j + 1) * dv;
			}
		}


		// ��ֵ
		test::�� = ��;
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
		fprintf(out, "\t\t\t\t��ԭʼ���ݡ�\n");
		fprintf(out, "�缫��Ȧ�=%lfmm;\n", ��);
		fprintf(out, "�缫����n=%ld;\n", n);
		fprintf(out, "���ڵ缫�����dz=");
		for (int i = 0; i < n; i++) {
			fprintf(out, "%lfmm;", dz[i]);
		}
		fprintf(out, "\n���ڵ缫��Ҫ���ֵĲ�����N=");
		for (int i = 0; i < n; i++) {
			fprintf(out, "%ld;", N[i]);
		}
		fprintf(out, "\n�缫��λV=");
		for (int i = 0; i < n; i++) {
			fprintf(out, "%lfV;", V[i]);
		}
		fprintf(out, "\n�缫�ڿ׾��뾶r1=%lfmm;\n", r1);
		fprintf(out, "r1��Χ�ڵȲ�������������M1=%ld;\n", M1);
		fprintf(out, "�缫�ڿ׾��ص���ձ߽紦�ľ������r2=%lfmm;\n", r2);
		fprintf(out, "r2��Χ�ڵȲ�������������M2=%ld;\n", M2);
		fprintf(out, "�������ƾ���e=%lfV;\n", e);
		fprintf(out, "�����ӡ�ռ��λʱ���������NST=%ld;\n", NST);
		fprintf(out, "���ϵ�λ���Ⱦ��ֵʱ�Ĳ�����INS=%ld;\n", INS);
		fprintf(out, "�ȵ�λ�ߵĵ�λ������λֵ=");
		for (int i = 0; i < m; i++) {
			fprintf(out, "%lfV;", Ve[i]);
		}

	}
	void init(void) {
		int i, j;
		int a;
		int k = 0;
		for (i = 0; i < n; i++)  //���������ʼ��
		{
			col = col + N[i] + 1;
		}
		row = M1 + M2 + 1;        //���������ʼ��
		printf("\n����������ϡ��� %d �� �� %d �У�\n", row, col);
		fprintf(PRINT, "\n����������ϡ��� %d �� �� %d �У�\n", row, col);
		E = (double**)calloc(row * sizeof(double*), sizeof(double));  //Ϊ���ĵ糡�����ڴ棨��ά��
		for (i = 0; i < row; i++)
			E[i] = (double*)calloc(col * sizeof(double), sizeof(double));

		for (i = 0; i < row; i++)//��E��ֵȫ��Ϊ��
		{
			for (j = 0; j < col; j++)
			{
				E[i][j] = 0;
			}
		}


		//�缫Ϊ�����壬�缫�����ѹӦ�����Ϊ��֪�缫��ѹ
		for (i = 0; i <= M2; i++)
		{
			for (j = 0; j < n; j++)
			{
				E[i][k + N[j]] = V[j];
				E[i][k + N[j] + 1] = V[j];
				k = k + N[j] + 1;
			}
		}



		for (i = 0; i < row; i++)//ӫ���������У�ȫ����ֵV[n-1]
		{
			E[i][col - 1] = V[n - 1];
		}
		k = 0;

		//Ϊ�缫֮���ѹ��ֵ����ֵ��
		for (i = 0; i < M2; i++)//�е缫�ķ�Χ�ڵ糡��ֵ
		{
			k = 0;
			a = 0;
			for (j = k + 1; j <= k + N[a] - 1; j++)
			{
				E[i][j] = E[i][j - 1] + V[a] / N[a];   //���ڵ缫�����Բ�ֵ	          
			}
			k = k + N[a] + 1;
			for (a = 1; a < n; a++)
			{
				for (j = k + 1; j <= k + N[a] - 1; j++)
				{
					E[i][j] = E[i][j - 1] + (V[a] - V[a - 1]) / N[a];   //���ڵ缫�����Բ�ֵ	          
				}
				k = k + N[a] + 1;
			}
			for (j = 0; j <= col; j++)
			{
				E[i + 1][j] = E[i][j];
			}
			//��һ�����ڵ缫��糡��ֵ����
		}


		//�����Ƿ�ɹ���ֵ
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

		��d = (double**)calloc(row * sizeof(double*), sizeof(double));  //Ϊÿ�����ӵĵ�λ�в�����ڴ棨��ά��
		for (i = 0; i < row; i++)
			��d[i] = (double*)calloc(col * sizeof(double), sizeof(double));

		for (i = 0; i < row; i++)
		{
			for (j = 0; j < col; j++)
			{
				int flag = 1;             //�߽���жϲ������߽�1���Ǳ߽�0
				int K = 0;                //������
				if ((i != 0) && (j != 0) && (j != col - 1))        //�ϡ����ұ߽���ж�
				{
					flag = 0;
					for (k = 0; k < n; k++)                       //�缫λ�ñ߽��Ҳ���������
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
				if (flag == 0)//����ķǱ߽��������
				{
					c1 = 2 / (zz[j - 1] * (zz[j - 1] + zz[j]));//SOR��cϵ������ c1=2/h1(h1+h2)
					c2 = 2 / (zz[j] * (zz[j - 1] + zz[j]));//c2=2/h2(h1+h2)

					if (i == (row - 1))//����c3,c4��Ҫ��������
					{
						c3 = 0;//������һ��c3=0
						c4 = 4 / (rr[i - 1] * rr[i - 1]);//������һ��c4=4/(h4*h4)
					}
					else
					{
						c3 = (2 * RR[i] - rr[i - 1]) / (RR[i] * rr[i] * (rr[i - 1] + rr[i]));//���������c3
						c4 = (2 * RR[i] + rr[i]) / (RR[i] * rr[i - 1] * (rr[i - 1] + rr[i]));//���������c4
					}
					c0 = c1 + c2 + c3 + c4;

					if (i == (row - 1))//���ڵ�λ��0��Ҫ��������
					{
						��0 = (1 - ��) * E[i][j] + �� * (c1 * E[i][j - 1] + c2 * E[i][j + 1] + c4 * E[i - 1][j]) / c0;
					}
					else
					{
						��0 = (1 - ��) * E[i][j] + �� * (c1 * E[i][j - 1] + c2 * E[i][j + 1] + c3 * E[i + 1][j] + c4 * E[i - 1][j]) / c0;
					}
					��d[i][j] = fabs(��0 - E[i][j]);//����ÿ����ĵ�λ�в�����ֵ��
					E[i][j] = ��0;
				}
			}
		}
		��ds = 0;    //��ʼ���в��
		��dm = 0;    //��ʼ�����в�
		for (i = 0; i < row; i++)
		{
			for (j = 0; j < col; j++)
			{
				��ds = ��ds + ��d[i][j];
				if (��d[i][j] > ��dm)
					��dm = ��d[i][j];
			}
		}
		��da = ��ds / (row * col);     //����ƽ���в�	
	};
	void calculate() {

		//��������ߴ�
		int i, j, k;
		//dj = (int*)calloc((n + 1) * sizeof(int), sizeof(int));   //Ϊ�缫�Ҳ��λ�ã���Ӧ�������������ڴ�
		dj = new int[n + 1];

		dj[0] = 0;                                 //��dj[0]=0
		for (i = 1; i < n + 1; i++)//Ϊdj��ÿ��Ԫ�ظ�ֵ
		{
			if (i == n)
				dj[i] = dj[i - 1] + N[i - 1];//ӫ���������൱�ڵ�n���缫�ұ߽磩���λ��
			else
				dj[i] = dj[i - 1] + N[i - 1] + 1;//һ��缫�ұ߽���λ��
		}
		//zz = (double*)calloc((col - 1) * sizeof(double), sizeof(double));//Ϊÿ�����ӵģ����򣩳��ȷ����ڴ�
		zz = new double[col - 1];

		k = 0;
		for (i = 0; i < n; i++)                                          //����ÿ�����ӵģ����򣩳���
		{
			for (j = k; j < (k + N[i]); j++)                             //����缫֮��ĸ��ӵģ����򣩳���
			{
				zz[j] = dz[i] / N[i];
			}
			k = k + N[i];
			if (k >= col - 1)                                             //���㵽�����ϼ�����
				break;
			zz[j] = ��;
			k = k + 1;
		}

		//ZZ = (double*)calloc(col * sizeof(double), sizeof(double));//Ϊÿ��λ�õģ����򣩺���������ڴ�
		ZZ = new double[col];


		ZZ[0] = 0;
		for (i = 1; i < col; i++)
		{
			ZZ[i] = ZZ[i - 1] + zz[i - 1];
			fprintf(PRINT, "%lf\t\t\t", ZZ[i]);
		}
		fprintf(PRINT, "\n");

		//rr = (double*)calloc((row - 1) * sizeof(double), sizeof(double));//Ϊÿ�����ӵģ����򣩿�ȷ����ڴ�
		rr = new double[row - 1];

		for (i = 0; i < row - 1; i++)                                 //����ÿ�����ӵģ����򣩿��
		{
			if (i < M2)
				rr[i] = r2 / M2;
			else
				rr[i] = r1 / M1;

		}

		//RR = (double*)calloc(row * sizeof(double), sizeof(double));//Ϊÿ��λ�õģ���������������ڴ�	
		RR = new double[row];
		RR[0] = 0;
		for (i = 1; i < row; i++)
		{
			RR[i] = RR[i - 1] + rr[i];
			fprintf(PRINT, "%lf\t\t\t", RR[i]);
		}


		//����糡�ֲ�
		//int i, j, k;             //����ʹ����ĸ
		int sign;                //������־
		�� = 1;                  //���ųڵ�������
		T = 1;                   //�����ִ�
		Time = 1;                //��������Ϊ1
		Time_total = 0;          //�ܵ�������Ϊ0

		//��һ�ֵ���
		printf("\n��%d�ֵ���\n", T);
		fprintf(PRINT, "\n��%d�ֵ���\n", T);
		SOR();           //��ʼ��һ�ֵ���
		Time_total = Time_total + 1;    //�ܵ�������+1
		printf("\t\t����������%d\t\t��������ֵ:%lf\t\t", Time_total, ��);
		printf("ƽ���в%lf\t\t���в%.7lf\n", ��da, ��dm);
		fprintf(PRINT, "\t\t����������%d\t\t��������ֵ:%lf\t\t", Time_total, ��);
		fprintf(PRINT, "ƽ���в%lf\t\t���в%.7lf\n", ��da, ��dm);

		//�ڶ��ֵ���
		�� = 1.375;       //�����Ϊ��������ֵ��
		T = T + 1;
		printf("\n��%d�ֵ���\n", T);
		fprintf(PRINT, "\n��%d�ֵ���\n", T);
		Time_total = 0;
		for (Time = 1; Time < 13; Time++)
		{
			SOR();                            //12�ε���
			Time_total = Time_total + 1;
			if (Time == 11)
				��da1 = ��da;
			else if (Time == 12)
				��da2 = ��da;

		}
		printf("\t\t����������%d\t\t��������ֵ:%lf\t\t", Time_total, ��);
		printf("ƽ���в%lf\t\t���в%.7lf\n", ��da, ��dm);
		fprintf(PRINT, "\t\t����������%d\t\t��������ֵ:%lf\t\t", Time_total, ��);
		fprintf(PRINT, "ƽ���в%lf\t\t���в%.7lf\n", ��da, ��dm);
		�� = ��da2 / ��da1;                        //��в�ƽ��ֵ֮�Ȧ�
		�̦� = (�� + �� - 1) / (sqrt(��) * ��);      //���м�����̦�
		�ئ� = 2 / (1 + sqrt(1 - �̦� * �̦�));    //��������Ӧئ�
		��m = 1.25 * �ئ� - 0.5;                     //�����ѵ�������

		//�����ֵ���
		do
		{
			�� = ��m;                             //�����Ϊ��һ����õ���ѵ�������
			T = T + 1;
			printf("\n��%d�ֵ���\n", T);
			fprintf(PRINT, "\n��%d�ֵ���\n", T);
			Time_total = 0;
			for (Time = 1; Time < 13; Time++)
			{
				SOR();                            //12�ε���
				Time_total = Time_total + 1;
				if (Time == 11)
					��da1 = ��da;
				else if (Time == 12)
					��da2 = ��da;

			}
			printf("\t\t����������%d\t\t��������ֵ:%lf\t\t", Time_total, ��);
			printf("ƽ���в%lf\t\t���в%.7lf\n", ��da, ��dm);
			fprintf(PRINT, "\t\t����������%d\t\t��������ֵ:%lf\t\t", Time_total, ��);
			fprintf(PRINT, "ƽ���в%lf\t\t���в%.7lf\n", ��da, ��dm);
			��m1 = ��m;
			�� = ��da2 / ��da1;
			�̦� = (�� + �� - 1) / (sqrt(��) * ��);
			�ئ� = 2 / (1 + sqrt(1 - �̦� * �̦�));
			��m = 1.25 * �ئ� - 0.5;
		} while (fabs((��m - ��m1) / (2 - ��m1)) >= 0.05);               //������������æ�m1�ͦ�m�����ж�ʽ����̶���m1
		�� = ��m;

		do
		{
			T = T + 1;
			printf("��%d�ֵ���", T);
			fprintf(PRINT, "\n��%d�ֵ���\n", T);
			SOR();//��ѵ�������ȷ�����������
			Time_total = Time_total + 1;
			printf("\t\t����������%d\t\t��������ֵ:%lf\t\t", Time_total, ��);
			printf("ƽ���в%lf\t\t���в%.7lf\n", ��da, ��dm);
			fprintf(PRINT, "\t\t����������%d\t\t��������ֵ:%lf\t\t", Time_total, ��);
			fprintf(PRINT, "ƽ���в%lf\t\t���в%.7lf\n", ��da, ��dm);
			sign = 0;
			for (i = 0; i < row; i++)
				for (j = 0; j < col; j++)
				{
					if (��d[i][j] >= e)//�����ﵽ���ƾ��Ȧ���ֹͣ
						sign = sign + 1;
				}
		} while (sign > 0);
		//������λ�������

		//�������õ��Ŀռ�糡�ֲ�
		printf("\n\t\t\t\t��������ϵĵ�λֵ��\n");
		fprintf(PRINT, "\n\t\t\t\t��������ϵĵ�λֵ��\n");
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
		if (m == 1)//���������Ǽ��ֵ��ֻ��һ�����������ȵ�λֵ
		{
			for (m = 1; a < (V[n - 1] - dv); m++)//�����λ����
				a = a + dv;
			Ve = (double*)calloc((m - 1) * sizeof(double), sizeof(double));
			for (i = 0; Ve[i] < (V[n - 1] - dv); i++)
				Ve[i + 1] = Ve[i] + dv;//�������λ�ߵĵ�λֵ
		}

		//���¶���ͼ�β������λ������
		initgraph(100 + int(10 * ZZ[col - 1]), 100 + int(10 * RR[row - 1]));//����ര�ڱ�Ե���ܾ�50������ֵ������10��
		setbkcolor(BLACK);//�趨��ɫ����
		cleardevice();
		rectangle(50, 50, 50 + int(10 * ZZ[col - 1]), 50 + int(10 * RR[row - 1]));//�ര�ڱ�Ե���ܾ�50������

		dj = (int*)calloc((n + 1) * sizeof(int), sizeof(int));//Ϊ�缫�Ҳ��λ�÷����ڴ�
		dj[0] = 0;//��dj[0]=0
		for (i = 1; i < n + 1; i++)//Ϊdj��ÿ��Ԫ�ظ�ֵ
		{
			if (i == n)
				dj[i] = dj[i - 1] + N[i - 1];//ӫ���������൱�ڵ�n���缫�ұ߽磩���λ��
			else
				dj[i] = dj[i - 1] + N[i - 1] + 1;//һ��缫�ұ߽���λ��
		}

		for (i = 1; i < n; i++)
			bar(50 + int(10 * ZZ[dj[i] - 1]), 50, 50 + int(10 * (ZZ[dj[i]])), 50 + int(10 * r2));   //���缫
		ze = (double*)calloc(row * col * sizeof(double), sizeof(double));   //�ȵ�λ��z����ֵ�����ڴ�
		re = (double*)calloc(row * col * sizeof(double), sizeof(double));   //�ȵ�λ��r����ֵ�����ڴ�
		zeb = (double*)calloc(row * col * sizeof(double), sizeof(double));  //����10���ĵȵ�λ��z����ֵ�����ڴ�
		reb = (double*)calloc(row * col * sizeof(double), sizeof(double));  //����10���ĵȵ�λ��r����ֵ�����ڴ�
		zlm = zz[0];
		for (i = 0; i < col - 1; i++)
		{
			if (zlm < zz[i])
				zlm = zz[i];        //����z�����������
		}
		rlm = rr[0];
		if (rlm < rr[row - 2])
			rlm = rr[row - 2];        //����r�����������
		zrlm = sqrt(zlm * zlm + rlm * rlm);	 //����Խ��߷��������
		for (k = 0; k < m; k++)
		{
			f1 = 0;    //��¼���еȵ�λ��������ʼ��
			l = 0;     //�ȵ�λ��������ʼ��
			for (i1 = 0; i1 < row; i1++)
			{
				for (j1 = 0; j1 < col - 1; j1++)   //z�᷽�����Բ�ֵ����ȵ�λ������
				{
					if (((E[i1][j1] < Ve[k]) & (Ve[k] < E[i1][j1 + 1])) || ((E[i1][j1] > Ve[k]) & (Ve[k] > E[i1][j1 + 1])) || (E[i1][j1] == Ve[k]))
					{
						re[l] = RR[i1];
						if (re[l] == RR[0])
							f1++;          //���б��ֵ��¼���еȵ�λ���������ڿ��Ƶ�λ�߻���
						if (E[i1][j1] == Ve[k])
							ze[l] = ZZ[j1];
						else
							ze[l] = ZZ[j1] + (Ve[k] - E[i1][j1]) * zz[j1] / (E[i1][j1 + 1] - E[i1][j1]);
						reb[l] = 50 + 10 * re[l];             //����10�����ͼ������ 
						zeb[l] = 50 + 10 * ze[l];             //����10�����ͼ������ 
						l++;
					}
				}

				for (j = 0; j < col; j++)        //r�᷽�����Բ�ֵ����ȵ�λ������
				{
					if (i1 < row - 1)
					{
						i = i1;
						if (((E[i][j] < Ve[k]) & (Ve[k] < E[i + 1][j])) || ((E[i][j] > Ve[k]) & (Ve[k] > E[i + 1][j])))
						{
							ze[l] = ZZ[j];
							re[l] = RR[i] + (Ve[k] - E[i][j]) * rr[i] / (E[i + 1][j] - E[i][j]);
							reb[l] = 50 + 10 * re[l];              //����10�����ͼ������
							zeb[l] = 50 + 10 * ze[l];              //����10�����ͼ������
							l++;
						}
					}
				}


			}


			//��λ�߻���
			f = 0;   //�����߶����λ�ñ�����ʼ��
			c = 0;   //�������ӹ��ĵȵ�λ��λ��ѭ��������ʼ��
			f6 = 0;  //���ڼ�¼��λ���ڵ缫��λʱ�ĵ缫λ�ó�ʼ��
			hi = (int*)calloc(l * sizeof(int), sizeof(int));   //�������ӹ��ĵȵ�λ��λ��
			do
			{

				fprintf(PRINT, "\n��λֵΪ%lf�ĵ�λ�ߵĸ�������ֵ��\n", Ve[k]);
				setcolor(20 * k * RED + 10 * k * GREEN + BLUE);    //ȷ����λ����ɫ
				for (i = 0; i < l; i++)
				{
					pan = 1;
					for (b = 0; b < c; b++)    //�ж��Ƿ�Ϊ�����ӹ��ĵ�
					{
						if (i == hi[b])
						{
							pan = 0;
							break;
						}
					}
					if ((pan == 1) & (re[i] == RR[0]))    //�����а���ʱ�ҳ���һ����λ�ߵ����λ��
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
						for (b = 0; b < c; b++)    //�ж��Ƿ�Ϊ�����ӹ��ĵ�
						{
							if (j == hi[b])
							{
								pan = 0;
								break;
							}
						}
						if ((j != i) & (pan == 1) & (sqrt((zeb[j] - zeb[i]) * (zeb[j] - zeb[i]) + (reb[j] - reb[i]) * (reb[j] - reb[i])) < 12 * zrlm) & (fabs(zeb[j] - zeb[i]) < 12 * zlm) & (fabs(reb[j] - reb[i]) < 12 * rlm))
						{    //����ڶ��㲻�ǵ�һ�� �� ���������ӹ��ĵ� �� ����֮���������Ҫ��
							f3 = 0;
							for (f2 = 0; f2 < n; f2++)    //������ڵ缫��λ���¼�缫λ�ò������
							{
								if (Ve[k] == V[f2])
								{
									f3 = 1;    //�����
									f6 = f2;   //��¼�缫λ��
									break;
								}
							}
							if ((f3 == 1) & ((ZZ[dj[f6 + 1]] == ze[j]) || (ZZ[dj[f6 + 1] - 1] == ze[j])) & (RR[0] == re[j]))
							{    //����ǵ缫��λ�ҵڶ���Ϊ�缫��������λ����������¡�Z���������ӵ缫��λ��
								for (f2 = 0; f2 < M2; f2++)
								{
									fprintf(PRINT, "(%lf\t,\t%lf\t)\n", re[i], ze[i]);   //�����һ������
									if (re[i] == RR[0])
										f1--;     //�����һ���ڶ�������ֵ��1
									line(int(zeb[i]), int(reb[i]), int(zeb[j]), int(reb[j]));
									hi[c] = i;      //��¼��һ��λ��
									c++;
									fprintf(PRINT, "(%lf\t,\t%lf\t)\n", re[j], ze[j]);   //����ڶ�������
									if (re[j] == RR[0])
										f1--;	  //����ڶ����ڶ�������ֵ��1			
									for (f4 = 0; f4 < l; f4++)
									{
										if ((ze[f4] == ze[i]) & (re[f4] == re[i] + r2 / M2))
										{
											i = f4;      //Ѱ����һ����һ��
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
											j = f5;      //Ѱ����һ���ڶ���
											break;
										}
									}
								}
								fprintf(PRINT, "(%lf\t,\t%lf\t)\n", re[i], ze[i]);  //����缫�׶�����λ������
								if (re[i] == RR[0])
									f1--;			//����˵��ڶ�������ֵ��1		
								line(int(zeb[i]), int(reb[i]), int(zeb[j]), int(reb[j]));   //���ӵ缫�׶�ˮƽ��
								hi[c] = i;            //��¼�缫�׶�����λ��λ��
								c++;
								break;
							}
							else if ((f3 == 1) & ((zz[dj[f6 + 1]] == ze[j]) || (ZZ[dj[f6 + 1] - 1] == ze[j])) & (RR[0] != re[j]))
							{   //����ǵ缫��λ�ҵڶ��㲻�ǵ缫��������λ����������ϵ���Z���������ӵ缫��λ��
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
							{   //����ڶ��㲻�ǵ缫����λ���������һ�����겢���ӵ�һ����
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
					{   //����ڶ���λ�ڶ��л����������ڶ������겢����ѭ���ж��Ƿ�����һ����λ��
						fprintf(PRINT, "(%lf\t,\t%lf\t)\n", re[j], ze[j]);
						if (re[j] == RR[0])
							f1--;
						hi[c] = j;
						c++;
						break;
					}
					else   //����ڶ��㲻λ�ڶ��л�����򽫵ڶ�����Ϊ��һ��ѭ���ĵ�һ��
						i = j - 1;
				}
			} while (f1 > 0);      //������б��ֵ��Ϊ0���������������һ����λ��
		}




	};



};