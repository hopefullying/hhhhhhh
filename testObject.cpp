#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <conio.h>
#include<graphics.h>
#include <stdbool.h>

#include "test.cpp"


using namespace std;
int main() {

	FILE* data;
	FILE* out;


	fopen_s(&data, "E:\\1120160947\\1120160947.dat", "r");
	fopen_s(&out, "E:\\1120160947\\2.res", "w");
	
	test hhh;

	hhh.Myset(data,out);
	hhh.init();

	hhh.calculate();
	hhh.painting();


	fclose(data);//关掉读入的目标文件//	
	fclose(out);
	system("pause");

}

