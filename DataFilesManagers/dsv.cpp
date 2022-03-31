#include <fstream>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <tuple>

double Mean(double data[], int k) {
    double sum = 0;

    for(int i = 0; i < k+1; i++) {
        sum += data[i];
    }

    double dsv = sum/(k+1);

    return dsv;
}

double Dsv(double data[], int k) {
    double sum = 0, mean, dsv = 0;

    for(int i = 0; i < k+1; ++i) {
        sum += data[i];
    }

    mean = sum / (k+1);

    for(int i = 0; i < k+1; ++i) {
        dsv += pow(data[i] - mean, 2);
    }

    return sqrt(dsv / (k+1));
}

int main(int argc,char** argv){

    switch (argc){
        case 3:
            
            break;
        default:
            std::cout << "\nUso: ./dsv [Arquivo_Entrada] [Arquivo_Saida]\n";
            return 0;
            break;
    }

    std::fstream filein, fileout;

    filein.open(argv[1], std::ios::in);
    fileout.open(argv[2], std::ios::out);

    double x[100000];
    double y[2];  
    double temp;
    int i = 0;

    filein >> y[0] >> y[1];
    temp = y[1];
    x[0] = y[0];

    while(!filein.eof()){

        filein >> y[0] >> y[1];

        if(y[1] == temp){
            i++;
            x[i] = y[0];
            if(filein.eof()) fileout << Mean(x, i) << "\t" << Dsv(x, i) << "\t" << temp <<"\n";
        }else{
            fileout << Mean(x,i) << "\t" << Dsv(x, i) << "\t" << temp <<"\n";
            i = 0;
            temp = y[1];
            x[0] = y[0];
            if(filein.eof()) fileout << Mean(x,i) << "\t" << Dsv(x,i) << "\t" << temp <<"\n";
        }
    }

    filein.close();
    fileout.close();

    return 0;
}