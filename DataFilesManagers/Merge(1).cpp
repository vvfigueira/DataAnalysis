// Code for merging various files in one. Must be used with `Merge.sh`

#include <iostream>
#include <fstream>

int main(int argv, char** argc){

    std::ifstream filein;
    std::fstream fileout;

    filein.open(argc[1], std::ios::in);

    fileout.open(argc[2], std::ios::app);

    double x, y;

    while(!filein.eof()){
    
        filein >> x >> y;

        fileout << x <<"\t"<< y << "\n";

    }

    fileout.close();

    filein.close(); 
    
    return 0;
}