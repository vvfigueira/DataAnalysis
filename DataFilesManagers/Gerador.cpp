#include <fstream>

int main(){
    
    std::ofstream ofs;
    ofs.open("Macro.mac", std::ofstream::out | std::ofstream::trunc);
    ofs.close();
    std::ofstream eFile ("Macro.mac" ,std::ofstream::app);

    double pressInit = 133.322;
    double passoP = 133.322/5;
    int repetP = 65;
    bool Pconst = true;

    double tempInit = 0;
    double passoT = 50;
    int repetT = 12;
    bool Tconst = true;

    double campoInit = 12500;
    double passoC = 500;
    int repetC = 25;
    bool Cconst = true;

    double enerInit = 2;
    double passoE = 2;
    int repetE = 199;
    bool Econst = true;

    int nPart = 200000;

    eFile << "/control/verbose 0\n"
          << "/process/em/verbose 0\n"
          << "/process/verbose 0\n"
          << "/process/had/verbose 0\n"
          << "/run/verbose 1\n\n" 
          << "/run/initialize\n\n"
          << "/gps/particle gamma\n"
          << "/gps/position -0.5 0 0 m\n"
          << "/gps/direction 1 0 0\n"
          << "/gps/ene/type Mono\n\n";

    for(int i = 0; i<= repetE; i++){
        eFile << "/gps/ene/mono " << enerInit +passoE*i<< " eV\n"
        
              //<< "/run/reinitializeGeometry\n"
              << "/run/beamOn " << nPart <<"\n\n";
    }
    eFile.close();

    return 0;
}