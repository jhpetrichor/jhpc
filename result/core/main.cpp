#include <iostream>
#include <fstream>
#include <sstream>
#include <set>

using namespace std;

int main(int argc, char *argv[]) {
    string filenome = argv[1];
    string ofile = argv[2];

    fstream file(filenome, ios::in);
    string lien;
    int i = 0;

    set<set<string>> complexes;
    set<string> complex;
    while (getline(file, lien)) {
        i += 1;
        cout << i <<endl;
        if(i % 4 == 0 || i % 4 == 1) continue;
        istringstream iss(lien);
        string protein;
        iss >> protein;
        while(iss >> protein) {
            complex.insert(protein);
        }
        if(i % 4 == 3) {
            complexes.insert(complex);
            complex.clear();
            cout <<  complexes.size() << "\n";
        }
    }

    ofstream of(ofile, ios::out);
    for(auto complex : complexes) {
        for(auto protein : complex) {
            of << protein << "\t";
        }
        of << "\n";
    }

}
