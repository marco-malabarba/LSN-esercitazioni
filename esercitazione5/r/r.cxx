#include <string>

using namespace std;

int main() {
string command = "python ";
string filename = "~/Scrivania/programmi/r/r.py";
command += filename;
system(command.c_str());

command="g++ main.cxx -o main.exe";
system(command.c_str());

command="./main.exe";
system(command.c_str());

command = "python ";
filename = "~/Scrivania/programmi/r/disegna.py";
command += filename;
system(command.c_str());

return 0;
}
