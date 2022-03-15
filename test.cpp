#include<cmath>
#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Dense>
#include<iostream>
using namespace std;
using namespace Eigen;
int main(){
    Vector3f P(-1,2,0),i(-2,1,0),l(1,1,0);
    cout<<P.dot(l)<<endl;
    cout<<endl;
    cout<<i.dot(l)<<endl;
    return 0;
}