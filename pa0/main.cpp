#include<cmath>
#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Dense>
#include<iostream>
using namespace std;
using namespace Eigen;
int main(){
    /*
    cout << "Example of cpp \n";
    float a = 1.0, b = 2.0;
    cout << a << endl;
    cout << a/b << endl;
    cout << sqrt(b) << endl;
    cout << acos(-1) << endl;
    cout << sin(30.0/180.0*acos(-1)) << endl;

    cout << "Example of vector \n";
    Vector3f v(1.0f,2.0f,3.0f);
    Vector3f w(1.0f,0.0f,0.0f);
    cout << "Example of output \n";
    cout << v << endl;
    // vector add
    cout << "Example of add \n";
    cout << v + w << endl;
    // vector scalar multiply
    cout << "Example of scalar multiply \n";
    cout << v * 3.0f << endl;
    cout << ""<< endl;
    cout << 2.0f * v << endl;

    // Example of matrix
    cout << "Example of matrix \n";
    // matrix definition
    Matrix3f i,j;
    i << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;
    j << 2.0, 3.0, 1.0, 4.0, 6.0, 5.0, 9.0, 7.0, 8.0;
    // matrix output
    cout << "Example of output \n";
    cout << i << endl;
    // matrix add i + j
    // matrix scalar multiply i * 2.0
    // matrix multiply i * j
    // matrix multiply vector i * v
    cout<<endl;*/
    Vector2f P(2,1),i(1,2);
    Matrix2f r;
    float theta=45;
    r<<cos(acos(-1)*theta/180),-sin(acos(-1)*theta/180),sin(acos(-1)*theta/180),cos(acos(-1)*theta/180);
    cout<<"P:"<<endl;
    cout<<P<<endl;
    cout<<"Rotated 45:"<<endl;
    P=r*P;
    cout<<P<<endl;
    cout<<"Paced (1,2):"<<endl;
    P=P+i;
    cout<<P<<endl;
    return 0;
}