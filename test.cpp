#include<iostream>
#include <chrono>
#include <iostream>
#include <vector>
using namespace std;
int main(){
    vector<int> a={1,2,3,4};
    for (int i = a.size()-1 ; i >0; i--)
    {
        for (int j = 0; j < i; j++)
        {
            a[j]=a[j]+a[j+1];
        }
    }
    cout<<a[0]<<endl;
    return 0;
}
