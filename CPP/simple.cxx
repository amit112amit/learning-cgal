#include <iostream>
#include <set>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <ctime>

int main(){
    clock_t t = clock();
    for(auto i = 0; i < 100000; ++i){
        std::set<unsigned> index{0,1,2};
        unsigned a = std::rand()%3;
        index.erase(a);
        unsigned b = std::rand()%3;
        index.erase(b);
        std::cout<< *index.begin() << std::endl;
    }
    std::cout<< "Time = " << ((float)clock() - (float)t)/CLOCKS_PER_SEC << std::endl;
}
