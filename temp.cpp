#include <tuple>

int f(int a, int b) { return a + b;}

int main() {
    //int g(tuple<int, int>);

    //std::tuple<int, int> tup(1, 2);

    //invoke(f, 1, 2); // calls f(1, 2)
    //invoke(g, tup);  // calls g(tup)
;
    std::apply(f, std::pair(1,2));   // also calls f(1, 2)
    return 0;
}

