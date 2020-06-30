#include <math.h>

#include <iostream>
#include <queue>
#include <vector>
#include <algorithm>
#include <list>


std::vector<unsigned int> foo(const std::vector<double> v, size_t k) {
    std::priority_queue<std::pair<double, size_t>> queue;
    for (auto i = v.size(); i; --i) {
        double fabs_v = fabs(v[i]);
        if (fabs_v >= 1e-4) {
            queue.push(std::pair<double, size_t>(fabs_v, i));
        }
    }

    size_t bd = std::min(k, queue.size());
    std::vector<unsigned int> idx(bd);
    for (auto i = 0; i < bd; ++i) {
        idx[i] = queue.top().second;
        queue.pop();
    }

    return idx;
}

int main() {

    std::vector<double> v(100);
    for(auto i = 0; i < 100; ++i) {
        v[i] = -i * (100.5 - i);
    }

    for (auto i: foo(v, 10)) {
        std::cout << i << std::endl;
    }
    

    std::vector<double> u{};
    std::cout << u.size() << std::endl;




    printf("%ld\n", sizeof(unsigned int));
    printf("%ld\n", sizeof(int));
    printf("%ld\n", sizeof(size_t));
    printf("%ld\n", sizeof(ssize_t));


    

    return 0;
}


