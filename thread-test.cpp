#include <iostream>
#include <thread>

int main() {
    int num_thr = std::thread::hardware_concurrency();

    std::cout << num_thr << std::endl;

    return 0;
}
