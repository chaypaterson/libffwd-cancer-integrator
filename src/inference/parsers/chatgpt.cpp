#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>

int main() {
    std::string str("[7, 48, 0, 13, 2]");
    std::vector<int> vect;

    // Remove brackets and spaces from the input string
    str.erase(std::remove_if(str.begin(), str.end(), 
                [](char c) { return c == '[' || c == ']' || c == ' '; }), 
                str.end());

    // Tokenize the string based on ','
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, ',')) {
        // Convert token to integer and store in vector
        vect.push_back(std::stoi(token));
    }

    // Print the vector elements
    for (int num : vect) {
        std::cout << num << " ";
    }

    return 0;
}
