#include "functions.h"

template <class T>
ostream& operator<<(ostream& os, const list<T> list);

int main() {
    1;

    return 0;
}

template <class T>
ostream& operator<<(ostream& os, const list<T> list) {
    for (auto& i: list) {
        os << i << " ";
    }

    return os;
}