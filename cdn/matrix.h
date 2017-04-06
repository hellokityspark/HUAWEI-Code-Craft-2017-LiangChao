#include <iostream>
#include <vector>
#include <cassert>

template <typename T>
class matrix {
public:
    matrix(int m, int n):_m(m), _n(n) {
        init();
    }
    matrix(int m):_m(m), _n(m) {
        init();
    }
    void init() {
        _size = _m * _n;
        _elem = new T[_size];
        for (int i = 0; i < _size; ++i) {
            _elem[i] = 0;
        }
    }
    void pageRank(std::vector<double>& pr) {
        std::vector<double> tmp(_m, 0);
        double loss;
        while (true) {
            loss = 0;
            for (int i = 0; i < _m; ++i) {
                tmp[i] = 0;
                for (int j = 0; j < _m; ++j) {
                    tmp[i] += (*this)(i, j) * pr[j];
                }
            }

            for (int i = 0; i < _m; ++i) {
                loss += (pr[i] - tmp[i]) * (pr[i] - tmp[i]);
                pr[i] = tmp[i];
            }
            if (loss < 1e-6 * _m) break;
        }
    }
    const T& operator()(int m, int n) const {
        return _elem[m * _n + n];
    }
    T& operator()(int m, int n) {
        return _elem[m * _n + n];
    }
    inline int rows() const { return _m;}
    inline int cols() const { return _n;}
    friend std::ostream& operator<<(std::ostream& out, const matrix<T>&  mat) {
        for (int i = 0; i < mat._m; ++i) {
            for (int j = 0; j < mat._n; ++j) {
                out <<mat(i, j) << " ";
            }
            out << std::endl;
        }
        return out;
    }
private:
    T* _elem;
    int _size;
    int _m, _n;
};

// template <typename T>
// std::ostream& operator<<(std::ostream& out, const matrix<T>& mat) {
//     for (int i = 0; i < mat._m; ++i) {
//         for (int j = 0; j < mat._n; ++j) {
//             out <<mat(i, j) << " ";
//         }
//         out << std::endl;
//     }
//     return out;
// }