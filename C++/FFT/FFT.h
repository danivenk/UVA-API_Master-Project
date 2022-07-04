#include "pocketfft/pocketfft_hdronly.h"
#include "../includes.h"

using namespace pocketfft;

template<class T, class K>
list<complex<K>> FFT(T in) {

    vector<K> in_vec(in.begin(), in.end());

    shape_t shape{in_vec.size()};
    stride_t stride_in(shape.size());
    stride_t stride_out(shape.size());
    size_t tmp_in = sizeof(K);
    size_t tmp_out = sizeof(complex<K>);

    for (int i=shape.size()-1; i>=0; --i)
    {
        stride_in[i]=tmp_in;
        stride_out[i]=tmp_out;
        tmp_in*=shape[i];
        tmp_out*=shape[i];
    }

    size_t ndata = in_vec.size();

    for (size_t i=0; i<shape.size(); ++i) {
        ndata*=shape[i];
    }

    shape_t axes;
    for (size_t i=0; i<shape.size(); ++i) {
        axes.push_back(i);
    }

    vector<complex<K>> out_vec(in_vec.begin(), in_vec.end());

    r2c(shape, stride_in, stride_out, axes, FORWARD, in_vec.data(), out_vec.data(), 1.);

    list<complex<K>> out;

    size_t i = 0;
    size_t lim = in_vec.size()/2 + 1;

    for (auto& out_val: out_vec) {
        if (i != 0 && i < lim) {
            out.push_back(out_val);
        }
        i++;
    }

    return out;
}