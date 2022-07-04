import numpy as np
from scipy.fft import fft
import mono_lags.um21_lagmodel as um
from matplotlib import pyplot

def main():
    lor_psd = um.lorentz_q(np.array([i/10 for i in range(100)]), np.array([.3 for _ in range(10)]), np.array([.6 for _ in range(10)]), np.array([.8 for _ in range(10)]))

    a = np.zeros(100)

    print(lor_psd + a)

    # print(np.square([i/10 for i in range(100)]))

    # b = 1+2j

    # print(b, np.conj(b) * b)


def test1():
    a = np.arange(0,10)
    b = []
    primes = [2]

    for value in a:
        if value > 2 and prime(value, primes):
            b.append(True)
            primes.append(value)
        elif value == 2:
            b.append(True)
        else:
            b.append(False)

    b = np.array(b)
    primes = np.array(primes)

    print(a)
    print(a[b])

    c = np.reshape(a, (len(a), 1)) * np.reshape(a[b], (1, len(a[b])))

    print(c)
    print(c[1][3])


def prime(number, primes):

    if number in primes:
        return False

    for prime in primes:
        if number % prime == 0:
            return False

    return True


def test2():
    x = np.linspace(0, 1, 1000)

    y = np.sin(x)
    fft_y = fft(y)


if __name__ == "__main__":
    main()