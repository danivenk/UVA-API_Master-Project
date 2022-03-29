import numpy as np
import mono_lags.um21_lagmodel as um

def main():
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

if __name__ == "__main__":
    main()