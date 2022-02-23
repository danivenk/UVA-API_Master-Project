import numpy as np
import mono_lags.um21_lagmodel as um

def main():
    um.sphere(2, [1, 10, 10])

    
    ntheta = 2
    nphi = 3
    dtheta = np.pi/(2.*ntheta)
    dphi = 2.*np.pi/nphi
    theta = np.arange(0,np.pi/2.,dtheta)
    phi = np.arange(0,2*np.pi,dphi)
    phi_arr, theta_arr = np.meshgrid(phi, theta)

    print("theta", theta_arr, sep="\n")
    print("phi", phi_arr, sep="\n")

    # a = np.arange(0,100)
    # b = []
    # primes = [2]

    # for value in a:
    #     if value > 2 and prime(value, primes):
    #         b.append(True)
    #         primes.append(value)
    #     else:
    #         b.append(False)

    # b = np.array(b)

    # a = a.reshape(10,10)
    # b = b.reshape(10,10)

    # print(primes)
    # print(a)
    # print(b)
    # print(a[b])

def prime(number, primes):

    if number in primes:
        return False

    for prime in primes:
        if number % prime == 0:
            return False

    return True

if __name__ == "__main__":
    main()