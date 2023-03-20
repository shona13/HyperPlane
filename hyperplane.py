from functools import reduce
import gmpy2 as gmp
from random import randint
import timeit


'''
    Author : Sonal Joshi
    Description : Hyperplane classification using Paillier for encryption scheme 
'''


# Generates arrays using list - Takes parameter l & d as user input
def arr_gen(l, d):
    # 1d array (list) / vector with d elements
    user = [gmp.mpz(randint(1, 100)) for j in range(d)]
    # 2d array (list) / matrix with l x d dimension
    cloud = [[gmp.mpz(randint(1, 100)) for j in range(d)] for i in range(l)]
    return user, cloud


# Key generation - Takes parameters prime p, generator g and common modulus n
def key_gen(p, q, n):
    # Calculates g
    g = n + 1
    # Calculates lambda
    lm = (p - 1) * (q - 1)
    return g, lm


# Randomly generates random number r - Takes parameter number of bits nbit
def random_num(nbit):
    state = gmp.random_state(hash(gmp.random_state()))
    r = gmp.mpz_random(state, nbit)
    # Checking for gcd of (r,n) = 1 i.e Check if r belongs to Zn*
    while gmp.gcd(r, n) != gmp.mpz(1):
        r = gmp.mpz_random(state, 2 ** nbit)
    return r


# Encryption function - Takes parameters  generator g, composite modulus n, and arrays m (user) and c (cloud) , l and d
def encrypt(g, n, m, c, l, d):
    # Initializing array new_cipher
    new_cipher = []
    # Creating an empty array mul_cipher with 0
    mul_cipher = [[0 for j in range(d)] for i in range(l)]
    # Creates a random number which is less than 2^256
    r = random_num(2 ** 256)
    # n^2 variable
    square_n = gmp.mpz(n * n)
    # User cipher E(x + r) which is a vector
    user_cipher = [gmp.mod(gmp.mul(gmp.powmod(g, gmp.mpz(m[j]+r), square_n), gmp.powmod(r, n, square_n)), square_n) for j
                   in
                   range(d)]
    # Cloud cipher E(W) which is a 2d matrix
    cloud_cipher = [
        [gmp.mod(gmp.mul(gmp.powmod(g, gmp.mpz(c[i][j]), square_n), gmp.powmod(r, n, square_n)), square_n) for j in
         range(d)] for i in range(l)]
    # Using constant multiplication property of homomorphic paillier encryption scheme
    mul_cipher = [[gmp.mpz(user_cipher[i] ** cloud_cipher[j][i]) for i in range(d)] for j in range(l)]
    # Multiplying every element in a paticular row and reducing it to 1d array
    for lists in mul_cipher:
        new_cipher.append(reduce((lambda x, y: x * y), lists))
    return user_cipher, cloud_cipher, new_cipher


# Decryption function - Takes parameters new_cipher as nc, composite modulus n & lambda function lm
def decrypt(nc, n, lm):
    # Initializing decm as array (list)
    decm = []
    # n^2 variable
    square_n = gmp.mpz(n * n)
    # Computes Inverse of lambda function
    inv_lm = gmp.invert(lm, n)
    # Decrypting every encrypted element in the 1d array & appending it to decm
    for cipher in nc:
        x = gmp.sub(gmp.powmod(cipher, lm, square_n), gmp.mpz(1))
        decm.append(gmp.mod(gmp.mul(gmp.f_div(x, n), inv_lm), n))
    # After decryption, calculating the argmax of the elements i.e index of maximum element in the array
    argmax = decm.index((max(decm)))
    return decm, argmax


# Main funcion
if __name__ == '__main__':
    #l = 2
    #d = 3
    # Taking user input
    l = int(input("Enter l (>=2): "))
    d = int(input("Enter d (>=2): "))
    # Defining the values for p,q,n
    p = 7
    q = 11
    n = 77

    '''
        Calling the functions and creating objects for each function
    '''
    # Array generation object
    user, cloud = arr_gen(l, d)
    # Key generation object
    g, lm = key_gen(p, q, n)
    # Encryption of user cipher, cloud cipher & new_cipher objects
    uc, cc, nc = encrypt(g, n, user, cloud, l, d)
    # Decryption & argmax (t) object
    decl, argmax = decrypt(nc, n, lm)

    '''
    Printing values on terminal    
    '''

    ### Encryption process ###
    print("-"*20, "Encryption", "-"*20)
    print(f"\nEncrypted vector E(x) = {uc}\nEncrypted result (C1,C2,..Cl) = {nc}\n\n")

    ### Decryption process ###
    print("-" * 20, "Decryption", "-" * 20)
    print(f"The input x belongs to the class t = {argmax}")

    ### Time taken to execute the program ###
    print("\nTime taken to run (seconds): ", timeit.timeit())

