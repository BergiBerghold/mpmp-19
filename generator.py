from sys import argv, stdout
from bitarray import bitarray
from time import time
import datetime

# Credit to https://codereview.stackexchange.com/users/210384/greg-ames for the very efficient prime number generator
# Counts and optionally prints all prime numbers no larger than 'n'  

#CUTOFF      = 10          # for debugging only
#SIEVE_SIZE  = 2           # for debugging only
CUTOFF      = 1e4
SIEVE_SIZE  = 2**20
GHz         = 3.5         # on my i5-6285U laptop

# mod 30 wheel constant arrays
modPrms     = [7,11,13,17,19,23,29,31]
modPrmsM30  = [7,11,13,17,19,23,29,1]
gaps        = [4,2,4,2,4,6,2,6,4,2,4,2,4,6,2,6] # 2 loops for overflow
ndxs        = [0,0,0,0,1,1,2,2,2,2,3,3,4,4,4,4,5,5,5,5,5,5,6,6,7,7,7,7,7,7]
rnd2wh      = [7,7,0,0,0,0,0,0,1,1,1,1,2,2,3,3,3,3,4,4,5,5,5,5,6,6,6,6,6,6]

def num2ix(n):
    """Return the wheel index for n."""
    n = n - 7              # adjust for wheel starting at 1st prime past 2,3,5 vs. 0
    return (n//30 << 3) + ndxs[n % 30]       

def ix2num(i):
    """Return a number matching i (a wheel index)."""
    return 30 * (i >> 3) + modPrms[i & 7]   

def progress(j, num_loops, enabled):
    """Display a progress bar on the terminal."""
    if enabled:
        size = 60
        x = size*j//num_loops
        print("%s[%s%s] %i/%i\r" % ("Sieving: ", "#"*x, "."*(size-x), j, num_loops), end=' ')
        stdout.flush()

def prime_gen_wrapper(n):
    """Decide whether to use the segmented sieve or a simpler version.  Stops recursion."""
    if n < CUTOFF:
        return smallSieve(n+1) # rwh1 returns primes < N.  We need sieving primes <= sqrt(limit)
    else:
        return segmentedSieve(n)

def smallSieve(n):
    """Returns a list of primes less than n."""
    # a copy of Robert William Hanks' rwh1 used to get sieving primes for smaller ranges
    # https://stackoverflow.com/questions/2068372/fastest-way-to-list-all-primes-below-n-in-python/3035188#3035188
    sieve = [True] * (n//2)
    for i in range(3,int(n**0.5)+1,2):
        if sieve[i//2]:
            sieve[i*i//2::i] = [False] * ((n-i*i-1)//(2*i)+1)
    return [2] + [2*i+1 for i in range(1,n//2) if sieve[i]]

def segmentedSieve(limit, statsOnly = False): 
    """
    Sieves potential prime numbers up to and including limit.

    statsOnly (default False) controls the return.
        when False, returns a list of primes found.
        when True, returns a count of the primes found.
    """
    # segmentation originally based on Kim Walisch's simple C++ example of segmantation found here 
    #     https://github.com/kimwalisch/primesieve/wiki/Segmented-sieve-of-Eratosthenes

    # mod 30 wheel factorization based on a non-segmented version found here in a comment by Willy Good
    # https://programmingpraxis.com/2012/01/06/pritchards-wheel-sieve/  

    sqrt = int(limit ** 0.5)
    lmtbf = SIEVE_SIZE * 8
    while (lmtbf >> 1) >= limit:
       lmtbf >>= 1         # adjust the sieve size downward for small N

    multiples = []; wx = []
    outPrimes = [2,3,5]    # the wheel skips multiples of these, but they may be needed as output
    count = len(outPrimes)
    lim_ix = num2ix(limit)
    buf = bitarray(lmtbf)
    show_progress = False
    if statsOnly:   # outer loop?
        print("sieve size:", end=' ')
        ss = len(memoryview(buf))
        if ss > 1024:
            print(ss//1024, "KB")
        else:
            print(ss, "bytes")
        if limit > 1e8:
            show_progress = True

    num_loops = (lim_ix + lmtbf - 1)//(lmtbf)   # round up

    # get sieving primes recursively, skipping those eliminated by the wheel
    svPrimes  = prime_gen_wrapper(sqrt)[count:]

    for lo_ix in range(0, lim_ix + 1, lmtbf):   # loop over all the segments
        low = ix2num(lo_ix)
        high = ix2num(lo_ix + lmtbf) - 1
        buf.setall(True)
        progress(lo_ix//(lmtbf), num_loops, show_progress)

        # generate new multiples of sieving primes and wheel indices needed in this segment
        for p in svPrimes[len(multiples):]:
            pSquared = p * p
            if pSquared > high:
                break
            multiples.append(pSquared)
            wx.append(num2ix(p) & 7)

        # sieve the current segment
        for x in range(len(multiples)):
            s  = multiples[x]
            if s <= high:
                p  = svPrimes[x]
                ci = wx[x]
                s -= 7
                p8 = p << 3
                for j in range(8):
                    c = (s//30 << 3) + ndxs[s % 30] - lo_ix
                    # buf[c::p8] = False * ((lmtbf - c) // p8 + 1)
                    buf[c::p8] = False              # much simpler with bitarray vs. pure python
                    s += p * gaps[ci]; ci += 1

        # calculate the next multiple of p to sieve in an upcoming segment and its wheel index
                f       = (high + p - 1)//p         # next factor of a multiple of p past this segment
                f_mod   = f % 30
                i = rnd2wh[f_mod]                   # round up to next wheel index to eliminate multiples of 2,3,5
                nxt = p * (f - f_mod + modPrmsM30[i])   # back to a normal multiple of p past this segment
                wx[x] = i                               # save wheel index
                multiples[x] = nxt                      #                  ... and next multiple of p

        # handle any extras in the last segment
        if high > limit:
            top = lim_ix - lo_ix
        else:
            top = lmtbf -1

        # collect results from this segment
        if statsOnly:
            count += buf[:top+1].count()
        else:
            for i in range(top + 1):
                if buf[i]:
                    x = i + lo_ix
                    p = 30 * (x >> 3) + modPrms[x & 7]   # ix2num(x) inlined, performance is sensitive here
                    outPrimes.append(p)

    if show_progress:
        progress(num_loops, num_loops, True)
        print()

    if statsOnly:
        return count
    else:
        return outPrimes

def segmentedSieveGenerator(limit, statsOnly = False):
    """
    Sieves potential prime numbers up to and including limit.

    statsOnly (default False) controls the return.
        when False, returns a list of primes found.
        when True, returns a count of the primes found.
    """
    # segmentation originally based on Kim Walisch's simple C++ example of segmantation found here
    #     https://github.com/kimwalisch/primesieve/wiki/Segmented-sieve-of-Eratosthenes

    # mod 30 wheel factorization based on a non-segmented version found here in a comment by Willy Good
    # https://programmingpraxis.com/2012/01/06/pritchards-wheel-sieve/

    sqrt = int(limit ** 0.5)
    lmtbf = SIEVE_SIZE * 8
    while (lmtbf >> 1) >= limit:
       lmtbf >>= 1         # adjust the sieve size downward for small N

    yield 2
    yield 3
    yield 5

    multiples = []; wx = []
    outPrimes = [2,3,5]    # the wheel skips multiples of these, but they may be needed as output
    count = len(outPrimes)
    lim_ix = num2ix(limit)
    buf = bitarray(lmtbf)
    show_progress = False
    if statsOnly:   # outer loop?
        print("sieve size:", end=' ')
        ss = len(memoryview(buf))
        if ss > 1024:
            print(ss//1024, "KB")
        else:
            print(ss, "bytes")
        if limit > 1e8:
            show_progress = True

    num_loops = (lim_ix + lmtbf - 1)//(lmtbf)   # round up

    # get sieving primes recursively, skipping those eliminated by the wheel
    svPrimes  = prime_gen_wrapper(sqrt)[count:]

    for lo_ix in range(0, lim_ix + 1, lmtbf):   # loop over all the segments
        low = ix2num(lo_ix)
        high = ix2num(lo_ix + lmtbf) - 1
        buf.setall(True)
        progress(lo_ix//(lmtbf), num_loops, show_progress)

        # generate new multiples of sieving primes and wheel indices needed in this segment
        for p in svPrimes[len(multiples):]:
            pSquared = p * p
            if pSquared > high:
                break
            multiples.append(pSquared)
            wx.append(num2ix(p) & 7)

        # sieve the current segment
        for x in range(len(multiples)):
            s  = multiples[x]
            if s <= high:
                p  = svPrimes[x]
                ci = wx[x]
                s -= 7
                p8 = p << 3
                for j in range(8):
                    c = (s//30 << 3) + ndxs[s % 30] - lo_ix
                    # buf[c::p8] = False * ((lmtbf - c) // p8 + 1)
                    buf[c::p8] = False              # much simpler with bitarray vs. pure python
                    s += p * gaps[ci]; ci += 1

        # calculate the next multiple of p to sieve in an upcoming segment and its wheel index
                f       = (high + p - 1)//p         # next factor of a multiple of p past this segment
                f_mod   = f % 30
                i = rnd2wh[f_mod]                   # round up to next wheel index to eliminate multiples of 2,3,5
                nxt = p * (f - f_mod + modPrmsM30[i])   # back to a normal multiple of p past this segment
                wx[x] = i                               # save wheel index
                multiples[x] = nxt                      #                  ... and next multiple of p

        # handle any extras in the last segment
        if high > limit:
            top = lim_ix - lo_ix
        else:
            top = lmtbf -1

        # collect results from this segment
        if statsOnly:
            count += buf[:top+1].count()
        else:
            for i in range(top + 1):
                if buf[i]:
                    x = i + lo_ix
                    p = 30 * (x >> 3) + modPrms[x & 7]   # ix2num(x) inlined, performance is sensitive here
                    yield p


# Driver Code
sieve_limit = 10**15
primes_counter = 0
start = time()

with open("generator_log.txt", "a") as f:
    f.write('\nNew Run starting at ' + str(datetime.datetime.now()) + '\n')

for idx, prime in enumerate(segmentedSieveGenerator(sieve_limit)):
    primes_counter += prime * prime

    if primes_counter % (idx+1) == 0:
        print(idx+1)
        with open("generator_log.txt", "a") as f:
            f.write(str(idx+1) + '\n')

    if (idx+1) % 10000000 == 0:
        print('Checked n to ' + str(idx+1) + ' in ' + str(time() - start))
        with open("generator_log.txt", "a") as f:
            f.write('Checked n to ' + str(idx+1) + ' in ' + str(time() - start) + '\n')
