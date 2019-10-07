Introduction
============
This project provides a C implementation (and a Python wrapper) encoding/recoding/decoding for network coding. Random linear network coding (RLNC) [1] and several of its sparse variants are implemented. The library supports three catagories of sparse network codes (SNC): random codes, band codes and BATS-like codes [2].

The following decoders are implemented: 

- (sub)generation-by-(sub)generation (GG) decoder has a linear-time complexity but exhibits higher decoding-induced overhead.

- overlap-aware (OA) decoder has optimized code overhead but exhibits higher complexity.

- band (BD) and compact band (CBD) decoders that apply to band code only.

- PP(perpetual) decoder applies to decode the window-wrapped codes only.

The library had been used for performance evaluation and comparison when authoring the following papers. If you find the lib useful, please cite the papers when appropriate.

- Ye Li, J. Zhu and Z. Bao, "Sparse Random Linear Network Coding With Precoded Band Codes," in IEEE Communications Letters, vol. 21, no. 3, pp. 480-483, March 2017.

- Ye Li, S. Zhang, J. Wang, X. Ji, H. Wu and Z. Bao, "A Low-Complexity Coded Transmission Scheme over Finite-Buffer Relay Links," in IEEE Transactions on Communications, vol. 66, no. 7, pp. 2873 - 2887, July. 2018.

- Ye Li, S. Zhang, J. Wang, X. Ji, H. Wu and Z. Bao, "A Low-Complexity Coded Transmission Scheme over Finite-Buffer Relay Links," in IEEE Transactions on Communications, vol. 66, no. 7, pp. 2873 - 2887, July. 2018.

- Y. Li, J. Wang, S. Zhang, Z. Bao and J. Wang, "Efficient Coastal Communications with Sparse Network Coding," in IEEE Network, vol. 32, no. 4, pp. 122-128, July/August 2018.

Usage
============
The library is available as a shared library which is compiled by

```shell
$ make libsparsenc.so
```

Accessing API is via `include/sparsenc.h`. 

Some examples are provided to test the codes and decoders (see examples/ directory). Run

```shell
$ make sncDecoders
```

and test using

```shell
usage: ./sncDecoders code_t dec_t datasize pcrate size_b size_g size_p bpc bnc sys
                       code_t   - RAND, BAND, WINDWRAP, BATS
                       dec_t    - GG, OA, BD, CBD, PP
                       datasize - Number of bytes
                       pcrate   - Precode rate (percentage of check packets)
                       size_b   - Subgeneration distance
                       size_g   - Subgeneration size
                       size_p   - Packet size in bytes
                       bpc      - Use binary precode (0 or 1)
                       bnc      - Use binary network code (0 or 1)
                       sys      - Systematic code (0 or 1)
```

Please note that currently only the OA decoder is tested for decoding the BATS code.

To test the code over example networks, run

```
$ make sncRecoders-n-Hop
```
for using random/band codes over an n-hop line network where intermediate nodes perform on-the-fly recoding, and 

```
$ make sncRecoderNhopBATS
```
for testing the BATS code over the n-hop line network.

Use

```
$ make sncRecoderFly
```
for sending the random/band codes over a butterfly network. Please see main functions under `examples/xxx.c` for details and `makefile` for other available examples.

Limitation
============
The library only supports coding against a given block of source packets, i.e., a *generation* of packets as termed in the network coding literature. Sliding-window mode is not supported.

Reference
============
[1] T. Ho, M. Medard, R. Koetter, and D. R. Karger, "A Random Linear Network Coding Approach to Multicast", in IEEE Transactions on Information Theory, Vol 52, No. 10, pp. 4413-4430, 2006.

[2] S. Yang and R. W. Yeung, "Batched Sparse Codes," in IEEE Transactions on Information Theory, vol. 60, no. 9, pp. 5322-5346, Sept. 2014.