Introduction
============
This project provides a set of APIs to encode and decode **sparse** **n**etwork **c**odes. Network coding [1][2] is known to be able to improve network throughput and combat packet loss by algebrically combining packets at intermediate nodes. This library provides several computationally efficient network codes to overcome the high complexity issue of the random linear network coding (RLNC). The basic idea is to encode against subsets of source packets. In the very early development of this library, the subsets were referred to as generations [3]. In the current version, the subsets are referred to as **subgenerations** because in network coding literature the term "generation" more often refers to the whole set of source packets to be coded by RLNC. Using subgenerations is to 1) avoid naming conflict and 2) to emphasize the subset concept. The RLNC is a special case that the number of subgenerations is one.

The library at present supports three catagories of sparse network codes (SNC): random codes, band codes and BATS-like codes. Both the random and band SNCs have a (pre-determined) fix number of subgenerations. The random SNC [4] encodes from subgenerations that are randomly overlapped whereas the band code [5][6] encodes from subgenerations that are overlapped consecutively (in terms of packet indices). The _band_ name comes from that its decoding matrix is a band matrix. A variant of the band codes called window-wrapped (or perpetual) codes [7] is also supported. The difference is that the window-wrapped code allows wrap-around in some encoding vectors and therefore the decoding matrix is not ideally banded. The third type of SNC supported is the fixed-degree BATS code [8] which encodes from randomly sampled subsets. Unlike the random and band codes, the number of sampled subsets is not pre-determined but depends on the channel quality. The benefits of the code is the buffer size required at intermediate nodes can be significantly reduced. In contrast to the original BATS code [9] which uses a degree distribution to detrmine the subset sizes, the BATS code implemented in the library has a fixed degree (i.e., subset size) and the number of coded packets sent from each batch (called batch transmission size, BTS) can be designed (analytically or heuristically) according to the link qualities.

Five decoders with different performance tradeoff considerations are implemented in the library. The (sub)generation-by-(sub)generation (GG) decoder has a linear-time complexity but exhibits higher decoding-induced overhead. On the other hand, the overlap-aware (OA) decoder [10] has optimized code overhead but exhibits higher complexity. These two decoders essentially can be used for all kinds of subgeneration-based code (i.e., not limited to the three codes currently provided in the library). The band (BD) and compact band (CBD) decoders can only be applied to the band code. The decoders have optimized code overheads as well and their complexities are between GG and OA decoders. CBD decoder uses compact matrix representation and therefore has lower memory usage. BD decoder, on the other hand, employes pivoting techniques and therefore its decoding complexity is lower. BD decoder uses full-size matrix representation for random access, as is heavily needed during pivoting. So BD decoder has higher memory usage. Note that since the RLNC can be viewed as a special band code with band-length being the total number of source packets, BD and CBD decoder can be used to decode RLNC as well, which we referred to as the *naive* mode. A PP(perpetual) decoder is also provided, which can only be used to decode the window-wrapped codes. It handles wrap-around encoding vectors more carefully than the CBD decoder, which would treat encoding vectors of the window-wrapped codes naively as dense vectors. PP decoder is merely included for comparisons in research.

Systematic coding is also supported, which generates coded packets from each subgeneration only after each source packet therein is sent once. The decoding cost can be significantly reduced when the code is used in networks with low packet loss rate. The price to pay is higher overhead if the number of subgenerations is greater than 1. It is recommended to use systematic coding for RLNC in many scenarios (e.g., [11]).

For more details about subgeneration-based codes and the decoder design, please refer to reference papers below.

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
[1] Ahlswede, Rudolf; N. Cai, Shuo-Yen Robert Li, and Raymond Wai-Ho Yeung. "Network Information Flow". IEEE Transactions on Information Theory, 46 (4): 1204–1216, 2000.

[2] T. Ho, M. Medard, R. Koetter, and D. R. Karger, "A Random Linear Network Coding Approach to Multicast", in IEEE Transactions on Information Theory, Vol 52, No. 10, pp. 4413–4430, 2006.

[3] Ye Li, W.-Y. Chan, and S. D. Blostein, "Network Coding with Unequal Size Overlapping Generations", in Proceedings of International Symposium on Network Coding (NetCod), pp. 161-166, Cambridge, MA, June, 2012.

[4] Yao Li, E. Soljanin and P. Spasojevic, "Effects of the Generation Size and Overlap on Throughput and Complexity in Randomized Linear Network Coding," in IEEE Transactions on Information Theory, vol. 57, no. 2, pp. 1111-1123, Feb. 2011.

[5] A. Fiandrotti, V. Bioglio, M. Grangetto, R. Gaeta and E. Magli, "Band Codes for Energy-Efficient Network Coding With Application to P2P Mobile Streaming," in IEEE Transactions on Multimedia, vol. 16, no. 2, pp. 521-532, Feb. 2014.

[6] Ye Li, J. Zhu and Z. Bao, "Sparse Random Linear Network Coding With Precoded Band Codes," in IEEE Communications Letters, vol. 21, no. 3, pp. 480-483, March 2017.

[7] J. Heide, M. V. Pedersen, F. H. P. Fitzek and M. Medard, "A Perpetual Code for Network Coding," 2014 IEEE 79th Vehicular Technology Conference (VTC Spring), Seoul, 2014, pp. 1-6.

[8] Ye Li, S. Zhang, J. Wang, X. Ji, H. Wu and Z. Bao, "A Low-Complexity Coded Transmission Scheme over Finite-Buffer Relay Links," in IEEE Transactions on Communications, 2018 [Early Access].

[9] S. Yang and R. W. Yeung, "Batched Sparse Codes," in IEEE Transactions on Information Theory, vol. 60, no. 9, pp. 5322-5346, Sept. 2014.

[10] Ye Li, W. Y. Chan and S. D. Blostein, "On Design and Efficient Decoding of Sparse Random Linear Network Codes," in IEEE Access, vol. 5, pp. 17031-17044, 2017.

[11] Ye Li, W.-Y. Chan, and S. D. Blostein, "Systematic Network Coding for Two-Hop Lossy Transmissions", EURASIP Journal on Advances in Signal Processing, pp.1-14, 2015:93 DOI: 10.1186/s13634-015-0273-3 
