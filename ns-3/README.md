Usage
============
As an example of using sparsenc library in ns-3, this folder provides a two-hop lossy relay link simulation script which uses network coding on top of UDP to send a number of source packets. The simulation relies on static sparsenc library libsparsenc.a, which can be obtained by 

```shell
$ make libsparsenc.a
```

Then, copy `../include/sparsenc.h`, `../libsparsenc.a`, `sparsenc-udp.cc`, `lossy-two-hop-tcp.cc`, and `wscript` to an ns-3 directory, say `$NS_FOLDER/ns-3.29/example/sparsenc/`. Edit `wscript` accordingly to reflect the actual relative path of the files. In the main folder of ns-3, i.e., `$NS_FOLDER/ns-3.29/`, run

```shell
$ ./waf --run sparsenc-udp
```

for simulations. You can customize `sparsenc-udp.cc` for other simulation parameters. Note that the `sparsenc-udp.cc` uses BATS_SNC code for the transmission, but other codes implemented in `libsparsenc.a` can also be used. You can follow C simulation scripts in `../examples/` and modify the ns-3 script accordingly. `lossy-two-hop-tcp.cc` is simply a baseline scheme which uses TCP for the lossy two-hop transmissions.
