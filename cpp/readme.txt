A C++ port of the Cython version of idX.

The purpose of this port was to experiment with ways of reducing the memory used during the identification process, using C++ containiners. RapidJSON was selected to provide access to the kernel JSON Lines files. parallel_hashmap was selected to provide cross-platform-equivalent performance for large containers.

Requires the RapidJSON header files in a folder called rapidjson:

Requires the parallel_hashmap files in a folder called parallel_hashmap:

https://github.com/greg7mdp/parallel-hashmap/tree/master/parallel_hashmap

This port was written to be compliant with C++ 14 and it has been compiled using GCC version Ubuntu 7.4.0 and Visual Studio C++ 2019.
