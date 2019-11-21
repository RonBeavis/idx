A C++ port of the Cython version of idX

This port is to experiment with ways of reducing the memory and other system resources used during the identification process. The goal is to reduce the amount of energy and computer resources required to do MS/MS PSM assignment.

RapidJSON was selected to provide access to the kernel JSON Lines files.

parallel_hashmap was selected to provide cross-platform-equivalent performance for large containers.

Requires the RapidJSON header files in a folder called rapidjson:

   https://github.com/Tencent/rapidjson/tree/master/include/rapidjson

Requires the parallel_hashmap files in a folder called parallel_hashmap:

   https://github.com/greg7mdp/parallel-hashmap/tree/master/parallel_hashmap

This project was written to be compliant with C++ 14.

It has been compiled on the following platforms:

   Ubuntu 18, GCC version Ubuntu 7.4.0; and 
   Windows 10, Visual Studio C++ 2019 (using both the Microsoft C++ compiler and the Clang compiler). </li>

The current target release version is "2020.1".
