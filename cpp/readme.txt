<h1>A C++ port of the Cython version of idX</h1>

<p>This port is to experiment with ways of reducing the memory and other system resources used during the identification process. The goal is to reduce the amount of energy and computer resources required to do MS/MS PSM assignment.</p>

<p>RapidJSON was selected to provide access to the kernel JSON Lines files. </p>

<p>parallel_hashmap was selected to provide cross-platform-equivalent performance for large containers.</p>

<p>Requires the RapidJSON header files in a folder called rapidjson:</p>
<ol>
   <li>https://github.com/Tencent/rapidjson/tree/master/include/rapidjson</li>
</ol>
<p>Requires the parallel_hashmap files in a folder called parallel_hashmap:</p>
<ol>
   <li>https://github.com/greg7mdp/parallel-hashmap/tree/master/parallel_hashmap</li>
</ol>
<p>This project was written to be compliant with C++ 14.</p>

<p>It has been compiled on the following platforms:</p>
<ol>
   <li>Ubuntu 18, GCC version Ubuntu 7.4.0; and 
   <li>Windows 10, Visual Studio C++ 2019 (using both the Microsoft C++ compiler and the Clang compiler). </li>
</ol>
<p>The current target release version is "2020.1".</p>
