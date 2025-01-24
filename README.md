# Engorgio: An Arbitrary-Precision Unbounded-Size Hybrid  Encrypted Database via Quantized Fully Homomorphic Encryption
This code is the implementation of the paper Engorgio: An Arbitrary-Precision Unbounded-Size Hybrid  Encrypted Database via Quantized Fully Homomorphic Encryption.

## Requirements

```
git 
gcc >= 10
cmake >= 3.16
GMP 6.2.0 
```

## Building Engorgio
You can build the Engorgio (out-of-source) for your machine by executing the following commands:
```
mkdir build
cd build
cmake ..
make
```

Following the instructions above, output binaries will be generated in the `build/bin/` directory. You can run these binaries by:
```
$./bin/comparison_test
$./bin/sort_test
$./bin/topk_test
$./bin/vectorized_query_test
$./bin/relational_query_test
$./bin/hybrid_query_test
```

## Examples
### Homomorphic Comparison
- codes `test/comparison_test.cpp`  
- output binary `build/bin/comparison_test` 
- This demo shows the homomorphic comparison operator `ArbQuantComp` in Engorgio. 

### Homomorphic Sorting
- codes `test/sort_test.cpp`  
- output binary `build/bin/sort_test` 
- This demo shows the homomorphic sorting operator `HomSort` in Engorgio.

### Homomorphic Top-k
- codes `test/topk_test.cpp`  
- output binary `build/bin/topk_test` 
- This demo shows the homomorphic Topk-k operator `HomTopK` in Engorgio.

### Relational Query Evaluation
- codes `test/relational_query_test.cpp`  
- output binary `build/bin/relational_query_test` 
- This demo evaluates relational queries on an encrypted database with varying numbers of rows.

### Vectorized Query Evaluation
- codes `test/vectorized_query_test.cpp`
- output binary `build/bin/vectorized_query_test`
- This demo evaluates vectorized queries on an encrypted database with varying numbers of rows.

### Hybrid Query Evaluation
- codes `test/hybrid_query_test.cpp`
- output binary `build/bin/hybrid_query_test`
- This demo evaluates hybrid queries on an encrypted database with varying numbers of rows.
