//===----------------------------------------------------------------------===//
//
//                         BusTub
//
// hyperloglog_presto.cpp
//
// Identification: src/primer/hyperloglog_presto.cpp
//
// Copyright (c) 2015-2025, Carnegie Mellon University Database Group
//
//===----------------------------------------------------------------------===//

#include "primer/hyperloglog_presto.h"

namespace bustub {

/** @brief Parameterized constructor. */
template <typename KeyType>
HyperLogLogPresto<KeyType>::HyperLogLogPresto(int16_t n_leading_bits) 
    : cardinality_(0), b(n_leading_bits) {
  size_t m = 1ULL << b;
  dense_bucket_.resize(m, std::bitset<DENSE_BUCKET_SIZE>(0));
}

/** @brief Function to count consecutive zeros in a bitset starting from a given index. */
template <typename KeyType>
auto HyperLogLogPresto<KeyType>::CountConsecutiveZeros(const std::bitset<BITSET_CAPACITY> &bset, int b_bits) -> uint64_t {
    uint64_t count = 0;
    // Start at the bit immediately after the index bits
    for (int i = BITSET_CAPACITY - 1 - b_bits; i >= 0; i--) {
        if (!bset[i]) {
            count++;
        } else {
            break; // Stop at the first '1'
        }
    }
    return count;
}

/** @brief Element is added for HLL calculation. */
template <typename KeyType>
auto HyperLogLogPresto<KeyType>::AddElem(KeyType val) -> void {
  // Use the provided CalculateHash
  hash_t hash = CalculateHash(val);
  std::bitset<BITSET_CAPACITY> bset(hash);

  // 1. Calculate Index (Leftmost b bits)
  uint64_t register_index = 0;
  for (int i = 0; i < this->b; ++i) {
    if (bset[BITSET_CAPACITY - 1 - i]) {
      register_index |= (1ULL << (this->b - 1 - i));
    }
  }

  // 2. Count Leading Zeros in Suffix
  uint64_t rho = CountConsecutiveZeros(bset, this->b);

  // 3. Update Buckets
  // If the total leading zeros (rho) is greater than what we have stored:
  uint64_t current_dense = dense_bucket_[register_index].to_ulong();
  uint64_t current_overflow = 0;
  if (overflow_bucket_.count(register_index) > 0) {
      current_overflow = overflow_bucket_[register_index].to_ulong();
  }
  
  uint64_t total_current = current_dense + current_overflow;

  if (rho > total_current) {
    if (rho <= 15) {
      dense_bucket_[register_index] = std::bitset<DENSE_BUCKET_SIZE>(rho);
      // If it was in overflow before but now fits in dense (unlikely in HLL but safe), 
      // you'd clear overflow. But HLL only increases.
    } else {
      dense_bucket_[register_index] = std::bitset<DENSE_BUCKET_SIZE>(15);
      overflow_bucket_[register_index] = std::bitset<OVERFLOW_BUCKET_SIZE>(rho - 15);
    }
  }
}

/** @brief Function to compute cardinality. */
template <typename T>
auto HyperLogLogPresto<T>::ComputeCardinality() -> void {
  /** @TODO(student) Implement this function! */
  // 1. m is the number of registers (2^b)
  double m = static_cast<double>(dense_bucket_.size());
  double harmonic_sum = 0.0;

  // 2. Iterate through every register to calculate the summation
  for (size_t i = 0; i < dense_bucket_.size(); ++i) {
    // Start with the 4-bit dense value
    uint64_t total_rho = dense_bucket_[i].to_ulong();
    
    // If an overflow exists for this index, add the 3-bit overflow value
    auto it = overflow_bucket_.find(static_cast<uint16_t>(i));
    if (it != overflow_bucket_.end()) {
      total_rho += it->second.to_ulong();
    }
    
    // Calculate 2^-Rj and add to the harmonic sum
    harmonic_sum += std::pow(2.0, -static_cast<double>(total_rho));
  }

  // 3. Apply the final estimate formula
  // CONSTANT (alpha_m) is defined as 0.79402 in your header
  double estimate = CONSTANT * m * m / harmonic_sum;
  
  // 4. Update the member variable
  cardinality_ = static_cast<uint64_t>(estimate);
}

template class HyperLogLogPresto<int64_t>;
template class HyperLogLogPresto<std::string>;
}  // namespace bustub
