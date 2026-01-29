//===----------------------------------------------------------------------===//
//
//                         BusTub
//
// hyperloglog.cpp
//
// Identification: src/primer/hyperloglog.cpp
//
// Copyright (c) 2015-2025, Carnegie Mellon University Database Group
//
//===----------------------------------------------------------------------===//

#include "primer/hyperloglog.h"

namespace bustub {

/** @brief Parameterized constructor. */
template <typename KeyType>
HyperLogLog<KeyType>::HyperLogLog(int16_t n_bits) : cardinality_(0) {
  this->b = std::max(static_cast<int16_t>(0), n_bits); // Ensure non-negative
  size_t m = 1ULL << this->b; 
  registers_.assign(m, 0);
}

/**
 * @brief Function that computes binary.
 *
 * @param[in] hash
 * @returns binary of a given hash
 */
template <typename KeyType>
auto HyperLogLog<KeyType>::ComputeBinary(const hash_t &hash) const -> std::bitset<BITSET_CAPACITY> {
  return std::bitset<BITSET_CAPACITY>(hash);
}

/**
 * @brief Function that computes leading zeros.
 *
 * @param[in] bset - binary values of a given bitset
 * @returns leading zeros of given binary set
 */
template <typename KeyType>
auto HyperLogLog<KeyType>::PositionOfLeftmostOne(const std::bitset<BITSET_CAPACITY> &bset) const -> uint64_t {
  for(int i=BITSET_CAPACITY-1-this->b;i>=0;i--){
      if(bset[i]){

        return (BITSET_CAPACITY-this->b)-i;

      }
  }

  return 0;
  // by convention if all bits are zero, return 64
}

/**
 * @brief Adds a value into the HyperLogLog.
 *
 * @param[in] val - value that's added into hyperloglog
 */
template <typename KeyType>
auto HyperLogLog<KeyType>::AddElem(KeyType val) -> void {
  /** @TODO(student) Implement this function! */
   std :: bitset<BITSET_CAPACITY> bset = ComputeBinary(CalculateHash(val));
  
   uint64_t register_index = 0;
    for (int i = 0; i < this->b; ++i) {
        if (bset[BITSET_CAPACITY - 1 - i]) {
            register_index |= (1ULL << (this->b - 1 - i));
        }
    }
   uint64_t p=PositionOfLeftmostOne(bset);
   
   auto &registers= this->registers_;
   registers[register_index]=std::max(registers[register_index],static_cast<uint8_t>(p));
   
}

/**
 * @brief Function that computes cardinality.
 */
template <typename KeyType>
auto HyperLogLog<KeyType>::ComputeCardinality() -> void {
  /** @TODO(student) Implement this function! */
  int m=this->registers_.size();
  double sum_of_powers=0;
  for(int i=0;i<m;i++){
    sum_of_powers+=pow(2,-1*registers_[i]);
  }

  this->cardinality_=floor(CONSTANT * m * m / sum_of_powers);

}

template class HyperLogLog<int64_t>;
template class HyperLogLog<std::string>;

}  // namespace bustub
