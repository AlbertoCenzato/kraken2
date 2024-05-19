#ifndef KRAKEN2_MERGE_SORT_FILES_H_
#define KRAKEN2_MERGE_SORT_FILES_H_

#include <filesystem>

#pragma pack(4)
struct IdxValue{
  uint64_t idx;
  uint32_t value;
};

// order by idx and then by value
bool operator<(const IdxValue &lhs, const IdxValue &rhs);
bool operator==(const IdxValue &lhs, const IdxValue &rhs);

void multiStepMerge(
	const std::filesystem::path& hash_dir, 
	const std::filesystem::path& merge_dir, 
	const std::filesystem::path& output_file);


#endif