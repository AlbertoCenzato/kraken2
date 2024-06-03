#ifndef KRAKEN2_MERGE_SORT_FILES_H_
#define KRAKEN2_MERGE_SORT_FILES_H_

#include "utils/streams.h"
#include <filesystem>
#include <vector>

#pragma pack(4)
struct IdxValue{
  uint64_t idx;
  uint32_t value;
};

// order by idx and then by value
bool operator<(const IdxValue &lhs, const IdxValue &rhs);
bool operator==(const IdxValue &lhs, const IdxValue &rhs);

void mergeSortStreams(
		std::vector<std::unique_ptr<IStreamReader>> input_streams,
    IStreamWriter &output_stream);

void multiStepMerge(
	const std::filesystem::path& hash_dir, 
	const std::filesystem::path& merge_dir, 
	const std::filesystem::path& output_file);

void makeHashTable(
	const std::filesystem::path& merged_hash_file, 
	const std::filesystem::path& output_file,
	size_t capacity, size_t key_bits, size_t value_bits);

#endif